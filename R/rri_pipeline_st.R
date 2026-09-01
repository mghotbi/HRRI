#' Exploratory domain-score integration (legacy interface)
#' @description Integrates plant, soil and microbial latent scores. Direction is
#' not biologically identifiable without justified anchors. All three domains
#' may be incomplete. Available positive domain weights are renormalized per row.
#' Use rri_reference_scores for externally anchored, fixed-reference comparisons.
#' @param ROS_flux Data frame of plant physiological variables (rows = samples).
#' @param Eh_stability Data frame of soil redox chemistry variables (rows = samples).
#' @param micro_data Optional data frame of microbial abundance or functional features.
#' @param graph Optional \code{igraph} object or list of \code{igraph} objects
#'   representing microbial network structure.
#' @param id Optional data frame describing experimental design (same number of rows as inputs).
#' @param time_col Optional character. Name of time column in \code{id}.
#' @param group_cols Optional character vector of grouping variables in \code{id}.
#' @param mode Character. One of \code{"snapshot"}, \code{"rolling"}, or \code{"event"}.
#' @param window Integer >= 2. Rolling window size (for mode = "rolling").
#' @param align Character. Alignment rule for rolling window:
#'   \code{"right"}, \code{"center"}, or \code{"left"}.
#' @param event_col Optional character. Column in \code{id} identifying event phases.
#' @param baseline_label Character. Label identifying baseline phase.
#' @param recovery_labels Character vector identifying recovery phases.
#' @param alpha_micro Numeric between 0 and 1 controlling blending of microbial
#'   abundance and network components.
#' @param method_phys Character. Reduction method for plant block.
#' @param method_soil Character. Reduction method for soil block.
#' @param method_micro Character. Reduction method for microbial block.
#' @param direction_phys Character. Orientation rule for plant latent dimension.
#' @param direction_soil Character. Orientation rule for soil latent dimension.
#' @param direction_micro Character. Orientation rule for microbial latent dimension.
#' @param direction_anchor_phys Optional character. Anchor variable for plant orientation.
#' @param direction_anchor_soil Optional character. Anchor variable for soil orientation.
#' @param direction_anchor_micro Optional character. Anchor variable for microbial orientation.
#' @param scale_by Optional character vector of grouping variables used for scaling.
#' @param network_agg Character. Network aggregation method: \code{"equation"} or \code{"mean"}.
#' @param w1 Numeric weight for plant domain.
#' @param w2 Numeric weight for soil domain.
#' @param w3 Numeric weight for microbial domain. Must sum with w1 and w2 to 1.
#' @param add_coupling Logical. If TRUE, adds cross-domain coherence term.
#' @param coupling_weight Numeric between 0 and 1 controlling weight of coupling term.
#' @param coupling_fun Character. Coupling function: \code{"geometric_mean"} or \code{"agreement"}.
#' @param norm_method Optional character. If provided, overrides block-specific methods.
#' @param reducer Character. Reduction strategy: \code{"per_domain"} or \code{"mfa"}.
#' @param scaling Character. Scaling rule: \code{"minmax_legacy"} or \code{"pnorm"}.
#' @param comp_space Character. Compositional projection method:
#'   \code{"closure_legacy"} or \code{"clr"}.
#' @param ref_stats Optional list of reference statistics used for scaling.
#' @param add_compensation Logical. If TRUE, includes covariance-based compensation term.
#' @param compensation_weight Numeric between 0 and 1 controlling compensation weight.
#'
#' @details MFA is disabled pending a validated implementation. Scaling statistics
#' do not freeze PCA/FA loadings, so ref_stats is not a trained prediction model.
#' The CLR round trip changes display coordinates only: inversion returns closure.
#' Grouping does not imply within-group scaling; request scale_by explicitly.
#' Missing data are median-imputed for exploratory reduction, not corrected for MNAR.
#' Event scores are descriptive products of resistance and reference proximity;
#' baseline/recovery label defaults must be matched to the supplied data.
#' @return RRI object; identifiers accompany scores and rolling output retains
#' original input order. Stochastic and advanced reducers need separate validation.
#' @importFrom stats median sd prcomp pnorm cov
#' @importFrom igraph simplify vcount ecount edge_attr is_directed cluster_fast_greedy modularity transitivity global_efficiency centr_degree
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline_st(sim$ROS_flux, sim$Eh_stability, id = sim$id)
#'   head(res$row_scores)
#' }
#' @export
rri_pipeline_st <- function(
  ROS_flux = NULL,
  Eh_stability = NULL,
  micro_data = NULL,
  graph = NULL,
  id = NULL,
  time_col = NULL,
  group_cols = NULL,
  mode = c("snapshot", "rolling", "event"),
  window = 3,
  align = c("right", "center", "left"),
  event_col = NULL,
  baseline_label = "pre",
  recovery_labels = "recovery",
  alpha_micro = 0.5,
  method_phys = "pca",
  method_soil = "pca",
  method_micro = "pca",
  direction_phys = c("auto", "higher_is_better", "lower_is_better"),
  direction_soil = c("auto", "higher_is_better", "lower_is_better"),
  direction_micro = c("auto", "higher_is_better", "lower_is_better"),
  direction_anchor_phys = NULL,
  direction_anchor_soil = NULL,
  direction_anchor_micro = NULL,
  scale_by = NULL,
  network_agg = c("equation", "mean"),
  w1 = 0.4,
  w2 = 0.35,
  w3 = 0.25,
  add_coupling = FALSE,
  coupling_weight = 0,
  coupling_fun = c("geometric_mean", "agreement"),
  norm_method = NULL,
  reducer = c("per_domain", "mfa"),
  scaling = c("minmax_legacy", "pnorm"),
  comp_space = c("closure_legacy", "clr"),
  ref_stats = NULL,
  add_compensation = FALSE,
  compensation_weight = 0
) {
  mode <- match.arg(mode)
  align <- match.arg(align)
  network_agg <- match.arg(network_agg)
  coupling_fun <- match.arg(coupling_fun)
  direction_phys <- match.arg(direction_phys)
  direction_soil <- match.arg(direction_soil)
  direction_micro <- match.arg(direction_micro)
  reducer <- match.arg(reducer)
  scaling <- match.arg(scaling)
  comp_space <- match.arg(comp_space)

  if (!is.null(norm_method)) {
    method_phys <- norm_method
    method_soil <- norm_method
    method_micro <- norm_method
  }

  if (isTRUE(add_compensation) && isTRUE(add_coupling))
    stop("Combining compensation and coupling is unsupported; legacy coupling overwrote compensation.")
  if (any(c(method_phys, method_soil, method_micro) %in% c("umap", "nmf", "wgcna")))
    stop("UMAP, NMF and WGCNA routes are disabled in this audit patch pending reproducibility and row-alignment validation.")
  if (is.null(direction_anchor_phys) || is.null(direction_anchor_soil) || is.null(direction_anchor_micro))
    warning("Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.")
  if (reducer == "mfa") stop("MFA is disabled pending a validated implementation.")
  if (any(!c(method_phys, method_soil, method_micro) %in% c("mean", "scale", "pca", "fa")))
    stop("Supported reducers are mean, scale, pca and fa.")
  # ---- validations --------------------------------------------------------
  # Assisted-by: OpenAI Codex. Do not use hidden simulator states as observations.
  ROS_flux <- .rri_observed_block(ROS_flux)
  Eh_stability <- .rri_observed_block(Eh_stability)
  micro_data <- .rri_observed_block(micro_data)
  blocks <- Filter(Negate(is.null), list(ROS_flux, Eh_stability, micro_data))
  if (!length(blocks)) stop("At least one numeric observation block is required.")
  ns <- vapply(blocks, nrow, integer(1))
  if (length(unique(ns)) != 1L || ns[1L] < 1L)
    stop("Observation blocks must have the same positive number of rows.")
  n <- ns[1L]
  observed <- lapply(list(ROS_flux, Eh_stability, micro_data), function(x)
    if (is.null(x)) rep(FALSE, n) else rowSums(is.finite(as.matrix(x))) > 0L)
  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (nrow(id) != n || anyDuplicated(names(id)))
      stop("id must have one row per observation and uniquely named columns.")
  }
  if (any(lengths(list(w1, w2, w3)) != 1L)) stop("Domain weights must be scalars.")
  if (any(!is.finite(c(w1, w2, w3))) || any(c(w1, w2, w3) < 0)) {
    stop("`w1`, `w2`, `w3` must be finite and >= 0.")
  }
  if (abs(w1 + w2 + w3 - 1) > 1e-8) stop("`w1 + w2 + w3` must sum to 1.")

  if (!is.finite(alpha_micro) || alpha_micro < 0 || alpha_micro > 1) {
    stop("`alpha_micro` must be in [0, 1].")
  }
  if (!is.finite(coupling_weight) || coupling_weight < 0 || coupling_weight > 1) {
    stop("`coupling_weight` must be in [0, 1].")
  }
  if (!is.finite(compensation_weight) || compensation_weight < 0 || compensation_weight > 1) {
    stop("`compensation_weight` must be in [0, 1].")
  }

  if (mode != "snapshot") {
    if (is.null(id) || is.null(time_col) || !time_col %in% names(id)) {
      stop("For mode != 'snapshot', provide `id` with a valid `time_col`.")
    }
    if (is.null(group_cols) || any(!group_cols %in% names(id))) {
      stop("For mode != 'snapshot', provide valid `group_cols` in `id`.")
    }
    if (!is.numeric(window) || length(window) != 1 || !is.finite(window) || window < 2 || window != floor(window)) {
      stop("`window` must be a single integer >= 2.")
    }
  }

  # Grouping trajectories must not silently change the normalization population.
  if (!is.null(scale_by) && (is.null(id) || any(!scale_by %in% names(id)))) {
    stop("If `scale_by` is provided, `id` must be provided and contain all `scale_by` columns.")
  }

  # ---- helpers ------------------------------------------------------------
  as_numeric_df <- function(x) .rri_numeric_df(x)

  impute_median_by <- function(df, id, by) {
    if (is.null(df)) return(NULL)
    df <- as_numeric_df(df)

    if (is.null(by) || is.null(id) || length(by) == 0) {
      for (j in seq_len(ncol(df))) {
        v <- df[[j]]
        v[!is.finite(v)] <- NA_real_
        if (all(is.na(v))) next
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]] <- v
      }
      return(df)
    }

    key <- .rri_key(id, by)
    for (k in unique(key)) {
      idx <- which(key == k)
      for (j in seq_len(ncol(df))) {
        v <- df[[j]][idx]
        v[!is.finite(v)] <- NA_real_
        if (all(is.na(v))) next
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]][idx] <- v
      }
    }

    for (j in seq_len(ncol(df))) {
      v <- df[[j]]
      if (anyNA(v) && !all(is.na(v))) {
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]] <- v
      }
    }

    df
  }

  scale_01 <- function(v) .rri_scale(v)

  scale_pnorm <- function(v, stats = NULL) {
    if (all(is.na(v))) {
      return(rep(NA_real_, length(v)))
    }
    mu <- if (is.list(stats) && !is.null(stats$mean)) stats$mean else mean(v, na.rm = TRUE)
    si <- if (is.list(stats) && !is.null(stats$sd)) stats$sd else stats::sd(v, na.rm = TRUE)
    if (!is.finite(si) || si == 0) si <- 1
    stats::pnorm((v - mu) / si)
  }

  scale_vec <- function(v, stats = NULL) {
    if (scaling == "minmax_legacy") {
      scale_01(v)
    } else {
      scale_pnorm(v, stats)
    }
  }

  scale_vec_by <- function(v, id, by, stats = NULL) {
    if (is.null(by) || is.null(id) || length(by) == 0) {
      return(scale_vec(v, stats))
    }
    key <- .rri_key(id, by)
    out <- v
    for (k in unique(key)) {
      idx <- which(key == k)
      out[idx] <- scale_vec(v[idx], stats)
    }
    out
  }

  orient_latent <- function(latent, df, direction, anchor) {
    if (is.null(df) || all(!is.finite(latent))) return(latent)
    # PCA/FA axis signs are arbitrary. Anchor before applying a desired direction.
    if (is.null(anchor)) {
      if (direction == "lower_is_better") return(-latent)
      return(latent)
    }
    if (!anchor %in% names(df)) stop("Orientation anchor is absent: ", anchor)
    cc <- .rri_cor(latent, df[[anchor]])
    if (is.finite(cc) && cc < 0) latent <- -latent
    if (direction == "lower_is_better") -latent else latent
  }

  latent_dimension <- function(df, method) {
    if (is.null(df)) return(rep(NA_real_, n))
    df <- as_numeric_df(df)
    varying <- vapply(df, function(x) {
      z <- x[is.finite(x)]
      length(z) > 1L && stats::sd(z) > 0
    }, logical(1))
    # All-NA columns are unmeasured; constants carry no latent variation.
    df <- df[, varying, drop = FALSE]
    if (ncol(df) > 1L && any(!is.finite(as.matrix(df))))
      stop("Non-finite values remain after imputation.")

    if (ncol(df) == 0L) return(rep(NA_real_, nrow(df)))
    if (ncol(df) == 1L) {
      return(df[[1L]])
    }

    if (method == "mean") {
      return(rowMeans(df, na.rm = TRUE))
    }
    if (method == "scale") {
      return(rowMeans(scale(df), na.rm = TRUE))
    }
    if (method == "pca") {
      return(stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1])
    }

    if (method == "fa") {
      if (!requireNamespace("psych", quietly = TRUE)) stop("`method = 'fa'` requires {psych}.")
      fa <- psych::fa(df, nfactors = 1, rotate = "none", scores = "regression")
      return(as.numeric(fa$scores[, 1]))
    }

    stop("Unknown method: ", method)
  }

  network_scalar <- function(g, agg) {
    if (!requireNamespace("igraph", quietly = TRUE)) stop("`graph` requires {igraph}.")
    if (!inherits(g, "igraph")) stop("Each graph must be an igraph object.")
    g <- igraph::simplify(g)
    if (igraph::vcount(g) < 2L || igraph::ecount(g) == 0L) return(NA_real_)
    gw <- igraph::edge_attr(g, "weight")
    if (!is.null(gw) && any(!is.finite(gw))) stop("Graph weights must be finite.")
    if (igraph::is_directed(g) || any(igraph::edge_attr(g, "weight") < 0))
      stop("Network summary requires an undirected graph without negative weights.")

    # Deterministic partition avoids modifying the caller RNG state.
    comm <- igraph::cluster_fast_greedy(g, weights = NA)
    Q <- igraph::modularity(comm)

    C_vec <- igraph::transitivity(g, type = "local", isolates = "zero")
    C <- mean(C_vec, na.rm = TRUE)

    Eglob <- igraph::global_efficiency(g, directed = FALSE, weights = NA)
    H <- igraph::centr_degree(g)$centralization
    Qs <- pmin(pmax(1 - Q, 0), 1)
    Cs <- pmin(pmax(C, 0), 1)
    Hs <- pmin(pmax(H, 0), 1)
    Es <- Eglob

    metrics <- c(Qs, Cs, Es, (1 - Hs))
    if (agg == "mean") {
      return(mean(metrics, na.rm = TRUE))
    }
    sqrt(mean(metrics, na.rm = TRUE))
  }

  roll_apply <- function(x, window, fun, align) {
    nloc <- length(x)
    out <- rep(NA_real_, nloc)

    for (i in seq_len(nloc)) {
      if (align == "right") {
        idx <- (i - window + 1):i
      } else if (align == "left") {
        idx <- i:(i + window - 1)
      } else {
        half <- floor(window / 2)
        idx <- (i - half):(i - half + window - 1)
      }

      idx <- idx[idx >= 1 & idx <= nloc]
      if (length(idx) < window) next
      out[i] <- fun(x[idx])
    }

    out
  }

  event_features <- function(z, event, baseline_label, recovery_labels) {
    ok <- is.finite(z) & !is.na(event)
    if (sum(ok) < 3) {
      return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    }

    z <- z[ok]
    e <- event[ok]
    base_idx <- which(e == baseline_label)
    rec_idx <- which(e %in% recovery_labels)
    stress_idx <- which(!e %in% c(baseline_label, recovery_labels))

    if (length(base_idx) < 1) {
      return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    }

    baseline <- mean(z[base_idx])
    pulse <- if (length(stress_idx) >= 1) max(abs(z[stress_idx] - baseline)) else NA_real_
    recovery <- if (length(rec_idx) >= 1) 1 - abs(mean(z[rec_idx]) - baseline) else NA_real_

    list(baseline = baseline, pulse = pulse, recovery = recovery)
  }

  coupling_term <- function(P, S, M, fun) {
    if (fun == "geometric_mean") {
      (P * S * M)^(1 / 3)
    } else {
      1 - stats::sd(c(P, S, M))
    }
  }

  compensation_index <- function(P, S, M) {
    X <- cbind(P, S, M)
    if (nrow(X) < 3) {
      return(NA_real_)
    }
    C <- stats::cov(X, use = "pairwise.complete.obs")
    if (any(!is.finite(C))) {
      return(NA_real_)
    }
    -sum(C[upper.tri(C)])
  }

  clr_transform <- function(X, eps = 1e-8) {
    X <- pmax(as.matrix(X), eps)
    gm <- exp(rowMeans(log(X)))
    log(X / gm)
  }

  clr_to_simplex <- function(CLR) {
    X <- exp(CLR)
    X / rowSums(X)
  }

  # ---- impute (before latent extraction) ----------------------------------
  ROS_flux <- impute_median_by(ROS_flux, id, scale_by)
  Eh_stability <- impute_median_by(Eh_stability, id, scale_by)
  if (!is.null(micro_data)) micro_data <- impute_median_by(micro_data, id, scale_by)

  # ---- microbial inputs: abundance + network -------------------------------
  micro_abund_raw <- rep(NA_real_, n)
  if (!is.null(micro_data)) {
    micro_abund_raw <- orient_latent(
      latent_dimension(micro_data, method_micro),
      micro_data,
      direction_micro,
      direction_anchor_micro
    )
  }

  micro_net_raw <- rep(NA_real_, n)
  if (!is.null(graph)) {
    if (inherits(graph, "igraph")) {
      micro_net_raw <- rep(network_scalar(graph, network_agg), n)
    } else if (is.list(graph)) {
      if (length(graph) != n) stop("If `graph` is a list, it must have length `nrow(ROS_flux)`.")
      micro_net_raw <- vapply(graph, function(g) network_scalar(g, network_agg), numeric(1))
    } else {
      stop("`graph` must be NULL, an igraph object, or a list of igraph objects.")
    }
  }

  # ---- domain raw scores (per-domain default) ------------------------------
  phys_raw <- orient_latent(
    latent_dimension(ROS_flux, method_phys),
    ROS_flux,
    direction_phys,
    direction_anchor_phys
  )

  soil_raw <- orient_latent(
    latent_dimension(Eh_stability, method_soil),
    Eh_stability,
    direction_soil,
    direction_anchor_soil
  )

  # Only defined/used for MFA; keep as NA otherwise.
  micro_mfa_raw <- rep(NA_real_, n)

  # Retained legacy fields; no unsupported MFA projections are calculated.
  mfa_global_dim1 <- rep(NA_real_, n)
  mfa_used_partials <- FALSE
  mfa_fallback_to_per_domain <- FALSE
  phys_raw[!observed[[1L]]] <- NA_real_
  soil_raw[!observed[[2L]]] <- NA_real_
  micro_abund_raw[!observed[[3L]]] <- NA_real_

  # ---- scale to [0, 1] -----------------------------------------------------
  phys <- scale_vec_by(phys_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$phys else NULL)
  soil <- scale_vec_by(soil_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$soil else NULL)

  micro_abund <- scale_vec_by(
    micro_abund_raw,
    id,
    scale_by,
    stats = if (is.list(ref_stats)) ref_stats$micro_abund else NULL
  )
  micro_net <- scale_vec_by(
    micro_net_raw,
    id,
    scale_by,
    stats = if (is.list(ref_stats)) ref_stats$micro_net else NULL
  )
  micro_mfa <- scale_vec_by(
    micro_mfa_raw,
    id,
    scale_by,
    stats = if (is.list(ref_stats)) ref_stats$micro_mfa else NULL
  )

  # ---- blend microbial domain ---------------------------------------------
  micro <- .rri_weighted(cbind(micro_abund, micro_net),
                         c(alpha_micro, 1 - alpha_micro))$score
  domains <- cbind(Physio = phys, Soil = soil, Micro = micro)
  scored <- .rri_weighted(domains, c(w1, w2, w3))
  base_rri <- scored$score
  if ((isTRUE(add_compensation) && compensation_weight > 0) ||
      (isTRUE(add_coupling) && coupling_weight > 0)) {
    if (any(!is.finite(domains)))
      stop("Coupling and compensation require all three domains on every row.")
  }

  if (isTRUE(add_compensation) && compensation_weight > 0) {
    comp_term <- rep(NA_real_, n)

    if (!is.null(id) &&
      !is.null(group_cols) &&
      length(group_cols) > 0 &&
      all(group_cols %in% names(id))) {
      keyc <- .rri_key(id, group_cols)
      for (k in unique(keyc)) {
        idx <- which(keyc == k)
        comp_term[idx] <- compensation_index(phys[idx], soil[idx], micro[idx])
      }
    } else {
      comp_term[] <- compensation_index(phys, soil, micro)
    }

    comp_term01 <- scale_vec(
      comp_term,
      stats = if (is.list(ref_stats)) ref_stats$comp else NULL
    )
    if (any(!is.finite(comp_term01)))
      stop("Compensation is not estimable for every requested group.")

    ww <- c(w1, w2, w3)
    ww <- ww / sum(ww) * (1 - compensation_weight)

    base_rri <- ww[1] * phys + ww[2] * soil + ww[3] * micro + compensation_weight * comp_term01
  }

  if (isTRUE(add_coupling) && coupling_weight > 0) {
    ww <- c(w1, w2, w3)
    ww <- ww / sum(ww) * (1 - coupling_weight)

    coup <- mapply(
      coupling_term,
      P = phys,
      S = soil,
      M = micro,
      MoreArgs = list(fun = coupling_fun)
    )

    base_rri <- ww[1] * phys + ww[2] * soil + ww[3] * micro + coupling_weight * coup
  }

  rri <- base_rri # weighted [0,1] domain mean; do not force extrema after aggregation

  row_scores <- data.frame(
    Physio = phys,
    Soil = soil,
    Micro = micro,
    RRI = rri,
    domain_coverage = scored$coverage,
    n_domains = scored$n_observed,
    Micro_abundance = micro_abund,
    Micro_network = micro_net,
    Micro_mfa = micro_mfa
  )

  if (!is.null(id)) {
    if (any(names(id) %in% names(row_scores))) stop("id names conflict with score columns.")
    row_scores <- cbind(id, row_scores)
  }
  # ---- compositional projection (ternary-ready) ---------------------------
  comp <- row_scores[, c("Physio", "Soil", "Micro")]

  if (comp_space == "clr") {
    clr <- clr_transform(comp)
    comp_simplex <- clr_to_simplex(clr)
    comp <- as.data.frame(comp_simplex)
    attr(comp, "clr") <- clr
  } else {
    s <- rowSums(comp)
    s[s == 0] <- NA_real_
    comp <- comp / s
  }

  row_scores_comp <- cbind(comp, RRI = row_scores$RRI)

  rri_index <- .rri_mean(row_scores$RRI)
  attr(row_scores_comp, "RRI_index") <- rri_index

  # ---- dynamic modes -------------------------------------------------------
  dyn_scores <- NULL

  if (mode != "snapshot") {
    id2 <- id
    if (!is.numeric(id2[[time_col]]) || any(!is.finite(id2[[time_col]])))
      stop("Time must be finite numeric values in explicitly stated units.")
    if (anyDuplicated(.rri_key(id2, c(group_cols, time_col))))
      stop("Duplicate times within trajectory: include all experimental-unit keys.")

    key <- .rri_key(id2, group_cols)
    ord <- order(key, id2[[time_col]])

    rs_ord <- row_scores[ord, , drop = FALSE]
    id_ord <- id2[ord, , drop = FALSE]
    key_ord <- key[ord]

    if (mode == "rolling") {
      stab_fun <- function(x) 1 - stats::sd(x)

      P_roll <- S_roll <- M_roll <- rep(NA_real_, nrow(rs_ord))
      P_stab <- S_stab <- M_stab <- rep(NA_real_, nrow(rs_ord))

      for (k in unique(key_ord)) {
        idx <- which(key_ord == k)
        P_roll[idx] <- roll_apply(rs_ord$Physio[idx], window, mean, align)
        S_roll[idx] <- roll_apply(rs_ord$Soil[idx], window, mean, align)
        M_roll[idx] <- roll_apply(rs_ord$Micro[idx], window, mean, align)

        P_stab[idx] <- roll_apply(rs_ord$Physio[idx], window, stab_fun, align)
        S_stab[idx] <- roll_apply(rs_ord$Soil[idx], window, stab_fun, align)
        M_stab[idx] <- roll_apply(rs_ord$Micro[idx], window, stab_fun, align)
      }

      P_dyn <- scale_vec(0.5 * P_roll + 0.5 * P_stab)
      S_dyn <- scale_vec(0.5 * S_roll + 0.5 * S_stab)
      M_dyn <- scale_vec(0.5 * M_roll + 0.5 * M_stab)

      rri_dyn <- .rri_weighted(cbind(P_dyn, S_dyn, M_dyn), c(w1, w2, w3))$score

      dyn_scores <- data.frame(
        P_level = P_roll,
        P_stability = P_stab,
        S_level = S_roll,
        S_stability = S_stab,
        M_level = M_roll,
        M_stability = M_stab,
        Physio_dyn = P_dyn,
        Soil_dyn = S_dyn,
        Micro_dyn = M_dyn,
        RRI_dyn = rri_dyn
      )
      dyn_scores <- cbind(id_ord, dyn_scores)
      dyn_scores <- dyn_scores[order(ord), , drop = FALSE]
      rownames(dyn_scores) <- NULL
    } else if (mode == "event") {
      if (is.null(event_col) || !event_col %in% names(id_ord)) {
        stop("mode = 'event' requires `event_col` in `id`.")
      }
      if (any(baseline_label %in% recovery_labels))
        stop("Baseline and recovery labels must not overlap.")
      e <- id_ord[[event_col]]

      out_list <- vector("list", length(unique(key_ord)))
      names(out_list) <- unique(key_ord)

      for (k in unique(key_ord)) {
        idx <- which(key_ord == k)

        P_feat <- event_features(rs_ord$Physio[idx], e[idx], baseline_label, recovery_labels)
        S_feat <- event_features(rs_ord$Soil[idx], e[idx], baseline_label, recovery_labels)
        M_feat <- event_features(rs_ord$Micro[idx], e[idx], baseline_label, recovery_labels)

        P_res <- pmax(0, 1 - P_feat$pulse) * pmax(0, P_feat$recovery)
        S_res <- pmax(0, 1 - S_feat$pulse) * pmax(0, S_feat$recovery)
        M_res <- pmax(0, 1 - M_feat$pulse) * pmax(0, M_feat$recovery)

        rri_e <- .rri_weighted(matrix(c(P_res, S_res, M_res), nrow = 1L),
                               c(w1, w2, w3))$score
        gdat <- id_ord[idx[1], group_cols, drop = FALSE]

        out_list[[k]] <- cbind(
          gdat,
          P_baseline = P_feat$baseline,
          P_pulse = P_feat$pulse,
          P_recovery = P_feat$recovery,
          S_baseline = S_feat$baseline,
          S_pulse = S_feat$pulse,
          S_recovery = S_feat$recovery,
          M_baseline = M_feat$baseline,
          M_pulse = M_feat$pulse,
          M_recovery = M_feat$recovery,
          Physio_event = P_res,
          Soil_event = S_res,
          Micro_event = M_res,
          RRI_event = rri_e
        )
      }

      dyn_scores <- do.call(rbind, out_list)
      rownames(dyn_scores) <- NULL
    }
  }

  out <- list(
    row_scores = row_scores,
    effective_weights = sweep(ifelse(is.finite(domains), 1, 0), 2L,
                              c(w1, w2, w3), "*") / ifelse(scored$coverage > 0, scored$coverage, NA_real_),
    row_scores_comp = row_scores_comp,
    dyn_scores = dyn_scores,
    meta = list(
      mode = mode,
      graph = graph,
      domain_observed_rows = observed,
      missingness_policy = "Unmeasured domains and wholly missing rows remain NA; available positive weights are renormalized",
      interpretation = "Exploratory feature composite; reference scaling does not freeze PCA loadings",
      time_col = time_col,
      group_cols = group_cols,
      scale_by = scale_by,
      method = list(phys = method_phys, soil = method_soil, micro = method_micro),
      direction = list(phys = direction_phys, soil = direction_soil, micro = direction_micro),
      alpha_micro = alpha_micro,
      reducer = reducer,
      scaling = scaling,
      comp_space = comp_space,
      add_compensation = add_compensation,
      compensation_weight = compensation_weight,
      weights = list(
        w1 = w1,
        w2 = w2,
        w3 = w3,
        add_coupling = add_coupling,
        coupling_weight = coupling_weight,
        coupling_fun = coupling_fun
      ),
      rri_index = rri_index,
      mfa_global_dim1 = mfa_global_dim1,
      mfa_used_partials = mfa_used_partials,
      mfa_fallback_to_per_domain = mfa_fallback_to_per_domain
    )
  )

  class(out) <- "RRI"
  out
}
