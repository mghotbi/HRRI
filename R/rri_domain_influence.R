#' @title Realised influence of each domain on the composite score
#'
#' @description
#' A declared weight is not the same thing as realised influence. Two domains
#' given equal weight contribute unequally to the composite whenever their
#' scores differ in dispersion, in how strongly they covary with the other
#' domains, or in how often they are missing. This function reports what each
#' domain actually contributed, so that a claim such as "the index is driven by
#' soil" can be checked rather than inferred from the pattern of a figure.
#'
#' @param res An `RRI` object from [rri_pipeline()] or [rri_pipeline_st()].
#' @param domains Character vector of domain score columns. Defaults to
#'   `c("Physio", "Soil", "Micro")`.
#' @param rri_col Name of the composite column. Default `"RRI"`.
#' @param weights Optional named numeric vector of the nominal weights used to
#'   build the composite. If `NULL` (default) the function tries
#'   `res$effective_weights`, then `res$meta$weights`, and otherwise reports
#'   realised influence without a nominal comparison.
#'
#' @return A list with three elements.
#'   \describe{
#'     \item{`influence`}{One row per domain: `mean`, `sd`, `n_missing`,
#'       `cor_with_rri`, `nominal_weight`, `realised_share` and `ratio`.}
#'     \item{`covariance`}{Pairwise correlations between domain scores.
#'       Strong cross-domain correlation means influence cannot be attributed
#'       cleanly to one domain.}
#'     \item{`notes`}{Character vector of diagnostics worth acting on.}
#'   }
#'
#' @details
#' **How realised share is computed.** For weights \eqn{w_d} and domain scores
#' \eqn{S_d}, the composite is \eqn{R = \sum_d w_d S_d}. Because
#' \eqn{\mathrm{Var}(R) = \sum_d w_d \mathrm{Cov}(S_d, R)}, the quantity
#'
#' \deqn{\phi_d = w_d \, \mathrm{Cov}(S_d, R) / \mathrm{Var}(R)}
#'
#' is an exact decomposition: the \eqn{\phi_d} sum to one. Each domain's share
#' therefore includes its own variance *and* its share of the covariance it has
#' with the other domains. This is the appropriate attribution when domains are
#' correlated, which they generally are under a shared forcing.
#'
#' **How to read `ratio`.** `ratio` is realised share divided by nominal
#' weight. A value near 1 means the domain influenced the composite about as
#' much as intended. Values above roughly 1.3 or below roughly 0.7 indicate
#' that the declared weights are not delivering the intended balance, usually
#' for one of three reasons: the domain score is more (or less) dispersed than
#' the others after scaling; it covaries strongly with the others, so it
#' absorbs shared variance; or it is missing for many rows, so per-row weight
#' renormalisation quietly redistributes its weight.
#'
#' **What this does not establish.** A high realised share is a statement about
#' the score, not about the ecosystem. It does not show that the domain is
#' mechanistically more important, and it is not evidence that the composite is
#' wrong. It shows only where the variance in this particular composite came
#' from, for this cohort, under these weights and this scaling.
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2, n_depth = 2, n_plant = 3, n_time = 40,
#'   p_micro = 25, seed = 2026
#' )
#'
#' res <- suppressWarnings(rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   reducer = "per_domain",
#'   scaling = "pnorm"
#' ))
#'
#' infl <- rri_domain_influence(res)
#' infl$influence
#' infl$notes
#'
#' @seealso [rri_sensitivity()] for the effect of alternative weight grids;
#'   [rri_compensation_index()] for cross-domain asynchrony.
#' @importFrom stats cor var complete.cases
#' @export
rri_domain_influence <- function(res,
                                 domains = c("Physio", "Soil", "Micro"),
                                 rri_col = "RRI",
                                 weights = NULL) {

  rs <- if (is.list(res) && !is.data.frame(res) && !is.null(res$row_scores)) {
    as.data.frame(res$row_scores)
  } else {
    as.data.frame(res)
  }

  missing_cols <- setdiff(c(domains, rri_col), names(rs))
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  ## ---- nominal weights -----------------------------------------------------
  if (is.null(weights)) {
    w <- NULL
    if (is.list(res) && !is.null(res$effective_weights)) {
      ew <- res$effective_weights
      if (is.data.frame(ew) && nrow(ew)) {
        ew <- colMeans(ew[, intersect(domains, names(ew)), drop = FALSE],
                       na.rm = TRUE)
      }
      if (is.numeric(ew) && length(ew)) w <- ew
    }
    if (is.null(w) && is.list(res) && !is.null(res$meta$weights)) {
      mw <- unlist(res$meta$weights[c("w1", "w2", "w3")])
      if (length(mw) == 3L) { w <- mw; names(w) <- domains }
    }
    weights <- w
  }
  if (!is.null(weights)) {
    weights <- weights[intersect(domains, names(weights))]
    if (length(weights) != length(domains)) weights <- NULL
  }

  R <- suppressWarnings(as.numeric(rs[[rri_col]]))
  varR <- stats::var(R, na.rm = TRUE)

  ## ---- per-domain summaries -----------------------------------------------
  rows <- lapply(domains, function(d) {
    S <- suppressWarnings(as.numeric(rs[[d]]))
    ok <- is.finite(S) & is.finite(R)
    covSR <- if (sum(ok) > 2) stats::cov(S[ok], R[ok]) else NA_real_
    wd <- if (!is.null(weights)) unname(weights[[d]]) else NA_real_
    share <- if (is.finite(covSR) && is.finite(varR) && varR > 0 &&
                 is.finite(wd)) wd * covSR / varR else NA_real_
    data.frame(
      domain        = d,
      mean          = mean(S, na.rm = TRUE),
      sd            = stats::sd(S, na.rm = TRUE),
      n_missing     = sum(!is.finite(S)),
      cor_with_rri  = if (sum(ok) > 2) stats::cor(S[ok], R[ok]) else NA_real_,
      nominal_weight = wd,
      realised_share = share,
      stringsAsFactors = FALSE
    )
  })
  infl <- do.call(rbind, rows)
  infl$ratio <- infl$realised_share / infl$nominal_weight
  rownames(infl) <- NULL

  ## ---- cross-domain correlation -------------------------------------------
  D <- rs[, domains, drop = FALSE]
  D[] <- lapply(D, function(x) suppressWarnings(as.numeric(x)))
  cm <- suppressWarnings(stats::cor(D, use = "pairwise.complete.obs"))

  ## ---- notes ---------------------------------------------------------------
  notes <- character(0)

  if (is.null(weights)) {
    notes <- c(notes,
      "Nominal weights could not be recovered from `res`; supply `weights` to compare realised influence against intent.")
  } else {
    tot <- sum(infl$realised_share, na.rm = TRUE)
    if (is.finite(tot) && abs(tot - 1) > 0.02) {
      notes <- c(notes, sprintf(
        "Realised shares sum to %.3f rather than 1. This happens when the composite includes a coupling or compensation term, or when rows differ in domain coverage.", tot))
    }
    for (i in seq_len(nrow(infl))) {
      r <- infl$ratio[i]
      if (is.finite(r) && (r > 1.3 || r < 0.7)) {
        notes <- c(notes, sprintf(
          "%s contributed %.0f%% of composite variance against a nominal weight of %.0f%% (ratio %.2f).",
          infl$domain[i], 100 * infl$realised_share[i],
          100 * infl$nominal_weight[i], r))
      }
    }
  }

  miss <- infl$domain[infl$n_missing > 0]
  if (length(miss)) {
    notes <- c(notes, sprintf(
      "Missing domain scores in: %s. Per-row weight renormalisation redistributes their weight to the observed domains, so realised influence will depart from nominal.",
      paste(miss, collapse = ", ")))
  }

  offdiag <- cm[upper.tri(cm)]
  if (length(offdiag) && any(abs(offdiag) > 0.7, na.rm = TRUE)) {
    notes <- c(notes,
      "Some domain scores correlate above |0.7|. Under shared forcing this is expected, but it means influence cannot be attributed cleanly to a single domain.")
  }

  sds <- infl$sd[is.finite(infl$sd)]
  if (length(sds) > 1 && max(sds) / min(sds) > 1.5) {
    notes <- c(notes, sprintf(
      "Domain score dispersion is uneven (sd range %.3f-%.3f). Equal weights do not give equal influence when dispersion differs.",
      min(sds), max(sds)))
  }

  if (!length(notes)) {
    notes <- "Realised influence is close to the declared weights; no imbalance detected."
  }

  list(influence = infl, covariance = cm, notes = notes)
}
