#' Fit a conditional capacity-recovery curve (legacy function name)
#' @description Fits y(t) = B * (1 - A * exp(-r * t)) after a specified
#' disturbance peak. B is fixed from the observed baseline. A is a fractional
#' recovery deficit and r is a trajectory recovery rate. They are NOT the
#' accessibility alpha and reservoir exchange k in the accessible-capacity model.
#' @param Eh_stability Data frame containing EAC and optionally EDC.
#' @param id Aligned identifiers containing time and grouping columns if absent.
#' @param perturb_time Disturbance peak time; NULL detects the EAC minimum.
#' Detection is descriptive and can select noise or the last observation.
#' @param tau_unit Time unit: day, hour or week.
#' @param group_cols Columns defining one trajectory. Include plant_id when
#' plant-level trajectories are separate; duplicate times are rejected.
#' @param fit_edc Also fit increasing EDC recovery toward its own baseline.
#' This is inappropriate for a decreasing EDC trajectory; default FALSE.
#' @param min_points Minimum distinct post-peak times (at least three).
#' This is a computational threshold, not proof of identifiability.
#' @param verbose Print fit status.
#' @param control Named list passed to stats::nls.control.
#' @param baseline_end Explicit end of the baseline window, strictly before
#' the peak. NULL uses only the first observation time as a stated assumption.
#' @return An hrri_arch object containing estimates, conditional fit status,
#' fit objects and optional summaries of supplied latent columns. Legacy columns
#' Q_eac, alpha_eac, k_eac and M_eac remain for compatibility. Prefer the aliases
#' baseline_eac, deficit_fraction_eac, recovery_rate_eac and terminal_ratio_eac.
#' M_eac is a terminal-to-baseline ratio, not a measure of causal memory.
#' k_eac_h converts the trajectory rate to reciprocal hours, not exchange kinetics.
#' @details A nonlinear fit converging does not establish structural or practical
#' identifiability. Baseline uncertainty is not propagated into coefficient SEs.
#' With Cacc(t) = Q * alpha * (1 - exp(-k*t)), Q and alpha cannot be separated
#' using that curve alone. This function does not solve that inverse problem.
#' Fit only a single recovery window. Nonmonotonic and multi-event trajectories
#' require a different model and residual diagnostics.
#' @seealso rri_recovery_metrics, rri_accessible_capacity
#' @examples
#' tt <- 0:12
#' yy <- ifelse(tt <= 2, 10, 10 * (1 - 0.7 * exp(-0.3 * (tt - 3))))
#' dat <- data.frame(plot = "P1", depth = "D1", time = tt, EAC = yy)
#' a <- hrri_infer_architecture(dat, perturb_time = 3, baseline_end = 2,
#'                              verbose = FALSE)
#' a$estimates
#' @importFrom stats coef lm nls nls.control vcov median setNames sd
#' @importFrom utils modifyList head
#' @export
hrri_infer_architecture <- function(Eh_stability,
                                     id           = NULL,
                                     perturb_time = NULL,
                                     tau_unit     = c("day", "hour", "week"),
                                     group_cols   = c("plot", "depth"),
                                     fit_edc      = FALSE,
                                     min_points   = 3L,
                                     verbose      = TRUE,
                                     control      = list(),
                                     baseline_end = NULL) {

  # Assisted-by: OpenAI Codex. Conditional recovery fitting, not mechanistic inversion.
  tau_unit <- match.arg(tau_unit)
  if (!is.numeric(min_points) || length(min_points) != 1L ||
      !is.finite(min_points) || min_points < 3 || min_points != floor(min_points))
    stop("min_points must be an integer >= 3.")
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is.null(perturb_time) && !scalar(perturb_time)) stop("Invalid perturb_time.")
  if (!is.null(baseline_end) && !scalar(baseline_end)) stop("Invalid baseline_end.")
  if (is.null(baseline_end))
    warning("baseline_end absent: using the earliest time as the baseline.", call. = FALSE)
  if (!is.logical(fit_edc) || length(fit_edc) != 1L || is.na(fit_edc))
    stop("fit_edc must be TRUE or FALSE.");
  # ── HIDDEN columns: strip from inference, save for validation ────────────────
  HIDDEN_COLS <- c("alpha_accept", "k_accept_h", "alpha_donate", "k_donate_h",
                   "Cacc_total",   "Cacc_raw",   "Cacc_eac",   "Cacc_edc",
                   "latent_Physio","latent_Soil","latent_Micro")
  GT_SAVE     <- c("alpha_accept", "k_accept_h", "alpha_donate", "k_donate_h",
                   "Cacc_eac", "Cacc_edc", "Cacc_EAC", "Cacc_EDC")

  eh <- as.data.frame(Eh_stability)
  if (!nrow(eh) || anyDuplicated(names(eh))) stop("Require nonempty data with unique columns.")

  # Capture ground truth BEFORE stripping (inference uses none of it)
  present_gt <- intersect(GT_SAVE, names(eh))
  has_gt     <- length(present_gt) > 0
  gt_raw     <- if (has_gt) eh[, present_gt, drop = FALSE] else NULL

  # Strip hidden columns for all subsequent computation
  eh <- eh[, setdiff(names(eh), HIDDEN_COLS), drop = FALSE]

  # ── Validate required columns ────────────────────────────────────────────────
  if (!"EAC" %in% names(eh))
    stop("Eh_stability must contain column 'EAC'.")

  # Merge grouping + time columns from id if not already in Eh_stability
  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (nrow(id) != nrow(eh))
      stop("'id' must have the same number of rows as 'Eh_stability'.")
    for (col in intersect(c(group_cols, "time"), intersect(names(eh), names(id))))
      if (!identical(as.character(eh[[col]]), as.character(id[[col]])))
        stop("Conflicting identifiers in id and observations: ", col)
    for (col in setdiff(c(group_cols, "time"), names(eh))) {
      if (col %in% names(id)) eh[[col]] <- id[[col]]
    }
  }

  if (!"time" %in% names(eh))
    stop("Column 'time' not found in Eh_stability or id.")

  missing_gc <- setdiff(group_cols, names(eh))
  if (length(missing_gc))
    stop("group_cols not found in data: ", paste(missing_gc, collapse = ", "))

  if (!length(group_cols)) stop("Supply trajectory grouping columns.")
  if (!is.numeric(eh$time) || any(!is.finite(eh$time))) stop("time must be finite numeric.")
  if (anyDuplicated(.rri_key(eh, c(group_cols, "time"))))
    stop("Duplicate times within trajectory: include all experimental-unit keys.")
  eh[, intersect(c("EAC", "EDC"), names(eh))] <-
    .rri_numeric_df(eh[, intersect(c("EAC", "EDC"), names(eh)), drop = FALSE])
  # ── k unit conversion: tau_unit^-1  -->  h^-1 ────────────────────────────────
  k_to_h <- switch(tau_unit, day = 1 / 24, hour = 1, week = 1 / 168)

  # ────────────────────────────────────────────────────────────────────────────
  # Internal fitting helpers
  # ────────────────────────────────────────────────────────────────────────────

  # Log-linear OLS for k starting value
  .k_start <- function(x_rec, t_rel, Q_pre) {
    y  <- suppressWarnings(log(ifelse(1 - x_rec / Q_pre > 0 &
        1 - x_rec / Q_pre < 1, 1 - x_rec / Q_pre, NA_real_)))
    ok <- is.finite(y) & is.finite(t_rel) & t_rel > 0
    if (sum(ok) < 2L) return(0.05)
    b <- tryCatch(-stats::coef(stats::lm(y[ok] ~ t_rel[ok]))[["t_rel[ok]"]],
                  error = function(e) 0.05)
    max(b, 0.001)
  }

  # ── Core: fit recovery model for one group × direction (EAC or EDC) ──────────
  .fit_direction <- function(x_vec,       # observed capacity time series
                              time_vec,    # corresponding time values
                              t_nadir,     # disturbance peak (pre/post split)
                              min_pts,     # minimum post-disturbance points
                              ctrl_nls) {

    blank <- list(Q = NA_real_, alpha = NA_real_, k = NA_real_, k_h = NA_real_,
                  M = NA_real_, alpha_se = NA_real_, k_se = NA_real_,
                  n_pre = 0L, n_rec = 0L,
                  method = NA_character_, identifiability = "not_attempted",
                  fit_obj = NULL)

    ok <- !is.na(x_vec) & !is.na(time_vec) & is.finite(x_vec)
    if (sum(ok) < 2L) {
      blank$identifiability <- "no_valid_data"
      return(blank)
    }

    xv <- x_vec[ok]; tv <- time_vec[ok]

    b_end <- if (is.null(baseline_end)) min(tv) else baseline_end
    if (b_end >= t_nadir) {
      blank$identifiability <- "baseline_not_before_peak"
      return(blank)
    }
    pre_mask <- tv <= b_end
    rec_mask <- tv >  t_nadir
    blank$n_pre <- sum(pre_mask)
    blank$n_rec <- sum(rec_mask)

    # ── Q: pre-disturbance median capacity ───────────────────────────────────
    Q_pre <- if (blank$n_pre >= 1L) stats::median(xv[pre_mask]) else NA_real_
    if (!is.finite(Q_pre) || Q_pre <= 0) {
      blank$identifiability <- "insufficient_baseline"
      return(blank)
    }
    blank$Q <- Q_pre

    # ── Identifiability check ─────────────────────────────────────────────────
    if (blank$n_rec < min_pts) {
      blank$identifiability <- if (blank$n_rec == 0L) "no_recovery_window"
                               else "too_few_points"
      return(blank)   # Q only
    }

    x_rec <- xv[rec_mask]
    t_rel <- tv[rec_mask] - t_nadir   # time since nadir

    if (!any(x_rec < Q_pre) || .rri_cor(t_rel, x_rec) <= 0 ||
        !is.finite(.rri_cor(t_rel, x_rec))) {
      blank$identifiability <- "not_increasing_recovery"
      return(blank)
    }
    # ── Starting values ───────────────────────────────────────────────────────
    x_nadir  <- min(x_rec, na.rm = TRUE)
    alpha_s  <- min(max(1 - x_nadir / Q_pre, 0.01), 0.99)
    k_s      <- .k_start(x_rec, t_rel, Q_pre)

    # ── NLS (port algorithm with parameter bounds) ────────────────────────────
    # Model: EAC(t) = Q_pre * (1 - alpha * exp(-k * (t - t_nadir)))
    # Q_pre is fixed (pre-disturbance observable baseline, not a free parameter)
    nls_ok <- tryCatch({
      ctrl <- do.call(stats::nls.control,
                      utils::modifyList(list(maxiter = 200, tol = 1e-5, warnOnly = FALSE),
                                 ctrl_nls))
      fit <- stats::nls(
        x_rec ~ Q_pre * (1 - alpha * exp(-k * t_rel)),
        start     = list(alpha = alpha_s, k = k_s),
        lower     = c(alpha = 0.001, k = 1e-4),
        upper     = c(alpha = 1.000, k = 500),
        algorithm = "port",
        control   = ctrl
      )
      co <- stats::coef(fit)
      se <- tryCatch(sqrt(diag(stats::vcov(fit))),
                     error = function(e) c(alpha = NA_real_, k = NA_real_))
      blank$alpha    <- unname(co["alpha"])
      blank$k        <- unname(co["k"])
      blank$k_h      <- blank$k * k_to_h
      blank$alpha_se <- unname(se[["alpha"]])
      blank$k_se     <- unname(se[["k"]])
      blank$fit_obj  <- fit
      blank$method   <- "nls_port"
      blank$identifiability <- "fitted_conditional"
      TRUE
    }, error = function(e) FALSE)

    # ── Fallback: log-linear OLS ───────────────────────────────────────────────
    if (!nls_ok) {
      y_ll <- suppressWarnings(log(ifelse(1 - x_rec / Q_pre > 0 &
        1 - x_rec / Q_pre < 1, 1 - x_rec / Q_pre, NA_real_)))
      fin  <- is.finite(y_ll) & is.finite(t_rel) & t_rel > 0
      if (sum(fin) >= min_pts) {
        ols <- tryCatch(stats::coef(stats::lm(y_ll[fin] ~ t_rel[fin])),
                        error = function(e) c(NA_real_, NA_real_))
        if (length(ols) >= 2L && !is.na(ols[2L]) && ols[2L] < 0 &&
            is.finite(exp(ols[1L])) && exp(ols[1L]) > 0 && exp(ols[1L]) <= 1) {
          blank$k     <- -ols[2L]
          blank$k_h   <- blank$k * k_to_h
          blank$alpha <- unname(exp(ols[1L]))
          blank$method         <- "log_linear_fallback"
          blank$identifiability <- "fitted_log_linear"
        } else {
          blank$identifiability <- "nls_failed"
        }
      } else {
        blank$identifiability <- "nls_failed"
      }
    }

    # ── M: memory (post-recovery plateau / pre-disturbance baseline) ──────────
    # Use last 20% of recovery observations as a plateau estimate.
    if (blank$n_rec >= 3L) {
      n_tail   <- max(1L, round(blank$n_rec * 0.2))
      tail_idx <- order(t_rel, decreasing = TRUE)[seq_len(n_tail)]
      plateau  <- stats::median(x_rec[tail_idx], na.rm = TRUE)
      if (is.finite(plateau) && Q_pre > 0)
        blank$M <- plateau / Q_pre
    }

    blank
  }

  # ── Group nadir detector ─────────────────────────────────────────────────────
  .detect_nadir <- function(x_vec, t_vec) {
    ok <- !is.na(x_vec) & is.finite(x_vec)
    if (!any(ok)) return(NA_real_)
    t_vec[ok][which.min(x_vec[ok])]
  }

  # ────────────────────────────────────────────────────────────────────────────
  # Split data into groups and fit each
  # ────────────────────────────────────────────────────────────────────────────
  group_key <- .rri_key(eh, group_cols)
  unique_groups <- unique(group_key)

  if (verbose)
    message("hrri_infer_architecture: ", length(unique_groups),
            " groups | tau_unit='", tau_unit,
            "' | min_points=", min_points,
            " | fit_edc=", fit_edc)

  ctrl_nls <- utils::modifyList(list(maxiter = 200, tol = 1e-5), control)

  raw_results <- lapply(unique_groups, function(gk) {
    idx <- which(group_key == gk)
    sub <- eh[idx, , drop = FALSE]

    # Group identifier values (first-row representatives)
    gv <- stats::setNames(lapply(group_cols, function(g) sub[[g]][1L]), group_cols)

    time_v <- sub$time
    eac_v  <- sub$EAC
    edc_v  <- if (fit_edc && "EDC" %in% names(sub)) sub$EDC else NULL

    t_nadir_g <- if (!is.null(perturb_time)) {
      perturb_time
    } else {
      .detect_nadir(eac_v, time_v)
    }

    if (is.na(t_nadir_g)) {
      na_fit <- list(Q = NA_real_, alpha = NA_real_, k = NA_real_, k_h = NA_real_,
                     M = NA_real_, alpha_se = NA_real_, k_se = NA_real_,
                     n_pre = 0L, n_rec = 0L, method = NA_character_,
                     identifiability = "no_valid_data", fit_obj = NULL)
      return(list(gk = gk, gv = gv, t_perturb = NA_real_,
                  eac = na_fit, edc = NULL, idx = idx))
    }

    list(gk        = gk,
         gv        = gv,
         t_perturb = t_nadir_g,
         eac       = .fit_direction(eac_v, time_v, t_nadir_g, min_points, ctrl_nls),
         edc       = if (!is.null(edc_v))
                       .fit_direction(edc_v, time_v, t_nadir_g, min_points, ctrl_nls)
                     else NULL,
         idx       = idx)
  })

  # ────────────────────────────────────────────────────────────────────────────
  # Assemble estimates data frame (one row per group)
  # ────────────────────────────────────────────────────────────────────────────
  est_list <- lapply(raw_results, function(r) {
    gv <- as.data.frame(r$gv, stringsAsFactors = FALSE)
    e  <- r$eac
    d  <- r$edc   # NULL when fit_edc=FALSE or EDC absent

    data.frame(
      gv,
      t_perturb    = r$t_perturb,
      # EAC estimates
      Q_eac        = e$Q,
      alpha_eac    = e$alpha,
      k_eac        = e$k,
      k_eac_h      = e$k_h,
      alpha_eac_se = e$alpha_se,
      k_eac_se     = e$k_se,
      M_eac        = e$M,
      n_pre        = e$n_pre,
      n_rec        = e$n_rec,
      method_eac   = if (!is.null(e$method)) e$method else NA_character_,
      ident_eac    = if (!is.null(e$identifiability)) e$identifiability
                     else NA_character_,
      # EDC estimates (NA when not fitted)
      Q_edc        = if (!is.null(d)) d$Q     else NA_real_,
      alpha_edc    = if (!is.null(d)) d$alpha else NA_real_,
      k_edc        = if (!is.null(d)) d$k     else NA_real_,
      k_edc_h      = if (!is.null(d)) d$k_h   else NA_real_,
      M_edc        = if (!is.null(d)) d$M     else NA_real_,
      ident_edc    = if (!is.null(d)) d$identifiability else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  estimates <- do.call(rbind, est_list)
  rownames(estimates) <- NULL
  estimates$baseline_eac <- estimates$Q_eac
  estimates$deficit_fraction_eac <- estimates$alpha_eac
  estimates$recovery_rate_eac <- estimates$k_eac
  estimates$terminal_ratio_eac <- estimates$M_eac

  # ────────────────────────────────────────────────────────────────────────────
  # Ground truth: one row per group, keyed on group_cols (for Fig 3 merge)
  # ────────────────────────────────────────────────────────────────────────────
  ground_truth <- if (has_gt) {
    gt_list <- lapply(raw_results, function(r) {
      gv <- as.data.frame(r$gv, stringsAsFactors = FALSE)
      # Hidden states vary through time; report group means, never a first-row oracle.
      cbind(gv, as.data.frame(lapply(gt_raw[r$idx, , drop = FALSE], .rri_mean)))
    })
    gt_out <- do.call(rbind, gt_list)
    rownames(gt_out) <- NULL
    gt_out
  } else NULL

  # ────────────────────────────────────────────────────────────────────────────
  # Identifiability summary
  # ────────────────────────────────────────────────────────────────────────────
  ident_tbl <- table(estimates$ident_eac)
  n_full    <- sum(estimates$ident_eac %in% c("fitted_conditional", "fitted_log_linear"),
                   na.rm = TRUE)
  n_q_only  <- sum(estimates$ident_eac %in%
                     c("too_few_points", "no_recovery_window"), na.rm = TRUE)
  n_total   <- nrow(estimates)
  n_failed  <- n_total - n_full - n_q_only
  ident_msg <- sprintf(
    "%d/%d groups: %d conditional curve fits | %d baseline-only | %d failed",
    n_full, n_total, n_full, n_q_only, n_failed
  )
  if (verbose) message("  ", ident_msg)
  if (has_gt && verbose)
    message("  Ground truth saved ($ground_truth): ",
            paste(present_gt, collapse = ", "))

  ident_summary <- list(
    n_groups = n_total,
    n_full   = n_full,
    n_q_only = n_q_only,
    n_failed = n_failed,
    table    = ident_tbl,
    message  = ident_msg
  )

  # ────────────────────────────────────────────────────────────────────────────
  # Return object
  # ────────────────────────────────────────────────────────────────────────────
  structure(
    list(
      estimates        = estimates,
      ground_truth     = ground_truth,
      has_ground_truth = has_gt,
      identifiability  = ident_summary,
      tau_unit         = tau_unit,
      group_cols       = group_cols,
      perturb_time     = perturb_time,   # NULL = auto-detected per group
      min_points       = min_points,
      fit_edc          = fit_edc,
      raw_results      = raw_results,
      baseline_end     = baseline_end,
      interpretation   = "Conditional recovery curve; not Q-alpha-k-memory identification"
    ),
    class = "hrri_arch"
  )
}


# ── print.hrri_arch ─────────────────────────────────────────────────────────────
#' @rdname hrri_infer_architecture
#' @param x An hrri_arch object.
#' @param ... Additional method arguments.
#' @export
print.hrri_arch <- function(x, ...) {
  cat("<hrri_arch>  Conditional capacity-recovery fit\n")
  cat(sprintf("  Groups  : %d  |  tau_unit: %s  |  fit_edc: %s  |  min_points: %d\n",
              x$identifiability$n_groups, x$tau_unit,
              x$fit_edc, x$min_points))
  cat(sprintf("  %s\n", x$identifiability$message))
  if (x$has_ground_truth)
    cat("  Ground truth available  -->  merge(arch$estimates, arch$ground_truth)\n")

  cols_show <- intersect(
    c(x$group_cols, "Q_eac", "alpha_eac", "k_eac_h", "M_eac", "ident_eac"),
    names(x$estimates)
  )
  cat("\nEstimates (first 6 groups):\n")
  print(utils::head(x$estimates[, cols_show, drop = FALSE], 6L),
        row.names = FALSE, digits = 4)
  invisible(x)
}


# ── summary.hrri_arch ───────────────────────────────────────────────────────────
#' @rdname hrri_infer_architecture
#' @param object An hrri_arch object.
#' @export
summary.hrri_arch <- function(object, ...) {
  cat("<hrri_arch>  Summary of conditional recovery fits\n")
  cat(sprintf("  tau_unit: '%s'  |  groups: %d  |  fit_edc: %s\n",
              object$tau_unit, object$identifiability$n_groups, object$fit_edc))
  cat("\nIdentifiability (EAC):\n")
  print(object$identifiability$table)

  e <- object$estimates
  cat("\nParameter distributions (EAC, conditionally fitted groups):\n")
  for (param in c("Q_eac", "alpha_eac", "k_eac_h", "M_eac")) {
    v <- e[[param]][is.finite(e[[param]]) &
                      e$ident_eac %in% c("fitted_conditional", "fitted_log_linear")]
    if (length(v) >= 2L)
      cat(sprintf("  %-12s  n=%-4d  mean=%8.4f  sd=%7.4f  [%7.4f, %7.4f]\n",
                  param, length(v), mean(v), stats::sd(v), min(v), max(v)))
  }

  if (object$has_ground_truth)
    cat("\nHidden-state summaries are not ground truth for these recovery-curve parameters.\n")

  invisible(object)
}


# ── as.data.frame.hrri_arch ─────────────────────────────────────────────────────
#' @rdname hrri_infer_architecture
#' @param row.names,optional Standard data-frame method arguments.
#' @export
as.data.frame.hrri_arch <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- x$estimates
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}


# ── validate_architecture ───────────────────────────────────────────────────────
#' Compare explicitly matched estimands
#' @rdname hrri_infer_architecture
#' @description Tabulates estimates and independently supplied targets.
#' There is no default mapping from recovery-curve amplitude to accessibility
#' or from recovery rate to reservoir exchange kinetics.
#' @param arch An hrri_arch object with a keyed ground_truth data frame.
#' @param params Explicit strings of the form estimated_column:true_column.
#' @return Data frame with finite_pair flags and descriptive aggregate statistics.
#' The caller must establish that paired columns have the same definition and units.
#' @export
validate_architecture <- function(arch, params = NULL) {
  if (is.null(params) || !is.character(params) || !length(params) || anyNA(params))
    stop("Supply explicit matching estimands in params; recovery parameters are not mechanistic alpha/k.")
  if (any(grepl(":(alpha_accept|alpha_donate|k_accept_h|k_donate_h)$", params)))
    stop("Recovery amplitude/rate cannot validate mechanistic accessibility/exchange rates.")
  if (!inherits(arch, "hrri_arch"))
    stop("`arch` must be an hrri_arch object.")
  if (!arch$has_ground_truth)
    stop("No ground truth available. Run hrri_infer_architecture() on ",
         "simulator output that contains hidden columns (alpha_accept, etc.).")

  val <- merge(arch$estimates, arch$ground_truth,
               by = arch$group_cols, all.x = TRUE)

  rows <- lapply(params, function(p) {
    parts <- strsplit(p, ":")[[1L]]
    if (length(parts) != 2L)
      stop("Each `params` entry must be 'estimated_col:true_col', got: ", p)
    est_col <- parts[1L]; true_col <- parts[2L]
    if (!est_col  %in% names(val))
      stop("Column '", est_col,  "' not in estimates.")
    if (!true_col %in% names(val))
      stop("Column '", true_col, "' not in ground_truth.")

    gv   <- val[, arch$group_cols, drop = FALSE]
    hat  <- val[[est_col]]; tru <- val[[true_col]]
    ok   <- is.finite(hat) & is.finite(tru)
    res  <- hat - tru
    cbind(
      gv,
      data.frame(
        parameter  = p,
        estimated  = hat,
        true       = tru,
        residual   = res,
        finite_pair = ok,
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Aggregate stats per parameter as attribute
  agg <- do.call(rbind, lapply(params, function(p) {
    sub <- out[out$parameter == p & out$finite_pair, ]
    if (nrow(sub) < 3L) return(NULL)
    data.frame(
      parameter = p,
      r         = .rri_cor(sub$estimated, sub$true),
      RMSE      = sqrt(mean(sub$residual^2)),
      bias      = mean(sub$residual),
      n         = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
  attr(out, "aggregate") <- agg
  class(out) <- c("hrri_arch_validation", "data.frame")
  out
}


# ── print.hrri_arch_validation ──────────────────────────────────────────────────
#' @rdname hrri_infer_architecture
#' @export
print.hrri_arch_validation <- function(x, ...) {
  cat("<hrri_arch_validation>  Estimated vs ground truth\n\n")
  agg <- attr(x, "aggregate")
  if (!is.null(agg)) {
    cat("Aggregate statistics:\n")
    print(agg, row.names = FALSE, digits = 4)
    cat("\n")
  }
  NextMethod()
}
