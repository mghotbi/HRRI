#' Descriptive recovery metrics for a single disturbance
#'
#' @description Summarises decline and return of a higher-is-better score.
#' A score decline does not identify pathway truncation; a displaced plateau
#' does not establish alternative electron routing. Hysteresis is only reported
#' for a sufficiently closed, reversing forcing-response path. Temporal deficit
#' asymmetry is a separate diagnostic. Analyse repeated events separately.
#' @param res RRI object or data frame.
#' @param id Optional aligned identifiers; common columns must agree.
#' @param time_col Numeric time column; time must be unique within each group.
#' @param group_cols Columns identifying one longitudinal experimental unit.
#' @param perturb_start,perturb_end Finite start and end of one disturbance.
#' @param rri_col Numeric score column.
#' @param forcing_col Optional measured external forcing column, not a response proxy.
#' @param min_pts Minimum finite baseline and recovery observations.
#' @param lag_threshold Fraction of observed decline defining recovery onset.
#' @param plateau_window Number of final observations for plateau assessment.
#' @param plateau_tol Fractional terminal displacement defining a plateau flag.
#' @return One row per group, including diagnostic fit status and observation counts.
#' k is a log-linear fit of positive baseline deficits after the observed minimum.
#' It is a conditional trajectory descriptor, not a mechanistic exchange rate.
#' Legacy alt_routing fields are retained as NA; use displaced_plateau_flag.
#' @importFrom stats lm coef sd
#' @importFrom utils head tail
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
#'                       id = sim$id)
#'   rri_recovery_metrics(res, time_col = "time",
#'     group_cols = c("plot", "depth", "plant_id"),
#'     perturb_start = 20, perturb_end = 30)
#' }
#' @export
rri_recovery_metrics <- function(res, id = NULL, time_col = "time", group_cols = NULL,
  perturb_start, perturb_end, rri_col = "RRI", forcing_col = NULL,
  min_pts = 3L, lag_threshold = 0.05, plateau_window = 3L, plateau_tol = 0.10) {
  df <- if (is.list(res) && !is.data.frame(res) && !is.null(res$row_scores))
    as.data.frame(res$row_scores) else as.data.frame(res)
  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (nrow(id) != nrow(df)) stop("id and scores must have equal row counts.")
    for (nm in intersect(names(id), names(df))) {
      if (!isTRUE(all.equal(id[[nm]], df[[nm]], check.attributes = FALSE)))
        stop("id and scores disagree in column: ", nm)
    }
    df <- cbind(df, id[, setdiff(names(id), names(df)), drop = FALSE])
  }
  if (!all(c(time_col, rri_col, group_cols, forcing_col) %in% names(df)))
    stop("Requested score, time, group or forcing columns are absent.")
  if (!is.numeric(df[[time_col]]) || !is.numeric(df[[rri_col]]) ||
      any(!is.finite(df[[time_col]]))) stop("Time must be finite numeric; score must be numeric.")
  if (length(perturb_start) != 1L || length(perturb_end) != 1L ||
      any(!is.finite(c(perturb_start, perturb_end))) || perturb_end <= perturb_start)
    stop("Require finite perturb_start < perturb_end.")
  for (v in list(min_pts, plateau_window)) if (length(v) != 1L || !is.finite(v) ||
      v < 2 || v != as.integer(v)) stop("min_pts and plateau_window must be integers >= 2.")
  if (length(lag_threshold) != 1L || !is.finite(lag_threshold) ||
      lag_threshold <= 0 || lag_threshold >= 1 || length(plateau_tol) != 1L ||
      !is.finite(plateau_tol) || plateau_tol <= 0) stop("Invalid detection thresholds.")
  key <- .rri_key(df, group_cols)
  auc <- function(tt, yy) if (length(tt) < 2L) NA_real_ else
    sum(diff(tt) * (utils::head(yy, -1) + utils::tail(yy, -1)) / 2)
  out <- lapply(unique(key), function(g) {
    sub <- df[key == g, , drop = FALSE]
    sub <- sub[order(sub[[time_col]]), , drop = FALSE]
    t <- sub[[time_col]]; y <- sub[[rri_col]]
    if (anyDuplicated(t)) stop("Duplicate times within a trajectory. Add unit identifiers or aggregate explicitly.")
    finite <- is.finite(y)
    pre <- which(t < perturb_start & finite)
    ptb <- which(t >= perturb_start & t <= perturb_end & finite)
    rec <- which(t > perturb_end & finite)
    post <- which(t >= perturb_start & finite)
    z <- data.frame(baseline_rri = NA_real_, min_rri = NA_real_, depth_min = NA_real_,
      depth_min_frac = NA_real_, tau_lag = NA_real_, k_recovery = NA_real_,
      t_half = NA_real_, overshoot = NA_real_, overshoot_frac = NA_real_,
      H_hysteresis = NA_real_, H_axis = if (is.null(forcing_col)) "unavailable" else "forcing",
      temporal_asymmetry = NA_real_, incomplete_return = NA_real_,
      incomplete_return_frac = NA_real_, displaced_plateau_flag = NA,
      displaced_plateau_level = NA_real_, alt_routing_flag = NA,
      alt_routing_level = NA_real_, n_pre = length(pre), n_perturb = length(ptb),
      n_recovery = length(rec), n_missing = sum(!finite), n_fit = 0L,
      fit_status = "insufficient_data", fit_r_squared = NA_real_,
      fit_start_time = NA_real_, final_observation_time = if (length(rec)) max(t[rec]) else NA_real_,
      hysteresis_status = "not_evaluated", stringsAsFactors = FALSE)
    if (length(pre) >= min_pts && length(rec) >= min_pts && length(post)) {
      b <- mean(y[pre]); mn <- min(y[post]); dep <- max(0, b - mn)
      z$baseline_rri <- b; z$min_rri <- mn; z$depth_min <- dep
      z$depth_min_frac <- if (b != 0) dep / abs(b) else NA_real_
      z$overshoot <- max(0, max(y[rec]) - b)
      z$overshoot_frac <- if (b != 0) z$overshoot / abs(b) else NA_real_
      final <- mean(y[utils::tail(rec, min(plateau_window, length(rec)))])
      z$incomplete_return <- final - b
      z$incomplete_return_frac <- if (b != 0) (final - b) / abs(b) else NA_real_
      z$fit_status <- "no_resolvable_decline"
      if (dep > sqrt(.Machine$double.eps) * max(1, abs(b))) {
        if (length(ptb)) {
          y_end <- y[utils::tail(ptb, 1)]
          onset <- rec[y[rec] >= y_end + lag_threshold * dep]
          if (length(onset)) z$tau_lag <- t[onset[1]] - perturb_end
        }
        trough_t <- t[post[which.min(y[post])]]
        z$fit_start_time <- max(perturb_end, trough_t)
        candidate <- rec[t[rec] >= z$fit_start_time]
        # Do not join separated declining segments across a baseline crossing.
        crossing <- which(y[candidate] >= b)
        if (length(crossing)) candidate <- utils::head(candidate, crossing[1] - 1L)
        fit_idx <- candidate[b - y[candidate] > .Machine$double.eps * max(1, abs(b))]
        z$n_fit <- length(fit_idx); z$fit_status <- "insufficient_positive_deficits"
        if (length(fit_idx) >= min_pts) {
          fit <- stats::lm(log(b - y[fit_idx]) ~ t[fit_idx])
          rate <- -unname(stats::coef(fit)[2])
          z$fit_r_squared <- summary(fit)$r.squared
          if (is.finite(rate) && rate > 0) {
            z$k_recovery <- rate; z$t_half <- log(2) / rate
            z$fit_status <- "fitted_conditional_exponential"
          } else z$fit_status <- "nonrecovering_or_flat"
        }
      }
      didx <- c(utils::tail(pre, 1), ptb)
      ridx <- c(utils::tail(ptb, 1), rec)
      ad <- auc(t[didx], pmax(0, b - y[didx]))
      ar <- auc(t[ridx], pmax(0, b - y[ridx]))
      if (is.finite(ad) && is.finite(ar) && ad + ar > 0)
        z$temporal_asymmetry <- (ad - ar) / (ad + ar)
      if (!is.null(forcing_col)) {
        if (!is.numeric(sub[[forcing_col]])) stop("Forcing must be numeric.")
        ix <- c(utils::tail(pre, 1), post)
        fx <- sub[[forcing_col]][ix]; ry <- y[ix]
        z$hysteresis_status <- "missing_or_degenerate_path"
        if (length(ix) >= 4L && all(is.finite(fx)) && all(is.finite(ry))) {
          xr <- diff(range(fx)); yr <- diff(range(ry))
          reversing <- any(diff(fx) > 0) && any(diff(fx) < 0)
          closed <- xr > 0 && yr > 0 && abs(utils::tail(fx, 1) - fx[1]) <= 0.05 * xr &&
            abs(utils::tail(ry, 1) - ry[1]) <= 0.05 * yr
          if (reversing && closed) {
            z$H_hysteresis <- 0.5 * sum(fx * c(ry[-1], ry[1]) -
              ry * c(fx[-1], fx[1])) / (xr * yr)
            z$hysteresis_status <- "closed_path_area_descriptive"
          } else z$hysteresis_status <- "open_or_nonreversing_path"
        }
      }
      if (length(rec) >= plateau_window && b != 0) {
        idx <- utils::tail(rec, plateau_window)
        slope <- unname(stats::coef(stats::lm(y[idx] ~ t[idx]))[2])
        stable <- stats::sd(y[idx]) < 0.05 * abs(b) &&
          abs(slope) * diff(range(t[idx])) < 0.05 * abs(b)
        z$displaced_plateau_flag <- isTRUE(stable && abs(final - b) / abs(b) > plateau_tol)
        if (z$displaced_plateau_flag) z$displaced_plateau_level <- final
      }
    }
    if (length(group_cols)) cbind(sub[1, group_cols, drop = FALSE], z) else z
  })
  if (!length(out)) stop("No trajectories supplied.")
  out <- do.call(rbind, out); rownames(out) <- NULL
  out$k <- out$k_recovery; out$H <- out$H_hysteresis
  out$I <- abs(out$incomplete_return); out$H_abs <- abs(out$H_hysteresis)
  out$H_norm <- pmin(1, out$H_abs)
  out$I_norm <- pmin(1, abs(out$incomplete_return_frac))
  out
}
