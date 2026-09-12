#' @title Agreement between a score and a reference target, respecting clustering
#'
#' @description
#' Quantifies how closely a score tracks a reference target, using statistics
#' appropriate to repeated observations of the same experimental units. A
#' pooled row-wise correlation over a longitudinal panel overstates precision,
#' because rows within a trajectory are not independent. This function reports
#' the naive row-level result alongside the cluster-aware one, so the
#' difference is visible rather than hidden.
#'
#' @param score Numeric vector of scores, or an `RRI` object. If an `RRI`
#'   object, `score_col` is taken from its `row_scores`.
#' @param target Numeric vector of reference values, the same length as
#'   `score`.
#' @param cluster Vector identifying the independent experimental unit for each
#'   observation, typically one trajectory. Rows sharing a value are treated as
#'   dependent. If `NULL`, every row is treated as independent and the function
#'   warns, because that assumption is rarely correct for time series.
#' @param score_col Column name used when `score` is an `RRI` object.
#' @param n_boot Number of cluster bootstrap resamples for confidence
#'   intervals. Set to 0 to skip.
#' @param n_perm Number of cluster permutations for the null test. Set to 0 to
#'   skip.
#' @param conf Confidence level for intervals.
#' @param seed Optional integer seed for reproducible resampling.
#'
#' @return An object of class `rri_accuracy`: a list with elements `agreement`,
#'   `calibration`, `decomposition`, `dependence`, `null_test` and `notes`,
#'   plus `data` (the complete cases actually used), `draws` (the resampled
#'   and permuted statistics) and `conf`. The last three exist so that
#'   [plot_rri_accuracy()] can draw the sampling distributions without
#'   repeating the resampling. See Details.
#'
#' @details
#' **Why correlation is not agreement.** Pearson's \eqn{r} is invariant to
#' location and scale: a score equal to \eqn{2 \times} the target plus a
#' constant correlates perfectly with it while agreeing with it nowhere. Lin's
#' concordance correlation coefficient
#'
#' \deqn{\rho_c = \frac{2 s_{xy}}{s_x^2 + s_y^2 + (\bar{x} - \bar{y})^2}}
#'
#' penalises departure from the 1:1 line and is reported alongside \eqn{r}. A
#' large gap between the two means the score is well correlated but
#' miscalibrated.
#'
#' **Why clustering matters.** With \eqn{m} observations per unit and
#' intra-cluster correlation \eqn{\rho}, the design effect is
#' \eqn{1 + (m - 1)\rho} and the effective sample size is \eqn{n / \mathrm{deff}}.
#' For a 40-point trajectory with \eqn{\rho = 0.3} this is roughly a twelvefold
#' reduction, so a confidence interval computed from the row count is far too
#' narrow. Intervals here come from a cluster bootstrap that resamples whole
#' trajectories with replacement, which preserves the within-unit dependence.
#'
#' **Error decomposition.** Mean squared error is split following Kobayashi and
#' Salam (2000) into squared bias, a difference in variability, and a lack of
#' correlation:
#'
#' \deqn{\mathrm{MSE} = (\bar{x} - \bar{y})^2 + (s_x - s_y)^2 + 2 s_x s_y (1 - r)}
#'
#' These say different things. Large squared bias means a systematic offset,
#' correctable by recentring. Large variance mismatch means the score is too
#' flat or too volatile. Large lack of correlation means the score does not
#' track the target's pattern, and no rescaling will fix it.
#'
#' **What the permutation null asks.** Whole trajectories are exchanged, so
#' each unit keeps its own temporal shape and only the pairing between score
#' and target is broken. This is deliberately the harder null. Permuting
#' individual rows would destroy the shared event-driven shape that every
#' trajectory has, making almost any score look significant; exchanging blocks
#' retains that shape and asks whether the score tracks *this* unit's target
#' beyond what the common disturbance imposes on all of them. A small p-value
#' under this null is therefore informative, and a large one is not evidence
#' that the score is uninformative about the disturbance itself.
#'
#' **What this does not establish.** If the score and the target come from the
#' same generative model, as they do for `latent_truth` from
#' [simulate_redox_holobiont()], high agreement is internal consistency and
#' nothing more. It is not predictive accuracy, not out-of-sample error, and
#' not evidence of ecological validity. Those require a target constructed
#' independently of the score, and replication at the level of independent
#' experimental units.
#'
#' @references
#' Lin, L.I. (1989) A concordance correlation coefficient to evaluate
#' reproducibility. *Biometrics*, **45**, 255--268.
#'
#' Kobayashi, K. & Salam, M.U. (2000) Comparing simulated and measured values
#' using mean squared deviation and its components. *Agronomy Journal*,
#' **92**, 345--352.
#'
#' @examples
#' ## Synthetic panel: 8 trajectories, 20 time points each.
#' ## The score is deliberately miscalibrated so the r-vs-CCC gap is visible.
#' set.seed(1)
#' k <- 8; m <- 20
#' unit   <- rnorm(k, 0, 0.30)
#' target <- unlist(lapply(unit, function(u) u + 0.5 + rnorm(m, 0, 0.05)))
#' score  <- 0.75 * target + 0.10 + rnorm(k * m, 0, 0.06)
#' traj   <- rep(seq_len(k), each = m)
#'
#' acc <- rri_accuracy(score, target, cluster = traj,
#'                     n_boot = 200, n_perm = 200, seed = 1)
#' acc
#'
#' ## Correlation is high, concordance is not: the score compresses the target.
#' acc$agreement
#' acc$calibration
#'
#' ## The row count is not the sample size.
#' acc$dependence
#'
#' \donttest{
#' ## Against the simulator's own prescribed target.
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2, n_depth = 2, n_plant = 3, n_time = 30,
#'   p_micro = 20, seed = 2026
#' )
#' res <- suppressWarnings(rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data, id = sim$id
#' ))
#' scored <- attach_hrri_ids(res$row_scores, sim$id)
#' tj <- interaction(scored$plot, scored$depth, scored$plant_id, drop = TRUE)
#' rri_accuracy(scored$RRI, sim$latent_truth, cluster = tj,
#'              n_boot = 500, n_perm = 500, seed = 1)
#' }
#'
#' @seealso [plot_rri_accuracy()] for the four-panel diagnostic figure;
#'   [benchmark_hrri()] for repeated-seed benchmarking;
#'   [rri_domain_influence()] for which domain drives the score.
#' @importFrom stats var sd cor cov lm coef quantile setNames complete.cases
#' @export
rri_accuracy <- function(score, target, cluster = NULL,
                         score_col = "RRI",
                         n_boot = 1000, n_perm = 1000,
                         conf = 0.95, seed = NULL) {

  if (is.list(score) && !is.data.frame(score) && !is.null(score$row_scores)) {
    rs <- as.data.frame(score$row_scores)
    if (!score_col %in% names(rs))
      stop("`score_col` not found in row_scores: ", score_col, call. = FALSE)
    score <- rs[[score_col]]
  }
  score  <- suppressWarnings(as.numeric(score))
  target <- suppressWarnings(as.numeric(target))
  if (length(score) != length(target))
    stop("`score` and `target` must have the same length.", call. = FALSE)

  notes <- character(0)
  if (is.null(cluster)) {
    cluster <- seq_along(score)
    warning("No `cluster` supplied; every row treated as independent. ",
            "For longitudinal data this overstates precision.", call. = FALSE)
    notes <- c(notes,
      "No clustering supplied. Intervals assume independent rows and are almost certainly too narrow for repeated observations.")
  }
  cluster <- as.character(cluster)

  ok <- is.finite(score) & is.finite(target)
  score <- score[ok]; target <- target[ok]; cluster <- cluster[ok]
  n <- length(score)
  if (n < 4L) stop("Fewer than four complete observations.", call. = FALSE)

  ## Seeding without restoring is a modification of the user's global
  ## environment: set.seed() writes .Random.seed into .GlobalEnv, so a call
  ## here would silently reset the stream a caller was relying on. Save the
  ## RNG state and put it back on exit, matching the pattern already used in
  ## simulate_redox_holobiont() and benchmark_hrri().
  if (!is.null(seed)) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  ## ---- core statistics ----------------------------------------------------
  .stats <- function(x, y) {
    ## x = score, y = target.
    ## Lin's CCC and the MSE decomposition are both defined on POPULATION
    ## moments (divisor n). Using sd() (divisor n - 1) leaves the three
    ## decomposition components failing to sum to MSE by O(1/n), which looks
    ## like an arithmetic error in the output table. Population moments are
    ## used throughout for these, and sd() only for the reported dispersion.
    if (length(x) < 3L) return(NULL)
    nn <- length(x)
    mx <- mean(x); my <- mean(y)
    px <- sqrt(mean((x - mx)^2))            # population sd of score
    py <- sqrt(mean((y - my)^2))            # population sd of target
    pxy <- mean((x - mx) * (y - my))        # population covariance
    r  <- if (px > 0 && py > 0) pxy / (px * py) else NA_real_
    ccc <- {
      den <- px^2 + py^2 + (mx - my)^2
      if (den > 0) 2 * pxy / den else NA_real_
    }
    e   <- x - y
    mse <- mean(e^2)
    c(r = r, ccc = ccc,
      rmse = sqrt(mse), mae = mean(abs(e)), bias = mean(e),
      ## R2 against the target's own mean as the reference predictor.
      ## Negative values are meaningful: worse than predicting that mean.
      r2 = if (py > 0) 1 - mse / py^2 else NA_real_,
      sb   = (mx - my)^2,               # squared bias
      sdsd = (px - py)^2,               # variance mismatch
      lcs  = 2 * px * py * (1 - r),     # lack of correlation
      sd_score = stats::sd(x), sd_target = stats::sd(y), n = nn
    )
  }

  point <- .stats(score, target)

  ## ---- dependence structure ------------------------------------------------
  ug <- unique(cluster); k <- length(ug)
  sizes <- as.numeric(table(cluster))
  icc <- NA_real_; deff <- NA_real_; n_eff <- n
  if (k > 1L && all(sizes > 1L)) {
    m   <- mean(sizes)
    gm  <- mean(score)
    gmu <- vapply(ug, function(g) mean(score[cluster == g]), numeric(1))
    msb <- sum(sizes * (gmu - gm)^2) / (k - 1)
    msw <- sum(vapply(seq_along(ug), function(i)
      sum((score[cluster == ug[i]] - gmu[i])^2), numeric(1))) / (n - k)
    if (is.finite(msb) && is.finite(msw) && (msb + (m - 1) * msw) > 0) {
      icc  <- (msb - msw) / (msb + (m - 1) * msw)
      icc  <- max(0, min(1, icc))
      deff <- 1 + (m - 1) * icc
      n_eff <- n / deff
    }
  }

  ## ---- cluster-mean (between-unit) analysis -------------------------------
  gs <- vapply(ug, function(g) mean(score[cluster == g]), numeric(1))
  gt <- vapply(ug, function(g) mean(target[cluster == g]), numeric(1))
  cluster_point <- if (k >= 3L) .stats(gs, gt) else NULL

  ## ---- cluster bootstrap ---------------------------------------------------
  ## Resampled values are retained so plot_rri_accuracy() can draw the
  ## sampling distributions rather than re-running the resampling itself.
  draws   <- list()
  boot_ci <- NULL
  if (n_boot > 0 && k > 2L) {
    B <- matrix(NA_real_, nrow = n_boot, ncol = length(point),
                dimnames = list(NULL, names(point)))
    for (b in seq_len(n_boot)) {
      pick <- sample(ug, k, replace = TRUE)
      idx  <- unlist(lapply(pick, function(g) which(cluster == g)),
                     use.names = FALSE)
      s <- .stats(score[idx], target[idx])
      if (!is.null(s)) B[b, ] <- s
    }
    a <- (1 - conf) / 2
    boot_ci <- t(apply(B, 2, stats::quantile, probs = c(a, 1 - a),
                       na.rm = TRUE))
    colnames(boot_ci) <- c("lower", "upper")
    draws$cluster <- B
  }

  ## ---- naive row bootstrap, for contrast only ------------------------------
  ## Reported so the cost of ignoring clustering is visible rather than
  ## asserted. This interval is the WRONG one for longitudinal data; it is
  ## shown only to quantify how much too narrow it is.
  naive_ci <- NULL
  if (n_boot > 0) {
    v <- rep(NA_real_, n_boot)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      s <- .stats(score[idx], target[idx])
      if (!is.null(s)) v[b] <- s[["r"]]
    }
    a <- (1 - conf) / 2
    naive_ci <- unname(stats::quantile(v, c(a, 1 - a), na.rm = TRUE))
    draws$naive_r <- v
  }

  ## ---- calibration ---------------------------------------------------------
  ## Regress target on score. Perfect calibration is intercept 0, slope 1.
  cal <- tryCatch({
    fit <- stats::lm(target ~ score)
    cf  <- stats::coef(fit)
    data.frame(intercept = unname(cf[1]), slope = unname(cf[2]),
               row.names = NULL)
  }, error = function(e) data.frame(intercept = NA_real_, slope = NA_real_))

  cal_ci <- NULL
  if (n_boot > 0 && k > 2L) {
    Bc <- matrix(NA_real_, n_boot, 2, dimnames = list(NULL, c("intercept","slope")))
    for (b in seq_len(n_boot)) {
      pick <- sample(ug, k, replace = TRUE)
      idx  <- unlist(lapply(pick, function(g) which(cluster == g)), use.names = FALSE)
      cf <- tryCatch(stats::coef(stats::lm(target[idx] ~ score[idx])),
                     error = function(e) c(NA, NA))
      Bc[b, ] <- cf[1:2]
    }
    a <- (1 - conf) / 2
    cal_ci <- t(apply(Bc, 2, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE))
    colnames(cal_ci) <- c("lower", "upper")
  }

  ## ---- cluster permutation null -------------------------------------------
  ## Permute whole trajectories, not rows: shuffling rows would destroy the
  ## dependence structure and give an anticonservative null.
  null_test <- NULL
  if (n_perm > 0 && k > 2L) {
    obs <- point[["r"]]
    null_r <- rep(NA_real_, n_perm)
    cnt <- 0L; tot <- 0L
    for (p in seq_len(n_perm)) {
      remap <- stats::setNames(sample(ug), ug)
      ## Build the permuted target by swapping trajectory blocks wholesale.
      perm_target <- target
      for (g in ug) {
        src <- which(cluster == remap[[g]])
        dst <- which(cluster == g)
        L <- min(length(src), length(dst))
        perm_target[dst[seq_len(L)]] <- target[src[seq_len(L)]]
      }
      s <- .stats(score, perm_target)
      if (!is.null(s) && is.finite(s[["r"]])) {
        tot <- tot + 1L
        null_r[p] <- s[["r"]]
        if (abs(s[["r"]]) >= abs(obs)) cnt <- cnt + 1L
      }
    }
    null_test <- list(
      statistic = "Pearson r",
      observed  = obs,
      p_value   = (cnt + 1) / (tot + 1),
      n_perm    = tot,
      method    = "cluster permutation (whole trajectories exchanged)"
    )
    draws$null_r <- null_r
  }

  ## ---- assemble ------------------------------------------------------------
  agreement <- data.frame(
    statistic = c("pearson_r", "lins_ccc", "rmse", "mae", "bias", "r2"),
    row_level = unname(point[c("r","ccc","rmse","mae","bias","r2")]),
    stringsAsFactors = FALSE
  )
  if (!is.null(cluster_point)) {
    agreement$cluster_mean_level <-
      unname(cluster_point[c("r","ccc","rmse","mae","bias","r2")])
  }
  if (!is.null(boot_ci)) {
    m <- boot_ci[c("r","ccc","rmse","mae","bias","r2"), , drop = FALSE]
    agreement$ci_lower <- unname(m[, "lower"])
    agreement$ci_upper <- unname(m[, "upper"])
  }

  mse <- point[["rmse"]]^2
  parts <- unname(point[c("sb","sdsd","lcs")])
  ## Guard the degenerate case of exact reproduction, where MSE is 0 and the
  ## percentage split is 0/0 rather than undefined in any interesting sense.
  pct <- if (is.finite(mse) && mse > 0) 100 * parts / mse else rep(0, 3)
  decomposition <- data.frame(
    component = c("squared_bias", "variance_mismatch", "lack_of_correlation"),
    value     = parts,
    percent   = pct,
    meaning   = c("systematic offset; removable by recentring",
                  "score too flat or too volatile relative to target",
                  "score does not track the target's pattern; rescaling cannot fix"),
    stringsAsFactors = FALSE
  )
  attr(decomposition, "mse") <- mse
  attr(decomposition, "residual") <- mse - sum(parts)   # should be ~0

  dependence <- data.frame(
    n_observations = n,
    n_clusters     = k,
    mean_cluster_size = mean(sizes),
    icc            = icc,
    design_effect  = deff,
    effective_n    = n_eff
  )
  if (!is.null(naive_ci) && !is.null(boot_ci)) {
    dependence$r_ci_naive_width <-
      naive_ci[2] - naive_ci[1]
    dependence$r_ci_cluster_width <-
      boot_ci["r", "upper"] - boot_ci["r", "lower"]
    attr(dependence, "r_ci_naive") <- naive_ci
  }

  if (!is.null(cal_ci)) {
    cal$intercept_lower <- cal_ci["intercept","lower"]
    cal$intercept_upper <- cal_ci["intercept","upper"]
    cal$slope_lower     <- cal_ci["slope","lower"]
    cal$slope_upper     <- cal_ci["slope","upper"]
  }

  ## ---- interpretation notes -----------------------------------------------
  if (is.finite(deff) && deff > 2) {
    notes <- c(notes, sprintf(
      "Rows are strongly clustered (ICC %.2f, design effect %.1f). The %d observations carry roughly the information of %.0f independent ones; quote the cluster bootstrap interval, not one based on n = %d.",
      icc, deff, n, n_eff, n))
  }
  if (!is.null(naive_ci) && !is.null(boot_ci)) {
    wn <- naive_ci[2] - naive_ci[1]
    wc <- boot_ci["r","upper"] - boot_ci["r","lower"]
    if (is.finite(wn) && is.finite(wc) && wn > 0 && wc / wn > 1.5) {
      notes <- c(notes, sprintf(
        "Ignoring clustering would give a %.0f%% interval for r of width %.3f; resampling whole trajectories gives width %.3f, %.1f times wider. Report the latter.",
        100 * conf, wn, wc, wc / wn))
    }
  }
  rr <- point[["r"]]; cc <- point[["ccc"]]
  if (is.finite(rr) && is.finite(cc) && (rr - cc) > 0.1) {
    notes <- c(notes, sprintf(
      "Correlation (%.3f) exceeds concordance (%.3f). The score tracks the target's pattern but does not agree with it in level or scale; see calibration.",
      rr, cc))
  }
  if (is.finite(cal$slope) && abs(cal$slope - 1) > 0.2) {
    notes <- c(notes, sprintf(
      "Calibration slope is %.2f rather than 1: the score compresses or exaggerates the target's range.",
      cal$slope))
  }
  if (is.finite(point[["r2"]]) && point[["r2"]] < 0) {
    notes <- c(notes,
      "R2 is negative: as an absolute predictor the score does worse than the target's own mean. It may still rank correctly; check pearson_r.")
  }
  big <- decomposition$component[which.max(decomposition$percent)]
  notes <- c(notes, sprintf("Error is dominated by %s (%.0f%% of MSE).",
                            big, max(decomposition$percent)))
  notes <- c(notes,
    "If score and target derive from the same generator, this is internal consistency, not validation.")

  out <- list(agreement     = agreement,
              calibration   = cal,
              decomposition = decomposition,
              dependence    = dependence,
              null_test     = null_test,
              notes         = notes,
              data          = data.frame(score = score, target = target,
                                         cluster = cluster,
                                         stringsAsFactors = FALSE),
              draws         = draws,
              conf          = conf)
  class(out) <- c("rri_accuracy", "list")
  out
}

#' @describeIn rri_accuracy Compact console summary of an accuracy assessment.
#' @param x An `rri_accuracy` object.
#' @param ... Ignored.
#' @export
print.rri_accuracy <- function(x, ...) {
  cat("Agreement with reference target\n")
  cat(strrep("-", 62), "\n")
  ag <- x$agreement
  num <- vapply(ag, is.numeric, logical(1))
  ag[num] <- lapply(ag[num], function(z) round(z, 4))
  print(ag, row.names = FALSE)

  d <- x$dependence
  cat("\nDependence structure\n")
  cat(sprintf("  %d observations in %d clusters (mean size %.1f)\n",
              d$n_observations, d$n_clusters, d$mean_cluster_size))
  if (is.finite(d$icc)) {
    cat(sprintf("  ICC %.3f, design effect %.1f, effective n %.0f\n",
                d$icc, d$design_effect, d$effective_n))
  }

  cat("\nError decomposition (percent of MSE)\n")
  dc <- x$decomposition
  for (i in seq_len(nrow(dc))) {
    cat(sprintf("  %-20s %5.1f%%\n", dc$component[i], dc$percent[i]))
  }

  cat(sprintf("\nCalibration: target = %.3f + %.3f x score  (ideal 0 and 1)\n",
              x$calibration$intercept, x$calibration$slope))

  if (!is.null(x$null_test)) {
    cat(sprintf("Cluster permutation test: p = %.4f (%d permutations)\n",
                x$null_test$p_value, x$null_test$n_perm))
  }

  cat("\nNotes\n")
  for (nn in x$notes) {
    cat(paste(strwrap(nn, width = 76, prefix = "  "), collapse = "\n"), "\n")
  }
  invisible(x)
}
