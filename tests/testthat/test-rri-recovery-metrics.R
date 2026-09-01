library(testthat)

## Companion to test-recovery-metrics.R. That file covers the metric maths on
## hand-built trajectories; this one checks the function against simulator
## output scored through a gene-augmented microbial block, which is the path
## users take in the vignette.

test_that("rri_recovery_metrics returns its documented column contract", {
  sim <- simulate_redox_holobiont(
    n_plot = 2,
    n_depth = 1,
    n_plant = 2,
    n_time = 12,
    p_micro = 20,
    seed = 1
  )

  ## Traits plus log-scaled gene abundances: two microbial views in one block.
  micro_block <- cbind(sim$micro_traits, log1p(sim$micro_gene_abundance))

  res <- suppressWarnings(rri_pipeline_st(
    ROS_flux = sim$ROS_flux,
    Eh_stability = sim$Eh_stability,
    micro_data = micro_block,
    id = sim$id,
    reducer = "per_domain",
    scaling = "pnorm"
  ))

  rec <- rri_recovery_metrics(
    res = res,
    id = sim$id,
    time_col = "time",
    group_cols = c("plot", "depth", "plant_id"),
    perturb_start = 5,
    perturb_end = 7
  )

  ## Grouping keys, one row per trajectory.
  expect_true(all(c("plot", "depth", "plant_id") %in% names(rec)))
  expect_equal(nrow(rec), 4)

  ## The six recovery signatures, plus their normalised aliases.
  signature_cols <- c(
    "depth_min", "depth_min_frac",     # pathway truncation
    "tau_lag",                         # response lag
    "overshoot", "overshoot_frac",     # overshoot
    "H_hysteresis", "H", "H_abs", "H_norm",
    "incomplete_return", "incomplete_return_frac", "I", "I_norm",
    "displaced_plateau_flag", "displaced_plateau_level"
  )
  expect_true(all(signature_cols %in% names(rec)))

  ## Fit diagnostics must accompany every rate estimate, so a reported k can
  ## always be traced to its fit quality and sample size.
  fit_cols <- c("k", "k_recovery", "t_half",
                "fit_status", "fit_r_squared", "fit_start_time", "n_fit")
  expect_true(all(fit_cols %in% names(rec)))
  expect_false(any(is.na(rec$fit_status)))

  ## Observation counts are always reported, never NA.
  count_cols <- c("n_pre", "n_perturb", "n_recovery", "n_missing")
  expect_true(all(count_cols %in% names(rec)))
  expect_false(any(is.na(rec[, count_cols])))

  ## Depth of decline is finite for at least one trajectory and non-negative
  ## wherever it is defined.
  expect_true(any(is.finite(rec$depth_min)))
  expect_true(all(rec$depth_min[is.finite(rec$depth_min)] >= 0))

  ## Normalised quantities stay bounded.
  for (nm in c("H_norm", "I_norm")) {
    v <- rec[[nm]][is.finite(rec[[nm]])]
    if (length(v)) expect_true(all(v >= 0 & v <= 1))
  }

  ## Legacy alt_routing fields are retained as NA by design; displaced_plateau
  ## is the supported successor.
  expect_true(all(is.na(rec$alt_routing_flag)))
  expect_true(all(is.na(rec$alt_routing_level)))

  ## Hysteresis is not evaluated without a measured forcing column.
  expect_true(all(rec$H_axis == "unavailable"))
  expect_true(all(rec$hysteresis_status == "not_evaluated"))
})
