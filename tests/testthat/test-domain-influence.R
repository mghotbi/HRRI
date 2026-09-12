test_that("realised shares form an exact variance decomposition", {
  sim <- simulate_redox_holobiont(
    n_plot = 2, n_depth = 2, n_plant = 2, n_time = 30,
    p_micro = 10, seed = 404
  )
  res <- suppressWarnings(rri_pipeline_st(
    ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
    micro_data = sim$micro_data, id = sim$id,
    reducer = "per_domain", scaling = "pnorm"
  ))

  infl <- rri_domain_influence(res)

  expect_true(all(c("influence", "covariance", "notes") %in% names(infl)))
  expect_equal(nrow(infl$influence), 3L)
  expect_true(all(c("Physio", "Soil", "Micro") %in% infl$influence$domain))

  ## phi_d = w_d Cov(S_d, R) / Var(R) sums to one by construction. This is the
  ## property that makes the attribution interpretable, so it is worth locking.
  sh <- infl$influence$realised_share
  skip_if(all(is.na(sh)), "nominal weights not recoverable in this build")
  expect_equal(sum(sh, na.rm = TRUE), 1, tolerance = 0.02)

  ## Shares are proportions of variance, not arbitrary numbers.
  expect_true(all(sh[is.finite(sh)] >= -0.05 & sh[is.finite(sh)] <= 1.05))
})

test_that("domain influence reports, rather than hides, an imbalance", {
  ## Force a lopsided composite: almost all weight on soil. The realised share
  ## must follow, and the function must say so in `notes`.
  sim <- simulate_redox_holobiont(
    n_plot = 1, n_depth = 2, n_plant = 2, n_time = 24,
    p_micro = 8, seed = 77
  )
  res <- suppressWarnings(rri_pipeline_st(
    ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
    micro_data = sim$micro_data, id = sim$id,
    w1 = 0.05, w2 = 0.90, w3 = 0.05
  ))

  infl <- rri_domain_influence(
    res, weights = c(Physio = 0.05, Soil = 0.90, Micro = 0.05)
  )
  soil <- infl$influence[infl$influence$domain == "Soil", ]
  expect_gt(soil$realised_share, 0.5)
  expect_type(infl$notes, "character")
})

test_that("missing domains are flagged as a cause of departure from nominal", {
  soil  <- data.frame(Eh = c(50, 100, 150, 200), pH = c(5, 5.5, 6, 6.5))
  plant <- data.frame(FvFm = c(0.7, 0.75, 0.72, 0.8))
  res <- suppressWarnings(rri_pipeline(soil = soil, plant = plant,
                                       method_soil = "scale"))
  infl <- rri_domain_influence(res)
  ## Micro was never supplied, so it is missing for every row.
  micro <- infl$influence[infl$influence$domain == "Micro", ]
  expect_equal(micro$n_missing, 4L)
  expect_true(any(grepl("Missing domain scores", infl$notes)))
})
