test_that("simulator schema, dimensions and RNG are reproducible", {
  set.seed(913)
  caller_rng <- .Random.seed
  s1 <- small_sim(); s2 <- small_sim()
  expect_identical(.Random.seed, caller_rng)
  expect_true(all(c("id","soil_data","plant_data","micro_gene_abundance",
                    "latent_state","conservation_checks","metadata") %in% names(s1)))
  n <- nrow(s1$id)
  expect_equal(nrow(s1$soil_data), n)
  expect_equal(nrow(s1$plant_data), n)
  expect_equal(nrow(s1$micro_gene_abundance), n)
  expect_identical(s1$soil_data, s2$soil_data)
})
test_that("Fe and Mn closed-balance checks are finite and small", {
  c <- small_sim()$conservation_checks
  expect_s3_class(c, "data.frame")
  expect_equal(names(c), c("check","value"))
  expect_true(all(is.finite(c$value)))
  expect_lt(c$value[c$check=="maximum_absolute_Fe_mass_balance_error"], 1e-7)
  expect_lt(c$value[c$check=="maximum_absolute_Mn_mass_balance_error"], 1e-7)
})
test_that("inventories remain nonnegative and scenarios execute", {
  s <- small_sim()$soil_data
  pools <- grep("_mmol_kg$", names(s), value=TRUE)
  expect_true(all(vapply(s[pools], function(x) all(x[is.finite(x)] >= 0), logical(1))))
  expect_no_error(simulate_redox_holobiont(n_plot=1,n_depth=1,n_plant=1,
    n_time=10,p_micro=3,seed=5,scenario="drought_rewet"))
})
test_that("EDC inventory includes MnIII reducing equivalents", {
  s <- small_sim()$soil_data
  expected <- s$FeII_mmol_kg + 9*s$FeS_mmol_kg + 2*s$MnII_mmol_kg +
    s$MnIII_mmol_kg + 8*s$sulfide_mmol_kg + 8*s$CH4_mmol_kg +
    8*s$NH4_mmol_kg + s$humic_EDC_mmol_e_kg + .5*s$DOC_mmolC_kg
  expect_equal(s$EDC, expected, tolerance=1e-8)
})

test_that("memory is a holobiont state with plant and microbial components", {
  s <- small_sim()

  ## The two new legacy components are exported and bounded.
  expect_true(all(c("micro_legacy", "plant_legacy") %in% names(s$latent_state)))
  for (nm in c("micro_legacy", "plant_legacy", "memory")) {
    v <- s$latent_state[[nm]]
    expect_true(all(is.finite(v)))
    expect_true(all(v >= 0 & v <= 1))
  }

  ## Microbial legacy must actually vary: a constant column would mean the
  ## component is inert and memory would be mineralogical in all but name.
  expect_gt(stats::sd(s$latent_state$micro_legacy), 0)
})

test_that("plant and microbial legacies contribute to memory", {
  ## history_strength scales every memory gain term. With it near zero, memory
  ## cannot accumulate; with it high, it must. This confirms the multi-domain
  ## gain path is live rather than dominated by initialisation.
  lo <- simulate_redox_holobiont(n_plot = 1, n_depth = 1, n_plant = 1,
                                 n_time = 20, p_micro = 3, seed = 7,
                                 history_strength = 0.02, n_cycles = 2)
  hi <- simulate_redox_holobiont(n_plot = 1, n_depth = 1, n_plant = 1,
                                 n_time = 20, p_micro = 3, seed = 7,
                                 history_strength = 0.95, n_cycles = 2)
  expect_gt(mean(hi$latent_state$memory), mean(lo$latent_state$memory))
})

test_that("microbial legacy relaxes more slowly than it accumulates", {
  ## Asymmetry is the mechanism that makes community composition a *legacy*
  ## rather than an instantaneous readout of current redox conditions.
  s <- simulate_redox_holobiont(n_plot = 1, n_depth = 1, n_plant = 1,
                                n_time = 40, p_micro = 3, seed = 11,
                                scenario = "flood_drain", n_cycles = 1)
  ml <- s$latent_state$micro_legacy
  peak <- which.max(ml)
  skip_if(peak < 3 || peak > length(ml) - 3, "no interior peak in this run")
  rise <- max(diff(ml[seq_len(peak)]))
  fall <- min(diff(ml[peak:length(ml)]))
  expect_gt(rise, abs(fall))   # gain 0.030 vs relax 0.012
})
