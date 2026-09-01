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
