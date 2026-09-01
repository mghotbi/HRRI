test_that("benchmark is reproducible and reports failures explicitly", {
  args <- list(n_plot=1,n_depth=1,n_plant=1,n_time=10,p_micro=5)
  b1 <- suppressWarnings(benchmark_hrri(n=1,seed_start=31,sim_args=args,
    domains="soil",verbose=FALSE,
    pipeline_args=list(method_soil="scale",direction_anchor_soil="Eh")))
  b2 <- suppressWarnings(benchmark_hrri(n=1,seed_start=31,sim_args=args,
    domains="soil",verbose=FALSE,
    pipeline_args=list(method_soil="scale",direction_anchor_soil="Eh")))
  expect_identical(b1$row_data,b2$row_data)
  expect_true(all(c("RMSE","Bias","r_truth","rank_truth","n_failed") %in%
                    names(b1$summary)))
  expect_false("Coverage95" %in% names(b1$summary))
  expect_equal(b1$settings$missingness,"MCAR")
})
test_that("seed cannot be duplicated through sim_args", {
  expect_error(benchmark_hrri(n=1,sim_args=list(seed=1),verbose=FALSE),
               "seed_start")
})
