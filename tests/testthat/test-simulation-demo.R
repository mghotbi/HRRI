test_that("simulation demonstration uses the current simulator contract", {
  d <- suppressWarnings(rri_simulation_demo(seed = 701L))
  expect_true(all(c("flood_drain", "drought_rewet") %in% names(d$scenarios)))
  expect_true(all(c("scores", "observations", "capacity", "metadata") %in% names(d)))
  expect_true(nrow(d$scores) > 0L)
  expect_true(any(is.finite(d$scores$RRI)))
  expect_true(all(is.finite(d$capacity$accessible_EAC)))
  expect_true(all(is.finite(d$capacity$accessible_EDC)))
})
