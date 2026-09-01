make_trajectory <- function(depth = .35) {
  t <- 1:30
  y <- c(rep(.7,7), seq(.7,.7-depth,length.out=11),
         seq(.7-depth,.68,length.out=8), rep(.68,4))
  data.frame(plot="P1", depth="D1", time=t, RRI=y)
}
metric <- function(x) rri_recovery_metrics(x, time_col="time",
  group_cols=c("plot","depth"), perturb_start=8, perturb_end=18)
test_that("recovery metrics run on aligned trajectories", {
  z <- metric(make_trajectory())
  expect_s3_class(z, "data.frame")
  expect_equal(nrow(z), 1L)
  expect_true(is.finite(z$baseline_rri))
  expect_true(z$depth_min_frac >= 0)
  expect_true(z$n_pre >= 3L)
  expect_true(z$n_recovery >= 3L)
})
test_that("deeper observed decline gives larger depth diagnostic", {
  expect_gt(metric(make_trajectory(.4))$depth_min_frac,
            metric(make_trajectory(.1))$depth_min_frac)
})
test_that("duplicate trajectory time is rejected", {
  x <- rbind(make_trajectory(), make_trajectory()[1,])
  expect_error(metric(x), "Duplicate times")
})
