test_that("conditional recovery curve recovers its own estimands", {
  t <- 0:20
  y <- ifelse(t <= 2, 10, 10*(1-.6*exp(-.25*pmax(t-3,0))))
  x <- data.frame(plot="P1",depth="D1",time=t,EAC=y)
  a <- hrri_infer_architecture(x,perturb_time=3,baseline_end=2,verbose=FALSE)
  expect_equal(a$estimates$deficit_fraction_eac,.6,tolerance=.03)
  expect_equal(a$estimates$recovery_rate_eac,.25,tolerance=.03)
  expect_match(a$interpretation,"not Q-alpha-k-memory")
})
test_that("mechanistic alpha/k are not accepted as matching estimands", {
  x <- data.frame(plot="P1",depth="D1",time=0:10,
    EAC=c(rep(10,3),seq(4,9,length.out=8)),alpha_accept=.5,k_accept_h=.1)
  a <- hrri_infer_architecture(x,perturb_time=3,baseline_end=2,verbose=FALSE)
  expect_error(validate_architecture(a,"alpha_eac:alpha_accept"),"cannot validate")
})
