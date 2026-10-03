test_that("all six paper figures have documented native R entry points", {
  funcs <- c("plot_rri_framework","plot_rri_identifiability","plot_rri_timeseries",
             "plot_rri_recovery_diagnostics","plot_rri_properties","plot_rri_accuracy")
  expect_true(all(funcs %in% getNamespaceExports("HRRI")))
  p <- plot_rri_framework()
  b <- ggplot2::ggplot_build(p)
  seg <- which(vapply(p$layers,function(x) inherits(x$geom,"GeomSegment"),logical(1)))
  expect_true(all(vapply(b$data[seg],nrow,integer(1))==1L))
})

test_that("analytical panels retain their stated equations", {
  p <- plot_rri_identifiability()
  d <- p[[1]]$data
  expected <- 100*.4*(-expm1(-.15*10))
  expect_equal(100*d$alpha*(-expm1(-d$k*10)),rep(expected,nrow(d)),tolerance=1e-10)
  cond <- p[[4]]$data
  expect_equal(tail(cond$sensitivity,1),1/abs(6/expm1(6)-12/expm1(12)),tolerance=1e-10)
  expect_error(plot_rri_identifiability(alpha=1),"strictly")
  expect_error(plot_rri_identifiability(duration_pairs=matrix(c(3,2),nrow=1)),"increasing")
})

test_that("availability counts actual finite values and rejects unmatched trajectories", {
  id <- data.frame(plot=rep(c("A","B"),each=3),time=rep(1:3,2))
  res <- list(row_scores=data.frame(RRI=c(.5,.2,.4,.6,.3,.5)))
  rec <- data.frame(plot=c("A","B"),k_recovery=c(.1,NA_real_),H_hysteresis=NA_real_)
  p <- plot_rri_recovery_diagnostics(res,id,rec,group_cols="plot")
  d <- p[[2]]$data
  expect_equal(d$n[d$metric=="k_recovery"],1L)
  expect_equal(d$n[d$metric=="H_hysteresis"],0L)
  expect_false(d$present[d$metric=="tau_lag"])
  expect_error(plot_rri_recovery_diagnostics(res,id,rec[1,,drop=FALSE],group_cols="plot"),"exactly one")
})

test_that("paper layout preserves outliers and degenerate bootstrap draws", {
  target <- seq(.1,.9,length.out=40); score <- target+.1;score[40]<-3
  a <- rri_accuracy(score,target,cluster=rep(1:4,each=10),n_boot=1,n_perm=0,seed=9)
  p <- plot_rri_accuracy(a,panels="calibration",style="paper")
  expect_equal(max(ggplot2::ggplot_build(p)$data[[1]]$y),3)
  expect_no_warning(ggplot2::ggplot_build(plot_rri_accuracy(a,panels="precision")))
})
