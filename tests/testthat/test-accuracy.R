## tests/testthat/test-accuracy.R
## Properties that must hold regardless of the data, so a regression in the
## arithmetic is caught rather than eyeballed.

## Inside the package (devtools::test(), R CMD check) rri_accuracy() is already
## on the search path and this block does nothing. Running this file on its own
## from RStudio's "Run Tests" button loads no package, so every test would fail
## with "could not find function" -- the guard below sources the definition in
## that case only.
for (fn in c("rri_accuracy", "plot_rri_accuracy", "theme_ems")) {
  if (exists(fn, mode = "function")) next
  cand <- file.path(c(".", "../../R", "R"), paste0(fn, ".R"))
  hit  <- cand[file.exists(cand)]
  if (length(hit)) {
    source(hit[1])
  } else if (fn == "rri_accuracy") {
    stop("rri_accuracy() not found. Either run devtools::load_all() first, ",
         "or put rri_accuracy.R beside this file.")
  }
}

test_that("perfect agreement gives r = ccc = 1 and zero error", {
  set.seed(1)
  y <- runif(200)
  g <- rep(1:20, each = 10)
  a <- rri_accuracy(y, y, cluster = g, n_boot = 0, n_perm = 0)
  ag <- stats::setNames(a$agreement$row_level, a$agreement$statistic)
  expect_equal(ag[["pearson_r"]], 1, tolerance = 1e-8)
  expect_equal(ag[["lins_ccc"]],  1, tolerance = 1e-8)
  expect_equal(ag[["rmse"]], 0, tolerance = 1e-8)
  expect_equal(ag[["bias"]], 0, tolerance = 1e-8)
})

test_that("MSE decomposition is exact", {
  set.seed(2)
  y <- runif(300)
  x <- 0.8 * y + 0.05 + rnorm(300, 0, 0.05)
  g <- rep(1:15, each = 20)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 0)
  expect_equal(attr(a$decomposition, "residual"), 0, tolerance = 1e-10)
  expect_equal(sum(a$decomposition$percent), 100, tolerance = 1e-6)
})

test_that("ccc detects miscalibration that r does not", {
  set.seed(3)
  y <- runif(200)
  x <- 2 * y + 0.3                     # perfectly correlated, badly calibrated
  g <- rep(1:20, each = 10)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 0)
  ag <- stats::setNames(a$agreement$row_level, a$agreement$statistic)
  expect_equal(ag[["pearson_r"]], 1, tolerance = 1e-8)
  expect_lt(ag[["lins_ccc"]], 0.7)     # agreement is poor despite r = 1
  expect_lt(a$calibration$slope, 0.6)  # target ~ score has slope 1/2
})

test_that("clustering widens the interval relative to the naive one", {
  set.seed(4)
  k <- 12; m <- 40
  u <- rnorm(k, 0, 0.3)
  y <- unlist(lapply(u, function(uu) uu + 0.5 + rnorm(m, 0, 0.05)))
  x <- 0.9 * y + rnorm(k * m, 0, 0.08)
  g <- rep(seq_len(k), each = m)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 400, n_perm = 0, seed = 4)
  expect_gt(a$dependence$icc, 0.1)
  expect_gt(a$dependence$design_effect, 2)
  expect_lt(a$dependence$effective_n, a$dependence$n_observations)
  expect_gt(a$dependence$r_ci_cluster_width, a$dependence$r_ci_naive_width)
})

test_that("independent rows give design effect near 1", {
  set.seed(5)
  y <- runif(300); x <- y + rnorm(300, 0, 0.1)
  g <- rep(1:30, each = 10)            # groups carry no signal
  a <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 0)
  expect_lt(a$dependence$design_effect, 2)
})

test_that("permutation p-value is calibrated under the null", {
  set.seed(6)
  k <- 10; m <- 20
  y <- unlist(lapply(rnorm(k), function(uu) uu + rnorm(m, 0, 0.1)))
  x <- unlist(lapply(rnorm(k), function(uu) uu + rnorm(m, 0, 0.1)))
  g <- rep(seq_len(k), each = m)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 300, seed = 6)
  expect_true(a$null_test$p_value >= 0 && a$null_test$p_value <= 1)
})

test_that("missing cluster argument warns rather than failing silently", {
  set.seed(7)
  y <- runif(50); x <- y + rnorm(50, 0, 0.1)
  expect_warning(rri_accuracy(x, y, n_boot = 0, n_perm = 0), "cluster")
})

test_that("length mismatch is an error", {
  expect_error(rri_accuracy(1:10, 1:9), "same length")
})

test_that("non-finite rows are dropped, not propagated", {
  set.seed(8)
  y <- runif(100); x <- y + rnorm(100, 0, 0.05)
  x[c(3, 9)] <- NA; y[50] <- NA
  g <- rep(1:10, each = 10)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 0)
  expect_equal(a$dependence$n_observations, 97)
  expect_true(all(is.finite(a$agreement$row_level)))
})

test_that("print method runs and returns its input invisibly", {
  set.seed(9)
  y <- runif(120); x <- 0.8 * y + 0.05 + rnorm(120, 0, 0.05)
  g <- rep(1:12, each = 10)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 50, n_perm = 50, seed = 9)
  expect_output(print(a), "Agreement with reference target")
  expect_identical(withVisible(print(a))$visible, FALSE)
})

test_that("resampled draws are retained for plotting", {
  set.seed(10)
  y <- runif(120); x <- 0.8 * y + rnorm(120, 0, 0.05)
  g <- rep(1:12, each = 10)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 60, n_perm = 60, seed = 10)
  expect_true(all(c("cluster", "naive_r", "null_r") %in% names(a$draws)))
  expect_equal(nrow(a$draws$cluster), 60L)
  expect_equal(length(a$draws$naive_r), 60L)
  expect_equal(nrow(a$data), 120L)
})

test_that("plot_rri_accuracy builds without error", {
  skip_if_not_installed("ggplot2")
  set.seed(11)
  k <- 6; m <- 15
  y <- unlist(lapply(rnorm(k, 0, 0.3), function(u) u + 0.5 + rnorm(m, 0, .05)))
  x <- 0.8 * y + 0.05 + rnorm(k * m, 0, 0.05)
  g <- rep(seq_len(k), each = m)
  a <- rri_accuracy(x, y, cluster = g, n_boot = 60, n_perm = 60, seed = 11)

  ## Single panel returns a bare ggplot regardless of patchwork.
  p1 <- plot_rri_accuracy(a, panels = "calibration")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p1), "ggplot_built")

  ## Every panel builds.
  for (nm in c("calibration", "agreement", "precision", "error")) {
    built <- ggplot2::ggplot_build(plot_rri_accuracy(a, panels = nm))
    expect_s3_class(built, "ggplot_built")
  }
})

test_that("plot_rri_accuracy rejects the wrong input and a missing bootstrap", {
  expect_error(plot_rri_accuracy(list(a = 1)), "rri_accuracy")

  skip_if_not_installed("ggplot2")
  set.seed(12)
  y <- runif(60); x <- y + rnorm(60, 0, 0.05)
  g <- rep(1:6, each = 10)
  a0 <- rri_accuracy(x, y, cluster = g, n_boot = 0, n_perm = 0)
  ## No draws: the precision panel must degrade to a placeholder, not fail.
  expect_s3_class(
    ggplot2::ggplot_build(plot_rri_accuracy(a0, panels = "precision")),
    "ggplot_built")
})

test_that("the ggtern guard refuses to load an incompatible ggtern", {
  skip_if_not_installed("ggplot2")
  guard <- tryCatch(get("hrri_ggtern_ok", envir = asNamespace("HRRI")),
                    error = function(e) NULL)
  skip_if(is.null(guard), "internal helper not reachable outside the package")

  ## Whether ggtern was already loaded by something else is not this test's
  ## business; only the delta across the call is.
  before <- "ggtern" %in% loadedNamespaces()
  ok <- guard()
  after  <- "ggtern" %in% loadedNamespaces()

  expect_type(ok, "logical")
  if (utils::packageVersion("ggplot2") >= "4.0.0") {
    expect_false(ok)
    ## The point of the guard: on an incompatible ggplot2 it must not be the
    ## thing that loads ggtern.
    expect_equal(after, before)
  }
})
