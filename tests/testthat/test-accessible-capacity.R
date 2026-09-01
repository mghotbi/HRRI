## tests/testthat/test-accessible-capacity.R
## Tests for rri_default_reservoirs() and rri_accessible_capacity()

# ── rri_default_reservoirs() ──────────────────────────────────────────────────

test_that("rri_default_reservoirs returns a valid named list", {
  res <- rri_default_reservoirs()
  expect_type(res, "list")
  expect_gt(length(res), 0)
  for (nm in names(res)) {
    r <- res[[nm]]
    expect_true("Q_col" %in% names(r), label = paste(nm, "has Q_col"))
    expect_true("alpha" %in% names(r), label = paste(nm, "has alpha"))
    expect_true("k"     %in% names(r), label = paste(nm, "has k"))
    expect_true("type"  %in% names(r), label = paste(nm, "has type"))
    expect_true(r$type %in% c("EAC", "EDC"),
                label = paste(nm, "type is EAC or EDC"))
  }
})

test_that("rri_default_reservoirs contains at least one EAC and one EDC reservoir", {
  res   <- rri_default_reservoirs()
  types <- vapply(res, `[[`, character(1), "type")
  expect_true("EAC" %in% types, label = "at least one EAC reservoir")
  expect_true("EDC" %in% types, label = "at least one EDC reservoir")
})

test_that("default reservoirs operate on the core EAC/EDC schema", {
  x <- data.frame(EAC = c(10, 20), EDC = c(8, 12))
  z <- rri_accessible_capacity(x, rri_default_reservoirs(), tau = 24,
                               normalise = FALSE)
  expect_true(all(is.finite(z$cacc_eac)))
  expect_true(all(is.finite(z$cacc_edc)))
  expect_equal(z$n_reservoirs_observed, c(2L, 2L))
})

# ── rri_accessible_capacity() ─────────────────────────────────────────────────

test_that("rri_accessible_capacity returns all required list slots", {
  sim <- small_sim()
  res <- rri_accessible_capacity(sim$soil_data, sim_reservoirs(), tau = 24)
  for (s in c("cacc", "cacc_raw", "cacc_eac", "cacc_edc",
              "net_oxidative_balance", "cacc_fraction", "total_inventory")) {
    expect_true(s %in% names(res), label = paste("missing slot:", s))
  }
})

test_that("cacc_raw is finite and non-negative", {
  sim  <- small_sim()
  res  <- rri_accessible_capacity(sim$soil_data, sim_reservoirs(), tau = 24)
  vals <- res$cacc_raw[!is.na(res$cacc_raw)]
  expect_true(all(is.finite(vals)), label = "cacc_raw finite")
  expect_true(all(vals >= 0),       label = "cacc_raw non-negative")
})

test_that("normalised cacc is in [0, 1]", {
  sim  <- small_sim()
  res  <- rri_accessible_capacity(sim$soil_data, sim_reservoirs(),
                                   tau = 24, normalise = TRUE)
  vals <- res$cacc[!is.na(res$cacc)]
  expect_true(all(vals >= 0 & vals <= 1), label = "cacc in [0,1]")
})

test_that("raw cacc equals formula value for a single scalar reservoir", {
  Q <- 80; alpha <- 0.75; k <- 1.2; tau <- 24
  df       <- data.frame(Q_col = Q, alpha_col = alpha, k_col = k)
  res_spec <- list(r1 = list(Q_col="Q_col", alpha="alpha_col",
                              k="k_col", type="EAC"))
  res      <- rri_accessible_capacity(df, res_spec, tau = tau, normalise = FALSE)
  expected <- Q * alpha * (1 - exp(-k * tau))
  expect_equal(res$cacc_raw, expected, tolerance = 1e-6)
})

test_that("cacc increases monotonically with tau", {
  df       <- data.frame(Q = 100, a = 0.8, k = 0.5)
  res_spec <- list(r1 = list(Q_col="Q", alpha="a", k="k", type="EAC"))
  tau_vals <- c(1, 6, 24, 72, 168)
  vals     <- sapply(tau_vals, function(tt)
    rri_accessible_capacity(df, res_spec, tau=tt, normalise=FALSE)$cacc_raw)
  expect_true(all(diff(vals) >= -1e-9), label = "cacc non-decreasing in tau")
})

test_that("cacc_eac and cacc_edc sum correctly to cacc_raw", {
  sim <- small_sim()
  res <- rri_accessible_capacity(sim$soil_data, sim_reservoirs(),
                                  tau = 24, normalise = FALSE)
  total <- res$cacc_eac + res$cacc_edc
  expect_equal(total, res$cacc_raw, tolerance = 1e-6)
})

test_that("return_components = TRUE attaches a components data frame", {
  sim <- small_sim()
  res <- rri_accessible_capacity(sim$soil_data, sim_reservoirs(),
                                  tau = 24, return_components = TRUE)
  expect_true("components" %in% names(res))
  comp <- res$components
  expect_s3_class(comp, "data.frame")
  expect_true("reservoir"         %in% names(comp))
  expect_true("mean_contribution" %in% names(comp))
  expect_equal(nrow(comp), length(sim_reservoirs()))
})

 

test_that("very large tau saturates cacc_raw at Q * alpha", {
  Q <- 100; alpha <- 0.75; k <- 2.0; tau <- 1e6
  df       <- data.frame(Q_col = Q, alpha_col = alpha, k_col = k)
  res_spec <- list(r1 = list(Q_col="Q_col", alpha="alpha_col",
                              k="k_col", type="EAC"))
  res      <- rri_accessible_capacity(df, res_spec, tau=tau, normalise=FALSE)
  expect_equal(res$cacc_raw, Q * alpha, tolerance = 1e-4)
})
