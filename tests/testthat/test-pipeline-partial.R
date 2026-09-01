test_that("partial domains remain explicit and weights are forwarded", {
  soil <- data.frame(Eh=1:5, pH=seq(5,6,length.out=5))
  z <- suppressWarnings(rri_pipeline(
    soil = soil, domain_weights = c(Physio = .1, Soil = .8, Micro = .1),
    method_soil = "scale", direction_anchor_soil = "Eh"))
  expect_true(all(is.na(z$row_scores$Physio)))
  expect_true(all(is.na(z$row_scores$Micro)))
  expect_equal(z$row_scores$RRI, z$row_scores$Soil)
  ## domain_coverage and n_domains are per-row columns (one value per
  ## observation), so assert the constant value across all rows.
  expect_equal(unique(z$row_scores$domain_coverage), .8)
  expect_equal(unique(z$row_scores$n_domains), 1)
  expect_length(z$row_scores$domain_coverage, nrow(soil))
})
test_that("hidden simulator states cannot enter observation scoring", {
  x <- data.frame(Eh=1:5, alpha_accept=seq(.1,.9,length.out=5),
                  Capacity=seq(.2,.8,length.out=5), Truth=seq(.8,.2,length.out=5))
  ## Collect every warning, then assert the hidden-column one was raised.
  ## (An unanchored-axis warning is also emitted and must not mask it.)
  seen <- character()
  z <- withCallingHandlers(
    rri_pipeline(soil = x, method_soil = "scale",
                 direction_anchor_soil = "Eh"),
    warning = function(cond) {
      seen <<- c(seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("hidden", seen, fixed = TRUE)))
  expect_equal(z$row_scores$Soil, .rri_scale(x$Eh))
})
test_that("nonnumeric data and duplicate dynamic times fail clearly", {
  ## Unanchored phys/micro axes emit an intentional warning before these
  ## errors are reached; suppress it so each expectation asserts on the error.
  expect_error(
    suppressWarnings(rri_pipeline(soil = data.frame(Eh = factor(1:4)))),
    "numeric"
  )
  id <- data.frame(plot = "P1", time = c(1, 1, 2, 3))
  expect_error(
    suppressWarnings(rri_pipeline(
      soil = data.frame(Eh = 1:4), id = id, mode = "rolling",
      time_col = "time", group_cols = "plot", method_soil = "scale",
      direction_anchor_soil = "Eh"
    )),
    "Duplicate times"
  )
})
