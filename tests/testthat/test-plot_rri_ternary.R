library(testthat)

## skip_if_not_installed() calls requireNamespace(), and loading ggtern is
## itself what breaks ggplot2 >= 4.0.0 for the rest of the session -- including
## every later test file. The ggplot2 version is therefore checked first, and
## ggtern is never touched when it cannot work.
test_that("plot_RRI_ternary returns a ggtern plot", {
  skip_if_not_installed("ggplot2")
  skip_if(utils::packageVersion("ggplot2") >= "4.0.0",
          "ggtern is incompatible with ggplot2 >= 4.0.0; loading it would break later tests")
  skip_if_not_installed("ggtern")
  skip_if_not_installed("viridis")

  ternary_df <- data.frame(
    Physio = c(0.3, 0.4, 0.2),
    Soil   = c(0.4, 0.3, 0.5),
    Micro  = c(0.3, 0.3, 0.3),
    RRI    = c(0.6, 0.7, 0.5)
  )

  attr(ternary_df, "RRI_index") <- mean(ternary_df$RRI)

  p <- suppressWarnings(plot_RRI_ternary(ternary_df))

  expect_s3_class(p, "ggplot")
})

test_that("plot_RRI_ternary rejects non-compositional input", {
  skip_if_not_installed("ggplot2")
  skip_if(utils::packageVersion("ggplot2") >= "4.0.0",
          "ggtern is incompatible with ggplot2 >= 4.0.0; loading it would break later tests")
  skip_if_not_installed("ggtern")

  bad_df <- data.frame(
    Physio = c(0.5, 0.5),
    Soil   = c(0.5, 0.5),
    Micro  = c(0.5, 0.5),
    RRI    = c(0.6, 0.7)
  )

  expect_error(
    plot_RRI_ternary(bad_df),
    "must sum to 1"
  )
})

test_that("plot_RRI_ternary errors when required columns are missing", {
  skip_if_not_installed("ggplot2")
  skip_if(utils::packageVersion("ggplot2") >= "4.0.0",
          "ggtern is incompatible with ggplot2 >= 4.0.0; loading it would break later tests")
  skip_if_not_installed("ggtern")

  bad_df <- data.frame(
    Physio = 0.4,
    Soil   = 0.3
  )

  expect_error(
    plot_RRI_ternary(bad_df),
    "Missing required columns"
  )
})

test_that("plot_RRI_ternary drops invalid rows but errors if none remain", {
  skip_if_not_installed("ggplot2")
  skip_if(utils::packageVersion("ggplot2") >= "4.0.0",
          "ggtern is incompatible with ggplot2 >= 4.0.0; loading it would break later tests")
  skip_if_not_installed("ggtern")

  df <- data.frame(
    Physio = c(0, NA),
    Soil   = c(0, NA),
    Micro  = c(0, NA),
    RRI    = c(0.5, 0.6)
  )

  expect_error(
    plot_RRI_ternary(df),
    "No valid compositional rows"
  )
})
