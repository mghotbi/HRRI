test_that("property profiles retain missingness and never silently clip", {
  sc <- c(Capacity = 0.2, Connectivity = NA_real_, Kinetics = 0.8, Memory = 0.1)
  p <- plot_rri_properties(sc)
  expect_s3_class(p, "ggplot")
  expect_equal(sum(is.na(p$data$value)), 1)
  expect_true(grepl("no common", p$labels$caption))
  b <- ggplot2::ggplot_build(p)
  expect_true("Missing" %in% b$data[[3]]$label)
  expect_error(plot_rri_properties(replace(sc, 1, 1.1)), "not clipped")
  expect_error(plot_rri_properties(replace(sc, 1, Inf)), "not clipped")
  expect_error(plot_rri_properties(sc, group_list = list(sc)), "unique")
  expect_no_warning(ggplot2::ggplot_build(plot_rri_properties(sc, type = "radar")))
  expect_no_warning(ggplot2::ggplot_build(plot_rri_properties(sc,
    group_list = list(A = sc, B = sc))))
})

test_that("landscape retains unavailable and constant-column missing cells", {
  rec <- data.frame(k = c(0.1, 0.1, NA_real_), t_half = NA_real_,
    I_norm = c(.1, .2, .3), incomplete_return_frac = c(0.1, -.2, NA_real_))
  p <- suppressMessages(plot_rri_recovery_landscape(rec,
    metrics = c("k", "t_half", "I_norm")))
  d <- p$data
  expect_equal(sum(is.na(d$value_scaled)), 4)
  expect_equal(d$value_scaled[d$metric == "k" & is.finite(d$value)], c(.5, .5))
  expect_true(all(c("not flagged", "incomplete return", "unclassified") %in% d$trajectory_class))
  expect_true(any(grepl("n = 0/3", as.character(d$metric_label))))
  expect_no_warning(ggplot2::ggplot_build(p))
  expect_error(suppressMessages(plot_rri_recovery_landscape(rec, metrics = character())), "nonempty")
  p2 <- suppressMessages(plot_rri_recovery_landscape(rec, metrics = c("k", "t_half"), drop_empty = TRUE))
  expect_match(p2$labels$caption, "Omitted: t_half")
})

test_that("map aligns keys and preserves exact event times and missing observations", {
  id <- data.frame(plot = rep(c("P1", "P2"), each = 3), time = rep(c(-2, 0, 5), 2))
  rs <- transform(id, RRI = c(.1, NA, .3, .4, .5, .6))
  p <- plot_rri_recovery_map(list(row_scores = rs[6:1, ]), id,
    group_cols = "plot", perturb_start = -1.5, perturb_end = 2.2)
  expect_equal(p$data$RRI, rs$RRI)
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$data[[2]]$xintercept, c(-1.5, 2.2))
  expect_equal(b$data[[1]]$xmax - b$data[[1]]$xmin, rep(1.8, nrow(id)))
  expect_equal(nrow(b$data[[3]]), 1)
  expect_error(plot_rri_recovery_map(list(row_scores = rs), id,
    group_cols = "plot", perturb_start = 0), "both")
  rs$RRI[1] <- 2
  expect_error(plot_rri_recovery_map(list(row_scores = rs), id, group_cols = "plot"), "in ")
  expect_error(plot_rri_recovery_map(list(row_scores = rs), id, group_cols = "plot", max_groups = 0), "positive")
  dup <- id[c(1, 1, 3:6), ]
  expect_error(plot_rri_recovery_map(list(row_scores = data.frame(RRI = rep(.5, 6))),
    dup, group_cols = "plot"), "Duplicate")
})

test_that("agreement plots distinguish mean limits and preserve outliers", {
  target <- seq(0, 1, length.out = 40)
  score <- target * .6 + .2
  score[40] <- 3
  acc <- rri_accuracy(score, target, cluster = rep(1:4, each = 10), n_boot = 0, n_perm = 0)
  p <- plot_rri_accuracy(acc, panels = "calibration", style="diagnostic")
  b <- ggplot2::ggplot_build(p)
  expect_equal(max(b$data[[1]]$y), 3)
  expect_equal(p$coordinates$limits$x, range(c(score, target)))
  p <- plot_rri_accuracy(acc, panels = "agreement", style="diagnostic")
  expect_match(p$labels$subtitle, "Descriptive cluster-mean limits")
  cm <- tapply(score - target, rep(1:4, each = 10), mean)
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$data[[4]]$yintercept, mean(cm) + c(-1, 1) * 1.96 * sd(cm))
  perfect <- rri_accuracy(target, target, cluster = rep(1:4, each = 10), n_boot = 0, n_perm = 0)
  expect_no_warning(ggplot2::ggplot_build(plot_rri_accuracy(perfect, panels = "error")))
})

test_that("aliased landscape columns and sparse bootstrap draws remain distinct", {
  rec <- data.frame(k = c(.1, .2), k_recovery = c(.1, .2))
  p <- suppressMessages(plot_rri_recovery_landscape(rec, metrics = names(rec)))
  expect_equal(nlevels(p$data$metric_label), 2)
  y <- seq(.1, .9, length.out = 20)
  a <- rri_accuracy(y + sin(1:20) * .02, y,
    cluster = rep(1:4, each = 5), n_boot = 1, n_perm = 0, seed = 1)
  expect_no_warning(ggplot2::ggplot_build(plot_rri_accuracy(a, panels = "precision")))
})
