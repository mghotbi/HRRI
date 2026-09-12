#' @title Four-panel diagnostic figure for an accuracy assessment
#'
#' @description
#' Draws the figure that accompanies [rri_accuracy()]. Each panel answers a
#' question a single correlation coefficient cannot: whether the score is
#' calibrated, whether disagreement grows with level, how much precision the
#' clustering costs, and which kind of error dominates.
#'
#' @param acc An object of class `rri_accuracy` from [rri_accuracy()], created
#'   with `n_boot > 0` so that the resampled statistics are available.
#' @param panels Character vector selecting panels, any of `"calibration"`,
#'   `"agreement"`, `"precision"` and `"error"`. Defaults to all four.
#' @param score_label,target_label Axis labels for the score and the reference
#'   target.
#' @param point_alpha Opacity of the individual observations. Lower it when
#'   trajectories overplot.
#' @param show_clusters Logical, or `NULL` to decide automatically. Colours
#'   observations by cluster. Set `FALSE` above roughly 20 clusters, where the
#'   colouring stops being informative.
#' @param base_size Base font size passed to [theme_ems()].
#' @param ncol Number of columns in the assembled figure. Ignored when
#'   **patchwork** is unavailable.
#'
#' @return If **patchwork** is installed, a single assembled `patchwork`
#'   object. Otherwise a named list of `ggplot` objects, so nothing is lost
#'   when the suggested package is absent.
#'
#' @details
#' **Panel A, calibration.** Score against target, with the 1:1 line dashed and
#' the fitted line solid. Perfect agreement puts the points on the dashed line;
#' a solid line flatter than it means the score compresses the target's range,
#' and one displaced from it means a systematic bias. Open points are cluster
#' means, the level at which these units are independent.
#'
#' **Panel B, agreement.** A Bland-Altman plot: the difference between score and
#' target against their mean, with the mean difference and the limits of
#' agreement. A scatter that fans out, or that slopes, shows that disagreement
#' depends on level, which a correlation coefficient cannot reveal. Because
#' observations are clustered, the limits come from cluster means; row-level
#' limits would be far too tight.
#'
#' **Panel C, precision.** The bootstrap sampling distribution of \eqn{r} under
#' row resampling and under trajectory resampling, with both intervals drawn
#' beneath. The difference in width is the cost of treating repeated
#' observations of one unit as independent observations of many. The
#' permutation null, when computed, sits behind them for reference.
#'
#' **Panel D, error.** Mean squared error split into squared bias, variance
#' mismatch and lack of correlation. The three sum to the mean squared error
#' exactly, so the panel is a partition rather than an approximation.
#'
#' Colours follow the package's chemistry-derived palette: teal for redox, rust
#' for iron, violet for manganese, ochre for cautionary annotation.
#'
#' @examples
#' set.seed(1)
#' k <- 8; m <- 20
#' unit   <- rnorm(k, 0, 0.30)
#' target <- unlist(lapply(unit, function(u) u + 0.5 + rnorm(m, 0, 0.05)))
#' score  <- 0.75 * target + 0.10 + rnorm(k * m, 0, 0.06)
#' traj   <- rep(seq_len(k), each = m)
#'
#' acc <- rri_accuracy(score, target, cluster = traj,
#'                     n_boot = 200, n_perm = 200, seed = 1)
#' p <- plot_rri_accuracy(acc)
#' \donttest{
#' print(p)
#' }
#'
#' @seealso [rri_accuracy()] for the statistics the figure displays.
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_hline geom_vline
#'   geom_col geom_text geom_density labs theme element_text element_blank
#'   element_rect scale_colour_manual scale_fill_manual scale_y_continuous
#'   coord_flip annotate after_stat
#' @importFrom stats sd quantile coef lm setNames
#' @export
plot_rri_accuracy <- function(acc,
                              panels = c("calibration", "agreement",
                                         "precision", "error"),
                              score_label = "Score",
                              target_label = "Reference target",
                              point_alpha = 0.45,
                              show_clusters = NULL,
                              base_size = 11,
                              ncol = 2) {

  if (!inherits(acc, "rri_accuracy")) {
    stop("`acc` must be an object returned by rri_accuracy().", call. = FALSE)
  }
  if (is.null(acc$data) || !nrow(acc$data)) {
    stop("`acc` carries no data; it may come from an older version of ",
         "rri_accuracy(). Re-run it to regenerate.", call. = FALSE)
  }
  panels <- match.arg(panels, several.ok = TRUE)

  ## ---- palette (mirrors vignettes/css/hrri.css) ---------------------------
  col_ink   <- "#1c2321"
  col_soft  <- "#4a5451"
  col_rule  <- "#c9d3d0"
  col_redox <- "#2f6b6b"   # deep teal
  col_fe    <- "#8c4a2f"   # rust
  col_mn    <- "#6b5b8a"   # violet
  col_warn  <- "#9a6a24"   # ochre

  d <- acc$data
  d$cluster <- as.factor(d$cluster)
  k <- nlevels(d$cluster)
  if (is.null(show_clusters)) show_clusters <- k <= 20L

  ## Cluster means: the level at which these units are independent.
  cm <- data.frame(
    cluster = levels(d$cluster),
    score   = as.numeric(tapply(d$score,  d$cluster, mean)),
    target  = as.numeric(tapply(d$target, d$cluster, mean)),
    stringsAsFactors = FALSE
  )

  ag <- stats::setNames(acc$agreement$row_level, acc$agreement$statistic)

  base_theme <- theme_ems(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = base_size,
                                              colour = col_ink),
      plot.subtitle   = ggplot2::element_text(size = base_size - 2.5,
                                              colour = col_soft),
      axis.title      = ggplot2::element_text(face = "plain", colour = col_soft),
      axis.text       = ggplot2::element_text(colour = col_soft),
      legend.position = "none"
    )

  out <- list()

  ## =========================================================================
  ## A. Calibration
  ## =========================================================================
  if ("calibration" %in% panels) {
    rng <- range(c(d$score, d$target), finite = TRUE)

    ## Drawn as score ~ target so the line lies in the same orientation as the
    ## points. acc$calibration reports target ~ score, the direction that
    ## should equal 1 for a calibrated score; both are labelled so the figure
    ## and the table cannot be read as contradicting each other.
    fit <- stats::coef(stats::lm(score ~ target, data = d))

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$target, y = .data$score))
    p <- if (show_clusters) {
      p + ggplot2::geom_point(ggplot2::aes(colour = .data$cluster),
                              size = 1.1, alpha = point_alpha)
    } else {
      p + ggplot2::geom_point(size = 1.1, alpha = point_alpha,
                              colour = col_redox)
    }

    p <- p +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                           colour = col_soft, linewidth = 0.5) +
      ggplot2::geom_abline(slope = fit[2], intercept = fit[1],
                           colour = col_fe, linewidth = 0.9) +
      ggplot2::geom_point(data = cm,
                          ggplot2::aes(x = .data$target, y = .data$score),
                          shape = 21, size = 2.8, stroke = 0.7,
                          fill = "white", colour = col_ink,
                          inherit.aes = FALSE) +
      ggplot2::annotate("text", x = rng[1], y = rng[2], hjust = 0, vjust = 1,
                        size = base_size / 3.4, colour = col_ink,
                        lineheight = 1.15,
                        label = sprintf("r = %.3f\nCCC = %.3f\nRMSE = %.3f",
                                        ag[["pearson_r"]], ag[["lins_ccc"]],
                                        ag[["rmse"]])) +
      ggplot2::labs(
        title = "A  Calibration",
        subtitle = sprintf(
          "dashed 1:1, fitted slope %.2f; open points are cluster means",
          fit[2]),
        x = target_label, y = score_label) +
      base_theme

    out$calibration <- p
  }

  ## =========================================================================
  ## B. Agreement (Bland-Altman)
  ## =========================================================================
  if ("agreement" %in% panels) {
    d$avg   <- (d$score + d$target) / 2
    d$diff  <- d$score - d$target
    cm$avg  <- (cm$score + cm$target) / 2
    cm$diff <- cm$score - cm$target

    ## Limits of agreement from cluster means, not rows. Row-level limits treat
    ## dependent observations as independent and come out far too tight.
    mu  <- mean(cm$diff)
    sg  <- if (nrow(cm) > 1L) stats::sd(cm$diff) else NA_real_
    loa <- c(mu - 1.96 * sg, mu + 1.96 * sg)

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$avg, y = .data$diff))
    p <- if (show_clusters) {
      p + ggplot2::geom_point(ggplot2::aes(colour = .data$cluster),
                              size = 1.1, alpha = point_alpha)
    } else {
      p + ggplot2::geom_point(size = 1.1, alpha = point_alpha,
                              colour = col_redox)
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, colour = col_rule, linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = mu, colour = col_fe, linewidth = 0.9)

    if (is.finite(sg)) {
      p <- p + ggplot2::geom_hline(yintercept = loa, colour = col_warn,
                                   linetype = "dotted", linewidth = 0.6)
    }

    p <- p +
      ggplot2::geom_point(data = cm,
                          ggplot2::aes(x = .data$avg, y = .data$diff),
                          shape = 21, size = 2.8, stroke = 0.7,
                          fill = "white", colour = col_ink,
                          inherit.aes = FALSE) +
      ggplot2::labs(
        title = "B  Agreement",
        subtitle = if (is.finite(sg)) {
          sprintf("bias %+.3f; 95%% limits [%.3f, %.3f] from cluster means",
                  mu, loa[1], loa[2])
        } else {
          "bias shown; too few clusters for limits of agreement"
        },
        x = "Mean of score and target",
        y = "Score minus target") +
      base_theme

    out$agreement <- p
  }

  ## =========================================================================
  ## C. Precision: what the clustering costs
  ## =========================================================================
  if ("precision" %in% panels) {
    dr <- acc$draws
    have <- !is.null(dr$naive_r) && !is.null(dr$cluster) &&
      any(is.finite(dr$naive_r)) && any(is.finite(dr$cluster[, "r"]))

    if (!have) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, size = base_size / 3,
                          colour = col_soft, lineheight = 1.2,
                          label = paste("Re-run rri_accuracy() with n_boot > 0",
                                        "to draw the sampling distributions",
                                        sep = "\n")) +
        ggplot2::labs(title = "C  Precision", x = NULL, y = NULL) +
        base_theme +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank())
    } else {
      cl_r <- dr$cluster[, "r"]
      bd <- rbind(
        data.frame(r = dr$naive_r[is.finite(dr$naive_r)],
                   src = "resampling rows (too narrow)",
                   stringsAsFactors = FALSE),
        data.frame(r = cl_r[is.finite(cl_r)],
                   src = "resampling trajectories",
                   stringsAsFactors = FALSE)
      )
      if (!is.null(dr$null_r) && any(is.finite(dr$null_r))) {
        bd <- rbind(bd, data.frame(r = dr$null_r[is.finite(dr$null_r)],
                                   src = "permutation null",
                                   stringsAsFactors = FALSE))
      }
      lev <- c("permutation null", "resampling rows (too narrow)",
               "resampling trajectories")
      bd$src <- factor(bd$src, levels = lev[lev %in% unique(bd$src)])

      pal <- c("permutation null"            = col_rule,
               "resampling rows (too narrow)" = col_mn,
               "resampling trajectories"     = col_redox)

      a    <- (1 - acc$conf) / 2
      ci_n <- unname(stats::quantile(dr$naive_r, c(a, 1 - a), na.rm = TRUE))
      ci_c <- unname(stats::quantile(cl_r,       c(a, 1 - a), na.rm = TRUE))

      p <- ggplot2::ggplot(bd, ggplot2::aes(x = .data$r, fill = .data$src,
                                            colour = .data$src)) +
        ggplot2::geom_density(
          ggplot2::aes(y = ggplot2::after_stat(.data$scaled)),
          alpha = 0.32, linewidth = 0.5) +
        ggplot2::scale_fill_manual(values = pal, drop = TRUE) +
        ggplot2::scale_colour_manual(values = pal, drop = TRUE) +
        ggplot2::geom_vline(xintercept = ag[["pearson_r"]],
                            colour = col_ink, linewidth = 0.6) +
        ggplot2::annotate("segment", x = ci_n[1], xend = ci_n[2],
                          y = -0.12, yend = -0.12, colour = col_mn,
                          linewidth = 1.5) +
        ggplot2::annotate("segment", x = ci_c[1], xend = ci_c[2],
                          y = -0.28, yend = -0.28, colour = col_redox,
                          linewidth = 1.5) +
        ggplot2::annotate("text", x = mean(ci_n), y = -0.04,
                          label = sprintf("rows: width %.3f", diff(ci_n)),
                          size = base_size / 3.8, colour = col_mn) +
        ggplot2::annotate("text", x = mean(ci_c), y = -0.38,
                          label = sprintf("trajectories: width %.3f",
                                          diff(ci_c)),
                          size = base_size / 3.8, colour = col_redox) +
        ggplot2::labs(
          title = "C  Precision",
          subtitle = sprintf(
            "ICC %.2f, design effect %.1f: effective n %.0f of %d rows",
            acc$dependence$icc, acc$dependence$design_effect,
            acc$dependence$effective_n, acc$dependence$n_observations),
          x = "Pearson r", y = "Bootstrap density (scaled)") +
        base_theme +
        ggplot2::theme(
          legend.position = "bottom",
          legend.title    = ggplot2::element_blank(),
          legend.key      = ggplot2::element_rect(fill = NA, colour = NA),
          legend.text     = ggplot2::element_text(size = base_size - 3.5))
    }
    out$precision <- p
  }

  ## =========================================================================
  ## D. Error decomposition
  ## =========================================================================
  if ("error" %in% panels) {
    dc  <- acc$decomposition
    lab <- c(squared_bias        = "Squared bias",
             variance_mismatch   = "Variance mismatch",
             lack_of_correlation = "Lack of correlation")
    dc$label <- factor(unname(lab[dc$component]), levels = rev(unname(lab)))
    pal2 <- stats::setNames(c(col_redox, col_warn, col_fe), levels(dc$label))

    ymax <- max(120, max(dc$percent, na.rm = TRUE) * 1.18)

    p <- ggplot2::ggplot(dc, ggplot2::aes(x = .data$label, y = .data$percent,
                                          fill = .data$label)) +
      ggplot2::geom_col(width = 0.62) +
      ggplot2::scale_fill_manual(values = pal2) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", .data$percent)),
        hjust = -0.18, size = base_size / 3.4, colour = col_ink) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
      ggplot2::labs(
        title = "D  Where the error is",
        subtitle = sprintf("MSE %.4g, partitioned exactly (residual %.1e)",
                           attr(dc, "mse"), attr(dc, "residual")),
        x = NULL, y = "Percent of mean squared error") +
      base_theme

    out$error <- p
  }

  out <- out[intersect(c("calibration", "agreement", "precision", "error"),
                       names(out))]

  if (length(out) == 1L) {
    out[[1]]
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    Reduce(`+`, out) + patchwork::plot_layout(ncol = ncol)
  } else {
    message("Install 'patchwork' to assemble these panels into one figure; ",
            "returning them as a list.")
    out
  }
}
