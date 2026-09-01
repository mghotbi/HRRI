#' @title Plot a recovery landscape from RRI perturbation-recovery metrics
#'
#' @description
#' Visualises user-level perturbation-recovery metrics from
#' `rri_recovery_metrics()`. Each row is one trajectory and each column is a
#' core resilience metric. Values are scaled within metric columns to make
#' amplitudes, overshoot, incomplete recovery, and recovery times comparable.
#'
#' @param rec A data frame returned by `rri_recovery_metrics()`.
#' @param group_cols Character vector identifying trajectory labels.
#' @param metrics Character vector of recovery metric columns to plot.
#' @param order_by Character scalar. Metric used to order trajectories.
#' @param base_size Numeric. Base font size.
#'
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_point geom_text
#' @importFrom ggplot2 scale_fill_gradientn scale_color_manual labs theme_minimal
#' @importFrom ggplot2 theme element_text element_blank margin
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#' @importFrom stats ave
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2,
#'   n_depth = 3,
#'   n_plant = 2,
#'   n_time = 12,
#'   p_micro = 20,
#'   seed = 109
#' )
#'
#' res <- suppressWarnings(rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   reducer = "per_domain",
#'   scaling = "pnorm"
#' ))
#'
#' rec <- rri_recovery_metrics(
#'   res = res,
#'   id = sim$id,
#'   time_col = "time",
#'   group_cols = c("plot", "depth", "plant_id"),
#'   perturb_start = 5,
#'   perturb_end = 7
#' )
#'
#' head(rec)
#' plot_rri_recovery_landscape(rec)
#'
#' @export
plot_rri_recovery_landscape <- function(
    rec,
    group_cols = c("plot", "depth", "plant_id"),
    metrics = c("depth_min_frac", "overshoot_frac", "I_norm", "k",
                "tau_lag", "t_half"),
    order_by = "I_norm",
    base_size = 12
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_rri_recovery_landscape()` requires {ggplot2}.", call. = FALSE)
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("`plot_rri_recovery_landscape()` requires {tidyr}.", call. = FALSE)
  }

  if (!requireNamespace("tidyselect", quietly = TRUE)) {
    stop("`plot_rri_recovery_landscape()` requires {tidyselect}.", call. = FALSE)
  }

  rec <- as.data.frame(rec)

  missing_metrics <- setdiff(metrics, names(rec))
  if (length(missing_metrics) > 0) {
    stop(
      "Missing metric columns: ",
      paste(missing_metrics, collapse = ", "),
      call. = FALSE
    )
  }

  ## rri_recovery_metrics() assigns no trajectory class: doing so needs
  ## thresholds the package deliberately leaves to the user. Derive a
  ## descriptive label when none is supplied. These describe the score
  ## trajectory only and identify no mechanism.
  if (!"trajectory_class" %in% names(rec)) {
    flag <- if ("displaced_plateau_flag" %in% names(rec))
      rec$displaced_plateau_flag %in% TRUE else rep(FALSE, nrow(rec))
    inc <- if ("incomplete_return_frac" %in% names(rec))
      rec$incomplete_return_frac else rep(NA_real_, nrow(rec))
    rec$trajectory_class <- ifelse(
      flag, "displaced plateau",
      ifelse(is.finite(inc) & inc < -0.10, "incomplete return", "returned"))
    message("`trajectory_class` not supplied; derived from ",
            "displaced_plateau_flag and incomplete_return_frac.")
  }

  group_cols <- intersect(group_cols, names(rec))

  if (length(group_cols) == 0) {
    rec$.trajectory <- paste0("Trajectory ", seq_len(nrow(rec)))
  } else {
    rec$.trajectory <- apply(
      rec[, group_cols, drop = FALSE],
      1,
      paste,
      collapse = " | "
    )
  }

  if (!order_by %in% names(rec)) {
    order_by <- metrics[1]
  }

  rec$.order_value <- suppressWarnings(as.numeric(rec[[order_by]]))
  rec <- rec[order(rec$.order_value, decreasing = TRUE, na.last = TRUE), ]

  rec$.trajectory <- factor(rec$.trajectory, levels = rev(rec$.trajectory))

  plot_df <- rec[, c(".trajectory", "trajectory_class", metrics), drop = FALSE]

  for (metric in metrics) {
    plot_df[[metric]] <- suppressWarnings(as.numeric(plot_df[[metric]]))
  }

  long_df <- tidyr::pivot_longer(
    plot_df,
    cols = tidyselect::all_of(metrics),
    names_to = "metric",
    values_to = "value"
  )

  long_df$value_scaled <- stats::ave(
    long_df$value,
    long_df$metric,
    FUN = function(x) {
      if (all(is.na(x))) {
        return(rep(NA_real_, length(x)))
      }

      r <- range(x, na.rm = TRUE)

      if (!all(is.finite(r)) || diff(r) == 0) {
        return(rep(0.5, length(x)))
      }

      (x - r[1]) / diff(r)
    }
  )

  metric_labels <- c(
    depth_min_frac = "Resistance\nloss",
    overshoot_frac = "Overshoot",
    I_norm         = "Incomplete\nrecovery",
    k              = "Recovery\nrate",
    tau_lag        = "Response\nlag",
    t_half         = "Half-\nrecovery",
    A_norm = "Resistance\nloss", O_norm = "Overshoot", tau_r = "Recovery\ntime"
  )

  long_df$metric_label <- metric_labels[long_df$metric]
  long_df$metric_label[is.na(long_df$metric_label)] <-
    long_df$metric[is.na(long_df$metric_label)]

  class_cols <- c(
    fast_recovery = "#2E7D32",
    slow_recovery = "#8C6D31",
    overshoot = "#2166AC",
    hysteresis = "#7B3294",
    incomplete_recovery = "#B2182B",
    unclassified = "grey50"
  )

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x = .data$metric_label,
      y = .data$.trajectory,
      fill = .data$value_scaled
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.45,
      width = 0.96,
      height = 0.9
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(is.na(.data$value), "", signif(.data$value, 2))
      ),
      size = base_size / 4,
      colour = "#111111"
    ) +
    ggplot2::geom_point(
      data = unique(long_df[, c(".trajectory", "trajectory_class")]),
      ggplot2::aes(
        x = 0.42,
        y = .data$.trajectory,
        colour = .data$trajectory_class
      ),
      inherit.aes = FALSE,
      size = 3.2
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 1),
      na.value = "grey93",
      name = "Scaled\nmetric"
    ) +
    ggplot2::scale_color_manual(
      values = class_cols,
      na.value = "grey60",
      name = "Trajectory\nclass"
    ) +
    ggplot2::labs(
      title = "Redox resilience recovery landscape",
      subtitle = paste0(
        "Trajectories ordered by ",
        order_by,
        "; colour intensity shows within-metric scaled magnitude"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        colour = "#111111",
        size = base_size,
        margin = ggplot2::margin(t = 8)
      ),
      axis.text.y = ggplot2::element_text(
        colour = "#222222",
        size = base_size * 0.68
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 6,
        colour = "#111111"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size,
        colour = "#444444",
        margin = ggplot2::margin(b = 12)
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.margin = ggplot2::margin(15, 20, 15, 20)
    )
}
