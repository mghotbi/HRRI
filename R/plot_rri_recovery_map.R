#' @title Plot RRI Recovery Map
#'
#' @description
#' Visualises per-group RRI trajectories through baseline, perturbation and
#' recovery phases as a tile-and-line map. Each row is one trajectory group;
#' time proceeds along the x-axis; tile fill encodes RRI magnitude; vertical
#' bands mark the perturbation window; and trajectory class is annotated on
#' the right margin.
#'
#' landscape shows cross-metric comparison per trajectory, while the recovery
#' map shows temporal RRI dynamics per group.
#'
#' @param res An object returned by \code{\link{rri_pipeline_st}}.
#' @param id  A data frame of experimental identifiers (same rows as
#'   \code{res$row_scores}), containing at minimum \code{time_col} and the
#'   columns in \code{group_cols}.
#' @param rec Optional data frame from \code{\link{rri_recovery_metrics}}.
#'   If supplied, trajectory class annotations are added to the right margin.
#' @param time_col Character. Name of the time column in \code{id}.
#' @param group_cols Character vector. Columns in \code{id} defining
#'   trajectory groups (e.g., \code{c("plot", "depth", "plant_id")}).
#' @param perturb_start Numeric. Start of perturbation phase (same units as
#'   \code{time_col}).
#' @param perturb_end Numeric. End of perturbation phase.
#' @param palette Character. Viridis palette option for RRI fill.
#' @param base_size Numeric. Base font size.
#' @param max_groups Integer. Maximum number of trajectory groups to display.
#'   Groups are sampled if the total exceeds this value.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom rlang .data
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2, n_depth = 2, n_plant = 3, n_time = 14,
#'   p_micro = 20, seed = 101
#' )
#'
#' res <- rri_pipeline_st(
#'   ROS_flux     = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data   = sim$micro_data,
#'   id           = sim$id,
#'   reducer      = "per_domain",
#'   scaling      = "pnorm"
#' )
#'
#' rec <- rri_recovery_metrics(
#'   res           = res,
#'   id            = sim$id,
#'   time_col      = "time",
#'   group_cols    = c("plot", "depth", "plant_id"),
#'   perturb_start = 5,
#'   perturb_end   = 8
#' )
#'
#' plot_rri_recovery_map(
#'   res           = res,
#'   id            = sim$id,
#'   rec           = rec,
#'   time_col      = "time",
#'   group_cols    = c("plot", "depth", "plant_id"),
#'   perturb_start = 5,
#'   perturb_end   = 8
#' )
#'
#' @export
plot_rri_recovery_map <- function(
  res,
  id,
  rec = NULL,
  time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = NULL,
  perturb_end = NULL,
  palette = "plasma",
  base_size = 11,
  max_groups = 40L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_rri_recovery_map()` requires {ggplot2}.", call. = FALSE)
  }

  # -- assemble long data frame --------------------------------------------------
  rs <- as.data.frame(res$row_scores)
  id_df <- as.data.frame(id)

  if (!"RRI" %in% names(rs)) {
    stop("`res$row_scores` must contain an RRI column.", call. = FALSE)
  }
  if (!time_col %in% names(id_df)) {
    stop(sprintf("time_col '%s' not found in `id`.", time_col), call. = FALSE)
  }

  if (nrow(id_df) != nrow(rs)) stop("id and row_scores must have equal row counts.")
  if (any(!group_cols %in% names(id_df))) stop("group_cols absent from id.")
  valid_groups <- group_cols

  df <- cbind(id_df, RRI = as.numeric(rs$RRI))

  if (length(valid_groups) == 0) {
    df$.group <- "all"
  } else {
    df$.group <- apply(
      df[, valid_groups, drop = FALSE], 1,
      paste,
      collapse = " | "
    )
  }

  df$.time <- as.numeric(df[[time_col]])

  # -- limit displayed groups ----------------------------------------------------
  unique_groups <- unique(df$.group)
  if (length(unique_groups) > max_groups) {
    n_total_groups <- length(unique_groups)
    unique_groups <- utils::head(sort(unique_groups), max_groups)
    df <- df[df$.group %in% unique_groups, , drop = FALSE]
    message(sprintf(
      "plot_rri_recovery_map: displaying %d of %d groups (set max_groups to override).",
      max_groups, n_total_groups
    ))
  }

  # -- order groups by mean RRI (descending) -------------------------------------
  group_means <- tapply(df$RRI, df$.group, mean, na.rm = TRUE)
  ordered_groups <- names(sort(group_means, decreasing = TRUE, na.last = TRUE))
  df$.group <- factor(df$.group, levels = rev(ordered_groups))

  # -- trajectory class annotation (right margin) --------------------------------
  class_cols <- c(
    fast_recovery       = "#2E7D32",
    slow_recovery       = "#8C6D31",
    overshoot           = "#2166AC",
    hysteresis          = "#7B3294",
    incomplete_recovery = "#B2182B",
    unclassified        = "grey50"
  )

  annot_df <- NULL
  if (!is.null(rec)) {
    rec_df <- as.data.frame(rec)
    if ("trajectory_class" %in% names(rec_df) && length(valid_groups) > 0) {
      rec_df$.group <- apply(
        rec_df[, intersect(valid_groups, names(rec_df)), drop = FALSE], 1,
        paste,
        collapse = " | "
      )
      # one row per group --- most severe class if multiple entries
      severity <- c(
        incomplete_recovery = 5, hysteresis = 4, overshoot = 3,
        slow_recovery = 2, fast_recovery = 1, unclassified = 0
      )
      rec_df$.sev <- severity[as.character(rec_df$trajectory_class)]
      rec_df$.sev[is.na(rec_df$.sev)] <- 0L

      annot_df <- do.call(rbind, lapply(split(rec_df, rec_df$.group), function(g) {
        g[which.max(g$.sev), c(".group", "trajectory_class"), drop = FALSE]
      }))
      annot_df$.group <- factor(annot_df$.group, levels = levels(df$.group))

      # x position for annotation = max time + small offset
      x_annot <- max(df$.time, na.rm = TRUE) * 1.02
      annot_df$.x <- x_annot
    }
  }

  # -- base tile map -------------------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$.time, y = .data$.group, fill = .data$RRI)
  ) +
    ggplot2::geom_tile(colour = NA, width = 0.9, height = 0.85) +
    ggplot2::scale_fill_viridis_c(
      option = palette,
      direction = -1,
      name = "RRI",
      limits = c(0, 1),
      na.value = "grey85"
    )

  # -- perturbation window bands -------------------------------------------------
  if (!is.null(perturb_start) && !is.null(perturb_end)) {
    band_df <- data.frame(
      xmin = as.numeric(perturb_start) - 0.5,
      xmax = as.numeric(perturb_end) + 0.5
    )
    p <- p +
      ggplot2::geom_rect(
        data = band_df,
        ggplot2::aes(
          xmin = .data$xmin, xmax = .data$xmax,
          ymin = -Inf, ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#D62728", alpha = 0.08
      ) +
      ggplot2::geom_vline(
        xintercept = c(
          as.numeric(perturb_start) - 0.5,
          as.numeric(perturb_end) + 0.5
        ),
        colour = "#D62728", linewidth = 0.5, linetype = "dashed"
      )
  }

  # -- trajectory class annotation dots -----------------------------------------
  if (!is.null(annot_df)) {
    p <- p +
      ggplot2::geom_point(
        data = annot_df,
        ggplot2::aes(
          x = .data$.x, y = .data$.group,
          colour = .data$trajectory_class
        ),
        inherit.aes = FALSE,
        size = 3, shape = 16
      ) +
      ggplot2::scale_colour_manual(
        values = class_cols,
        na.value = "grey60",
        name = "Trajectory\nclass"
      )
  }

  p <- p +
    ggplot2::labs(
      title = "RRI Recovery Map",
      subtitle = if (!is.null(perturb_start)) {
        sprintf(
          "Red band = perturbation window [%s, %s]",
          perturb_start, perturb_end
        )
      } else {
        "Tile fill = per-sample RRI; groups ordered by mean RRI"
      },
      x = time_col,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        size = base_size * 0.65,
        colour = "#333333"
      ),
      axis.text.x = ggplot2::element_text(size = base_size * 0.8),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 4
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size * 0.85,
        colour = "#555555"
      ),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(12, 24, 12, 12)
    )

  p
}
