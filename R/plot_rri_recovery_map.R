#' @title Plot RRI Recovery Map
#'
#' @description
#' Visualises per-group RRI trajectories through baseline, perturbation and
#' recovery phases as an observation tile map. Each row is one trajectory group;
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
#' @param order_by Order trajectories by mean RRI (default) or sorted keys.
#' @param tile_gap Fraction of each sampling cell left as a gap, in [0, 1).
#' @param event_colour Colour of event boundary lines.
#' @param direction Viridis colour direction, either 1 or -1.
#' @param max_groups Integer. Maximum number of trajectory groups to display.
#'   The first groups in sorted label order are displayed when the total exceeds
#'   this value; the displayed fraction is reported. No random sampling occurs.
#' @details Tile width is 90 percent of the smallest distinct observed time
#' spacing. This prevents tiles implying continuous observation across long gaps.
#' Explicit missing scores are grey with a cross; unsampled times remain blank.
#' Event boundaries are shown at the exact supplied times. Scores must be in
#' `[0, 1]` or missing. Identifier keys are aligned when supplied in row_scores;
#' scores-only output must retain input row order. Duplicate group-time keys
#' and ambiguous recovery annotations are rejected.
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
  max_groups = 40L,
  order_by = c("mean", "key"), direction = -1, tile_gap = .1, event_colour = "#D62728"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_rri_recovery_map()` requires {ggplot2}.", call. = FALSE)
  }

  order_by <- match.arg(order_by)
  if (length(direction) != 1L || !direction %in% c(-1,1)) stop("direction must be 1 or -1.")
  if(length(tile_gap)!=1L || !is.finite(tile_gap) || tile_gap<0 || tile_gap>=1) stop("tile_gap must be in [0, 1).")
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

  if (!nrow(id_df)) stop("id has no rows.", call. = FALSE)
  if (length(max_groups) != 1L || !is.finite(max_groups) || max_groups < 1 || max_groups != floor(max_groups))
    stop("max_groups must be a positive integer.", call. = FALSE)
  if (!is.numeric(id_df[[time_col]]) || any(!is.finite(id_df[[time_col]])))
    stop("Time must be finite numeric.", call. = FALSE)
  if (anyNA(id_df[group_cols])) stop("Group identifiers must not be missing.", call. = FALSE)
  if (xor(is.null(perturb_start), is.null(perturb_end)))
    stop("Supply both disturbance interval endpoints.", call. = FALSE)
  if (!is.null(perturb_start)) {
    interval <- c(perturb_start, perturb_end)
    if (!is.numeric(interval) || length(interval) != 2L || any(!is.finite(interval)) || interval[1] >= interval[2])
      stop("Disturbance endpoints must be finite and increasing.", call. = FALSE)
  }
  df <- .rri_align_scores(rs, id_df, c(group_cols, time_col))
  if (anyDuplicated(.rri_key(id_df, c(group_cols, time_col))))
    stop("Duplicate group-time observation keys.", call. = FALSE)
  if (!is.numeric(df$RRI) || any(is.infinite(df$RRI)) || any(df$RRI < 0 | df$RRI > 1, na.rm = TRUE))
    stop("RRI must be numeric in [0, 1] or missing.", call. = FALSE)

  if (length(valid_groups) == 0) {
    df$.group <- "all"
  } else {
    df$.group <- do.call(paste, c(lapply(df[valid_groups], as.character), sep = " | "))
    if (length(unique(df$.group)) != length(unique(.rri_key(df, valid_groups))))
      stop("Group labels are ambiguous; avoid ' | ' within identifiers.", call. = FALSE)
  }

  df$.time <- as.numeric(df[[time_col]])

  # -- limit displayed groups ----------------------------------------------------
  times <- sort(unique(df$.time))
  spacing <- if (length(times) > 1L) min(diff(times)) else 1
  tile_width <- (1-tile_gap) * spacing
  unique_groups <- unique(df$.group)
  n_total_groups <- length(unique_groups)
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
  if (order_by == "key") ordered_groups <- sort(unique(as.character(df$.group)))
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
    if ("trajectory_class" %in% names(rec_df)) {
      if (!all(valid_groups %in% names(rec_df)))
        stop("Recovery annotations require all group_cols.", call. = FALSE)
      rec_df$.group <- if (length(valid_groups))
        do.call(paste, c(lapply(rec_df[valid_groups], as.character), sep = " | ")) else "all"
      if (anyNA(rec_df[valid_groups]) || anyDuplicated(rec_df$.group))
        stop("Recovery annotations require one unambiguous row per group.", call. = FALSE)
      annot_df <- rec_df[rec_df$.group %in% levels(df$.group), c(".group", "trajectory_class"), drop = FALSE]
      annot_df$.group <- factor(annot_df$.group, levels = levels(df$.group))
      annot_df$.x <- max(df$.time) + max(spacing, diff(range(df$.time)) * 0.03)
      extra <- setdiff(unique(as.character(annot_df$trajectory_class)), names(class_cols))
      extra <- extra[!is.na(extra)]
      class_cols <- c(class_cols, stats::setNames(rep("grey50", length(extra)), extra))
    }
  }

  # -- base tile map -------------------------------------------------------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$.time, y = .data$.group, fill = .data$RRI)
  ) +
    ggplot2::geom_tile(colour = NA, width = tile_width, height = 1-tile_gap) +
    ggplot2::scale_fill_viridis_c(
      option = palette,
      direction = direction,
      name = "RRI",
      limits = c(0, 1),
      na.value = "grey85"
    )

  # -- perturbation window bands -------------------------------------------------
  if (!is.null(perturb_start) && !is.null(perturb_end)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = c(
          as.numeric(perturb_start),
          as.numeric(perturb_end)
        ),
        colour = event_colour, linewidth = 0.5, linetype = "dashed"
      )
  }

  p <- p + ggplot2::geom_point(data = df[is.na(df$RRI), , drop = FALSE],
    ggplot2::aes(x = .data$.time, y = .data$.group), inherit.aes = FALSE,
    shape = 4, size = 1.8, colour = "grey35")

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
          "Dashed lines = supplied event boundaries [%s, %s]",
          perturb_start, perturb_end
        )
      } else {
        "Tile fill = per-sample RRI; groups ordered by mean RRI"
      },
      caption = sprintf("Showing %d/%d groups; grey cross = missing score; blank = unsampled time.\nTile width follows observed spacing; order = %s.",
                        length(unique_groups), n_total_groups, order_by),
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
      plot.caption = ggplot2::element_text(size = base_size * 0.65, hjust = 0),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(12, 24, 12, 12)
    )

  p
}
