#' @title Radar Chart of Available HRRI Diagnostics
#'
#' @description
#' Displays available diagnostic summaries labelled Capacity, Connectivity,
#' Kinetics and Memory alongside the composite RRI score. These axes are
#' operational descriptors returned by \code{\link{rri_property_scores}};
#' they are not direct measurements or identified estimates of the theoretical
#' mechanisms bearing the same names.
#'
#' @param props A list returned by \code{\link{rri_property_scores}}, or a
#'   named numeric vector with elements \code{Capacity}, \code{Connectivity},
#'   \code{Kinetics}, \code{Memory} (values in `[0, 1]`).
#' @param rri_value Optional numeric. Composite RRI to display in the chart
#'   centre annotation. Defaults to \code{props$rri_summary} if available.
#' @param group_list Optional named list of property score vectors, one per
#'   group (e.g., per thaw stage or treatment). If supplied, multiple
#'   overlapping polygons are drawn, one per group.
#' @param fill_alpha Numeric in `[0, 1]`. Polygon fill transparency.
#' @param colours Character vector of polygon outline/fill colours, recycled
#'   across groups.
#' @param show_values Logical. Annotate each axis tip with the numeric score.
#' @param title Character. Plot title.
#' @param base_size Numeric. Base font size.
#'
#' @details
#' The chart uses Cartesian coordinates constructed with \code{ggplot2}; no
#' external radar-chart package is required. Each available axis runs from 0
#' (centre) to 1 (rim). Polygon area has no quantitative meaning, and axes based
#' on different transformations are not necessarily commensurable.
#'
#' \strong{Axis meanings:}
#' \itemize{
#'   \item \strong{Capacity} --- oxidative-oriented soil feature composite; not
#'     accessible capacity unless the caller calculates and supplies it explicitly.
#'   \item \strong{Connectivity} --- association or network-topology descriptor;
#'     not demonstrated electron transfer.
#'   \item \strong{Kinetics} --- cohort- or timescale-relative recovery-speed descriptor;
#'     not a mineral exchange rate.
#'   \item \strong{Memory} --- loop-area and persistent-displacement descriptor;
#'     not an identified causal memory state.
#'   \item \strong{RRI} --- composite score under the declared scaling and weights.
#' }
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom rlang .data
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 3, n_depth = 2, n_plant = 3, n_time = 14,
#'   p_micro = 30, seed = 99
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
#' props <- rri_property_scores(
#'   res       = res,
#'   rec       = rec,
#'   soil_df   = sim$Eh_stability,
#'   eac_col   = "EAC",
#'   edc_col   = "EDC",
#'   humic_col = "dissolved_organic_matter_redox"
#' )
#'
#' plot_rri_properties(props)
#'
#' @export
plot_rri_properties <- function(
  props,
  rri_value = NULL,
  group_list = NULL,
  fill_alpha = 0.20,
  colours = c("#1A3A5C", "#E07B39", "#2E7D32", "#7B3294", "#B2182B"),
  show_values = TRUE,
  title = "HRRI Diagnostic Profile",
  base_size = 13
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`plot_rri_properties()` requires {ggplot2}.", call. = FALSE)
  }

  # -- extract scores -------------------------------------------------------------
  .extract_scores <- function(x, rri_val = NULL) {
    if (is.list(x) && "property_scores" %in% names(x)) {
      sc <- as.numeric(x$property_scores)
      names(sc) <- names(x$property_scores)
      if (is.null(rri_val)) rri_val <- x$rri_summary
    } else {
      sc <- stats::setNames(as.numeric(x), names(x))
    }
    props_needed <- c("Capacity", "Connectivity", "Kinetics", "Memory")
    missing_p <- setdiff(props_needed, names(sc))
    if (length(missing_p) > 0) {
      sc[missing_p] <- NA_real_
    }
    if (!is.null(rri_val) && is.finite(rri_val)) {
      sc["RRI"] <- as.numeric(rri_val)
    } else if (!"RRI" %in% names(sc)) {
      sc["RRI"] <- NA_real_
    }
    sc
  }

  # -- axis order -----------------------------------------------------------------
  axis_order <- c("Capacity", "Connectivity", "Kinetics", "Memory", "RRI")
  n_axes <- length(axis_order)

  # -- build polygon data ---------------------------------------------------------
  # angles: n_axes spokes evenly spaced; close polygon by repeating first
  angles <- seq(0, 2 * pi, length.out = n_axes + 1)[seq_len(n_axes)]
  # The first spoke is placed at angle zero; labels describe diagnostics only.

  build_polygon <- function(scores, group_name) {
    vals <- pmax(0, pmin(1, as.numeric(scores[axis_order])))
    # Missing properties remain missing, not zero.
    # close polygon
    vals_c <- c(vals, vals[1])
    angles_c <- c(angles, angles[1])
    data.frame(
      group = group_name,
      axis = c(axis_order, axis_order[1]),
      angle = angles_c,
      value = vals_c,
      x = vals_c * cos(angles_c),
      y = vals_c * sin(angles_c),
      stringsAsFactors = FALSE
    )
  }

  # grid rings
  grid_rings <- seq(0.25, 1, by = 0.25)
  grid_dfs <- lapply(grid_rings, function(r) {
    a_seq <- seq(0, 2 * pi, length.out = 100)
    data.frame(
      x     = r * cos(a_seq),
      y     = r * sin(a_seq),
      ring  = r
    )
  })
  grid_df <- do.call(rbind, grid_dfs)

  # axis lines (spokes)
  spoke_df <- data.frame(
    x0    = 0,
    y0    = 0,
    x1    = cos(angles),
    y1    = sin(angles),
    axis  = axis_order
  )

  # axis labels (slightly beyond rim)
  label_r <- 1.18
  label_df <- data.frame(
    axis = axis_order,
    x = label_r * cos(angles),
    y = label_r * sin(angles),
    stringsAsFactors = FALSE
  )

  # -- polygon(s) ----------------------------------------------------------------
  if (!is.null(group_list) && is.list(group_list)) {
    poly_list <- lapply(seq_along(group_list), function(i) {
      sc <- .extract_scores(group_list[[i]])
      build_polygon(sc, names(group_list)[i])
    })
    poly_df <- do.call(rbind, poly_list)
    colours <- rep(colours, length.out = length(group_list))
  } else {
    sc <- .extract_scores(props, rri_value)
    poly_df <- build_polygon(sc, "system")
    colours <- colours[1]
    if (!is.null(rri_value)) {
      rri_val <- rri_value
    } else if (is.list(props) && !is.null(props$rri_summary)) {
      rri_val <- props$rri_summary
    } else {
      rri_val <- NA_real_
    }
  }

  groups <- unique(poly_df$group)
  colour_map <- stats::setNames(colours[seq_along(groups)], groups)

  # -- tip value labels ----------------------------------------------------------
  tip_df <- NULL
  if (isTRUE(show_values) && is.null(group_list)) {
    sc <- .extract_scores(props, rri_value)
    tip_df <- data.frame(
      axis = axis_order,
      x = 1.05 * cos(angles),
      y = 1.05 * sin(angles),
      label = sprintf("%.2f", pmax(0, pmin(1, as.numeric(sc[axis_order])))),
      stringsAsFactors = FALSE
    )
  }

  # -- assemble plot -------------------------------------------------------------
  p <- ggplot2::ggplot() +
    # grid rings
    ggplot2::geom_path(
      data = grid_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$ring),
      colour = "grey80", linewidth = 0.3
    ) +
    # grid ring labels at top (angle ~ pi/2)
    ggplot2::geom_text(
      data = data.frame(
        r     = grid_rings,
        x     = grid_rings * cos(pi / 2 + 0.08),
        y     = grid_rings * sin(pi / 2 + 0.08)
      ),
      ggplot2::aes(
        x = .data$x, y = .data$y,
        label = .data$r
      ),
      size = base_size / 5,
      colour = "grey55"
    ) +
    # spokes
    ggplot2::geom_segment(
      data = spoke_df,
      ggplot2::aes(
        x = .data$x0, y = .data$y0,
        xend = .data$x1, yend = .data$y1
      ),
      colour = "grey70",
      linewidth = 0.4
    )

  # polygons (one per group)
  for (i in seq_along(groups)) {
    g <- groups[i]
    col <- colour_map[[g]]
    p_sub <- poly_df[poly_df$group == g, ]
    if (all(is.finite(p_sub$value))) p <- p +
      ggplot2::geom_polygon(
        data = p_sub,
        ggplot2::aes(x = .data$x, y = .data$y),
        fill = col,
        colour = col,
        alpha = fill_alpha,
        linewidth = 1.1
      )
    p <- p + ggplot2::geom_point(
        data = p_sub[seq_len(nrow(p_sub) - 1), ], # exclude closing duplicate
        ggplot2::aes(x = .data$x, y = .data$y),
        colour = col, size = 3, shape = 21,
        fill = col, stroke = 1.2
      )
  }

  # axis labels
  p <- p +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$axis),
      fontface = "bold",
      size = base_size / 3.2,
      colour = "#1A3A5C"
    )

  # value annotations at axis tips
  if (!is.null(tip_df)) {
    p <- p +
      ggplot2::geom_label(
        data = tip_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = base_size / 4.5,
        fill = "white",
        colour = "#1A3A5C",
        label.size = 0.2,
        label.padding = grid::unit(0.12, "lines")
      )
  }

  # group legend if multiple groups
  if (length(groups) > 1) {
    leg_df <- data.frame(
      group = groups,
      colour = as.character(colour_map[groups]),
      stringsAsFactors = FALSE
    )
    p <- p +
      ggplot2::geom_point(
        data = leg_df,
        ggplot2::aes(
          x      = 0.88,
          y      = seq(from = -0.85, by = -0.13, length.out = nrow(leg_df)),
          colour = .data$group
        ),
        size = 4
      ) +
      ggplot2::geom_text(
        data = leg_df,
        ggplot2::aes(
          x     = 0.95,
          y     = seq(from = -0.85, by = -0.13, length.out = nrow(leg_df)),
          label = .data$group
        ),
        hjust = 0,
        size = base_size / 4.5,
        colour = "#333333"
      ) +
      ggplot2::scale_colour_manual(values = colour_map, guide = "none")
  }

  # RRI annotation in centre (single-group mode)
  if (is.null(group_list)) {
    centre_val <- if (!is.null(rri_value) && is.finite(rri_value)) {
      rri_value
    } else if (is.list(props) && !is.null(props$rri_summary) &&
      is.finite(props$rri_summary)) {
      props$rri_summary
    } else {
      NA_real_
    }

    if (!is.na(centre_val)) {
      p <- p +
        ggplot2::annotate(
          "label",
          x = 0, y = 0,
          label = sprintf("RRI\n%.3f", centre_val),
          size = base_size / 3.5,
          fontface = "bold",
          fill = "white",
          colour = "#1A3A5C",
          label.size = 0.3
        )
    }
  }

  p <- p +
    ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
    ggplot2::labs(
      title = title,
      subtitle = "Operational descriptors; polygon area is not a mechanistic quantity"
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 3, hjust = 0.5,
        colour = "#1A3A5C"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size * 0.80, hjust = 0.5, colour = "#555555"
      ),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )

  p
}
