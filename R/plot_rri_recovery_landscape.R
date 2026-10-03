#' @title Plot a recovery landscape from RRI perturbation-recovery metrics
#'
#' @description
#' Visualises trajectory-level recovery metrics from [rri_recovery_metrics()].
#' Each row is one trajectory and each column one recovery signature. Cell
#' colour encodes the within-column scaled magnitude; the printed number is
#' always the unscaled value, so nothing is hidden by the scaling.
#'
#' @param rec A data frame returned by [rri_recovery_metrics()].
#' @param group_cols Character vector of columns identifying one trajectory.
#' @param metrics Character vector of recovery metric columns to plot. Defaults
#'   to the columns returned by [rri_recovery_metrics()]. Legacy names
#'   (`A_norm`, `O_norm`, `tau_r`) are still labelled if supplied.
#' @param order_by Character scalar. Metric used to order trajectories.
#' @param orient Controls what darker colour means. `"concern"`
#'   inverts metrics for which a *smaller* value is the more concerning
#'   outcome. This is a descriptive polarity convention, not a shared scale of
#'   ecological concern; overshoot remains neutral. `"raw"` (default) scales every column upward, meaning dark is high-valued
#'   regardless of interpretation. See Details.
#' @param drop_empty Logical. Drop metric columns that are `NA` for every
#'   trajectory rather than drawing a blank column. `k` and `t_half` are `NA`
#'   unless a recovery rate was fitted, so an all-`NA` column is common and
#'   means "not estimable here", not "zero".
#' @param base_size Numeric. Base font size.
#'
#' @return A `ggplot` object.
#'
#' @details
#' **Why orientation matters.** The recovery signatures do not share a polarity.
#' A large `depth_min_frac` (deep decline), a large `tau_lag` (slow onset) and a
#' large `t_half` (slow return) are all unfavourable, but a large `k` is a
#' *fast* recovery rate and therefore favourable. Scaling every column upward
#' and applying one colour ramp would make dark mean "bad" in some columns and
#' "good" in others. With `orient = "concern"` the `k` column is inverted before
#' scaling so the ramp is interpretable across the panel. `overshoot_frac` is
#' treated as neutral and never inverted, because overshoot is not
#' unambiguously favourable or unfavourable.
#'
#' Scaling is min-max **within each column, within this cohort**. A dark cell
#' means "high relative to the other trajectories in this run", not high in any
#' absolute sense. Two datasets cannot be compared cell by cell.
#'
#' If `rec` has no `trajectory_class` column, one is derived from
#' `displaced_plateau_flag` and `incomplete_return_frac`. The derived labels
#' describe the score trajectory only and identify no mechanism. A negative
#' final displacement below -0.10 is labelled incomplete return; otherwise a
#' finite displacement is labelled not flagged, not evidence of equivalence.
#' All-missing columns are retained by default. Counts report finite values;
#' grey cells remain missing even in constant-valued columns. I_norm is the
#' capped absolute final displacement and does not encode its direction.
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
#' plot_rri_recovery_landscape(
#'   rec,
#'   metrics = c("depth_min_frac", "overshoot_frac", "I_norm",
#'               "k", "tau_lag", "t_half")
#' )
#'
#' @seealso [rri_recovery_metrics()], [plot_rri_recovery_map()]
#' @export
plot_rri_recovery_landscape <- function(
    rec,
    group_cols = c("plot", "depth", "plant_id"),
    metrics = c("depth_min_frac", "overshoot_frac", "I_norm", "k",
                "tau_lag", "t_half"),
    order_by = "I_norm",
    orient = c("raw", "concern"),
    drop_empty = FALSE,
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
  orient <- match.arg(orient)

  rec <- as.data.frame(rec)
  if (!nrow(rec)) stop("`rec` has no rows.", call. = FALSE)

  if (!length(metrics) || anyDuplicated(metrics))
    stop("metrics must be nonempty and unique.", call. = FALSE)
  dropped <- character()
  missing_metrics <- setdiff(metrics, names(rec))
  if (length(missing_metrics) > 0) {
    stop("Missing metric columns: ", paste(missing_metrics, collapse = ", "),
         ". Available: ", paste(names(rec), collapse = ", "), call. = FALSE)
  }

  ## ---- trajectory class ---------------------------------------------------
  ## rri_recovery_metrics() assigns no class: doing so needs thresholds the
  ## package deliberately leaves to the user. Derive a descriptive label when
  ## none is supplied. These describe the score trajectory, not a mechanism.
  derived_class <- FALSE
  if (!"trajectory_class" %in% names(rec)) {
    flag <- if ("displaced_plateau_flag" %in% names(rec)) {
      rec$displaced_plateau_flag %in% TRUE
    } else rep(FALSE, nrow(rec))
    inc <- if ("incomplete_return_frac" %in% names(rec)) {
      suppressWarnings(as.numeric(rec$incomplete_return_frac))
    } else rep(NA_real_, nrow(rec))
    rec$trajectory_class <- ifelse(
      flag, "displaced plateau",
      ifelse(is.finite(inc) & inc < -0.10, "incomplete return",
             ifelse(is.finite(inc), "not flagged", "unclassified"))
    )
    derived_class <- TRUE
    message("`trajectory_class` not supplied; derived from ",
            "displaced_plateau_flag and incomplete_return_frac.")
  }

  ## ---- trajectory labels --------------------------------------------------
  group_cols <- intersect(group_cols, names(rec))
  if (length(group_cols) == 0) {
    rec$.trajectory <- paste0("Trajectory ", seq_len(nrow(rec)))
  } else {
    ## do.call(paste) rather than apply(): apply() coerces the whole frame to a
    ## character matrix and reformats numerics unpredictably.
    rec$.trajectory <- do.call(
      paste, c(lapply(rec[group_cols], as.character), sep = " | ")
    )
    ## Duplicate labels would make factor(levels = ...) error. Disambiguate
    ## rather than fail: incomplete group_cols is a user mistake worth flagging,
    ## not a crash.
    if (anyDuplicated(rec$.trajectory)) {
      warning("`group_cols` do not uniquely identify trajectories; ",
              "appending row numbers. Supply the full key to avoid this.",
              call. = FALSE)
      rec$.trajectory <- paste0(rec$.trajectory, " #", seq_len(nrow(rec)))
    }
  }

  if (!order_by %in% names(rec)) order_by <- metrics[1]
  rec$.order_value <- suppressWarnings(as.numeric(rec[[order_by]]))
  rec <- rec[order(rec$.order_value, decreasing = TRUE, na.last = TRUE), ]
  rec$.trajectory <- factor(rec$.trajectory, levels = rev(rec$.trajectory))

  plot_df <- rec[, c(".trajectory", "trajectory_class", metrics), drop = FALSE]
  for (metric in metrics) {
    if (!is.numeric(plot_df[[metric]])) stop("Metric columns must be numeric.", call. = FALSE)
    plot_df[[metric]][!is.finite(plot_df[[metric]])] <- NA_real_
  }

  ## A metric that is NA for every trajectory renders as a blank grey stripe,
  ## which reads as a rendering fault rather than as missingness. k and t_half
  ## are the usual cases: both are NA unless fit_status is a fitted status.
  ## Report the omission instead of drawing an empty column.
  if (isTRUE(drop_empty)) {
    all_na <- vapply(plot_df[metrics], function(x) all(is.na(x)), logical(1))
    if (any(all_na)) {
      dropped <- metrics[all_na]
      message("Dropping metric(s) with no finite values: ",
              paste(dropped, collapse = ", "),
              ". Check `fit_status` - k and t_half are NA unless a rate was ",
              "fitted. Set drop_empty = FALSE to keep the empty columns.")
      metrics <- metrics[!all_na]
      if (!length(metrics)) stop("No metric has any finite value.", call. = FALSE)
      plot_df <- plot_df[, c(".trajectory", "trajectory_class", metrics),
                         drop = FALSE]
    }
  }

  long_df <- tidyr::pivot_longer(
    plot_df,
    cols      = tidyselect::all_of(metrics),
    names_to  = "metric",
    values_to = "value"
  )

  ## ---- direction-aware scaling -------------------------------------------
  ## Metrics for which a SMALLER value is the more concerning outcome. Only k
  ## qualifies among the defaults: a small recovery rate is a slow recovery.
  ## overshoot is deliberately absent - it is not unambiguously good or bad.
  lower_is_concerning <- c("k", "k_recovery")

  long_df$value_scaled <- stats::ave(
    long_df$value,
    long_df$metric,
    FUN = function(x) {
      out <- rep(NA_real_, length(x))
      ok <- is.finite(x)
      if (!any(ok)) return(out)
      r <- range(x[ok])
      out[ok] <- if (diff(r) == 0) 0.5 else (x[ok] - r[1]) / diff(r)
      out
    }
  )
  if (orient == "concern") {
    flip <- long_df$metric %in% lower_is_concerning
    long_df$value_scaled[flip] <- 1 - long_df$value_scaled[flip]
  }

  ## Labels are kept short deliberately: with six columns in a 7-inch figure
  ## each header has roughly one inch, and longer strings collide.
  metric_labels <- c(
    depth_min_frac = "Resistance\nloss",
    overshoot_frac = "Overshoot",
    I_norm         = "Absolute final\ndisplacement",
    k              = "Recovery\nrate",
    k_recovery     = "Recovery\nrate",
    tau_lag        = "Response\nlag",
    t_half         = "Half-life",
    H_norm         = "Hysteresis",
    ## legacy aliases so older metric tables still render
    A_norm = "Resistance\nloss", O_norm = "Overshoot", tau_r = "Recovery\ntime"
  )
  labels <- unname(metric_labels[metrics])
  labels[is.na(labels)] <- metrics[is.na(labels)]
  counts <- vapply(metrics, function(m) sum(is.finite(plot_df[[m]])), integer(1))
  if (orient == "concern") labels[metrics %in% lower_is_concerning] <-
    paste0(labels[metrics %in% lower_is_concerning], "\n(inverted)")
  labels <- paste0(labels, "\nn = ", counts, "/", nrow(rec))
  duplicate_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate_labels] <- paste0(labels[duplicate_labels], "\n", metrics[duplicate_labels])
  long_df$metric_label <- factor(long_df$metric, levels = metrics, labels = labels)

  ## ---- colours ------------------------------------------------------------
  ## Palette keys must match the labels actually present, including the derived
  ## ones. Anything unmatched falls back to grey and is legended as such.
  class_cols <- c(
    "not flagged"         = "#555555",
    "returned"            = "#2E7D32",
    "incomplete return"   = "#B2182B",
    "displaced plateau"   = "#7B3294",
    ## keys used by user-supplied classifications
    fast_recovery         = "#2E7D32",
    slow_recovery         = "#8C6D31",
    overshoot             = "#2166AC",
    hysteresis            = "#7B3294",
    incomplete_recovery   = "#B2182B",
    unclassified          = "grey50"
  )

  ## Dark tiles need light text: fixed dark text is unreadable at the top of
  ## the ramp.
  long_df$.label_col <- ifelse(
    !is.na(long_df$value_scaled) & long_df$value_scaled > 0.62,
    "white", "#111111"
  )

  subtitle <- paste0(
    "Ordered by ", order_by, "; colour = within-column scaled magnitude",
    if (orient == "concern")
      "; rate inverted, overshoot neutral" else "; darker is higher-valued"
  )

  ggplot2::ggplot(
    long_df,
    ggplot2::aes(x = .data$metric_label, y = .data$.trajectory,
                 fill = .data$value_scaled)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45,
                       width = 0.96, height = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(.data$value), "-",
                                  signif(.data$value, 2)),
                   colour = .data$.label_col),
      size = base_size / 4, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = unique(long_df[, c(".trajectory", "trajectory_class")]),
      ggplot2::aes(x = 0.35, y = .data$.trajectory,
                   colour = .data$trajectory_class),
      inherit.aes = FALSE, size = 3.2
    ) +
    ggplot2::scale_fill_gradientn(
      colours   = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits    = c(0, 1),
      na.value  = "grey93",
      name      = "Scaled\n(within column)"
    ) +
    ggplot2::scale_colour_manual(
      values   = c(class_cols, "white" = "white", "#111111" = "#111111"),
      breaks   = intersect(names(class_cols),
                           unique(long_df$trajectory_class)),
      na.value = "grey60",
      name     = if (derived_class) "Trajectory\n(derived)" else "Trajectory\nclass"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title    = "Redox resilience recovery landscape",
      subtitle = subtitle,
      caption = paste0("Numbers are raw values; grey / dash = unavailable; n = finite trajectories.\n",
        if (derived_class) "Not flagged is not evidence of return or equivalence. " else "",
        if (length(dropped)) paste0("Omitted: ", paste(dropped, collapse = ", "), ".") else ""),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      ## Bold at full base_size is what makes six headers collide; 0.78x
      ## plain with a little lineheight fits comfortably.
      axis.text.x = ggplot2::element_text(colour = "#111111",
                                          size = base_size * 0.78,
                                          lineheight = 0.95,
                                          margin = ggplot2::margin(t = 6)),
      axis.text.y = ggplot2::element_text(colour = "#222222",
                                          size = base_size * 0.68),
      plot.title  = ggplot2::element_text(face = "bold", size = base_size + 2,
                                          colour = "#111111"),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.85,
                                            colour = "#444444",
                                            margin = ggplot2::margin(b = 12)),
      legend.title    = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.caption = ggplot2::element_text(size = base_size * 0.65, hjust = 0),
      plot.margin     = ggplot2::margin(15, 20, 15, 34)
    )
}
