#' Aligned time series with separate physical units
#' @description Shows Eh, EAC and the observation-derived score in separate
#' panels sharing time. It does not place unlike units on a common axis.
#' @param sim Simulator output containing id and soil_data.
#' @param res Pipeline output containing row_scores.
#' @param plot_id,depth_id,plant_id Identifiers for a single trajectory.
#' @param perturb_start,perturb_end Optional disturbance interval in input time units.
#' @param base_size Base font size.
#' @return A ggplot with three vertically aligned panels.
#' @importFrom ggplot2 ggplot aes annotate geom_line geom_point facet_wrap scale_colour_manual labs
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
#'                       id = sim$id)
#'   plot_rri_timeseries(sim, res)
#' }
#' @export
plot_rri_timeseries <- function(sim, res, plot_id = "P1", depth_id = "D1",
  plant_id = "Plant1", perturb_start = NULL, perturb_end = NULL, base_size = 9) {
  # Assisted-by: OpenAI Codex. Preserve units, identifiers and missing observations.
  id <- as.data.frame(sim$id)
  keys <- c("plot", "depth", "plant_id", "time")
  if (!all(keys %in% names(id))) stop("sim$id requires plot, depth, plant_id and time.")
  sdf <- as.data.frame(sim$soil_data)
  if (nrow(sdf) != nrow(id)) stop("Soil observations are not aligned with id.")
  if (!"Eh" %in% names(sdf)) stop("soil_data requires Eh.")
  rs <- .rri_align_scores(res$row_scores, id, keys)
  if (!"RRI" %in% names(rs)) stop("row_scores requires RRI.")
  idx <- which(id$plot == plot_id & id$depth == depth_id & id$plant_id == plant_id)
  if (!length(idx)) stop("No matching trajectory.")
  if (!is.numeric(id$time) || any(!is.finite(id$time[idx]))) stop("Time must be finite numeric.")
  idx <- idx[order(id$time[idx])]
  if (anyDuplicated(id$time[idx])) stop("Duplicate times in selected trajectory.")
  vals <- list(sdf$Eh[idx],
    if ("EAC" %in% names(sdf)) sdf$EAC[idx] else rep(NA_real_, length(idx)),
    rs$RRI[idx])
  if (!all(vapply(vals, is.numeric, logical(1)))) stop("Plotted measurements must be numeric.")
  labels <- c("Eh (mV)", "EAC (mmol electrons / kg dry soil)", "HRRI score (dimensionless)")
  df <- data.frame(time = rep(id$time[idx], 3L), value = unlist(vals),
    variable = factor(rep(labels, each = length(idx)), levels = labels))
  df$value[!is.finite(df$value)] <- NA_real_
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                        colour = .data$variable))
  if (xor(is.null(perturb_start), is.null(perturb_end)))
    stop("Supply both disturbance interval endpoints.")
  if (!is.null(perturb_start)) {
    interval <- c(perturb_start, perturb_end)
    if (!is.numeric(interval) || length(interval) != 2L ||
        any(!is.finite(interval)) || interval[1] >= interval[2])
      stop("Disturbance endpoints must be finite and increasing.")
    p <- p + ggplot2::annotate("rect", xmin = interval[1], xmax = interval[2],
      ymin = -Inf, ymax = Inf, fill = .OI$orange, alpha = 0.12)
  }
  p + ggplot2::geom_line(linewidth = 0.65, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
    ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y", drop = FALSE) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c(.OI$blue, .OI$green, .OI$red), labels), guide = "none") +
    ggplot2::labs(x = "Time (input units; simulator uses days)", y = NULL,
      title = "Soil observations and HRRI",
      subtitle = paste(plot_id, depth_id, plant_id, sep = " / ")) +
    .theme_hrri(base_size)
}
