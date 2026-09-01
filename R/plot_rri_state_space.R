#' Plot domain-score space with correctly matched trajectory diagnostics
#'
#' @param res RRI result with aligned identifiers in row_scores.
#' @param rec Optional one-row-per-trajectory diagnostic table.
#' @param x_property,y_property Domain scores, or group association/kinetics.
#' @param colour_by RRI, Memory, or an explicitly supplied trajectory_class.
#' @param group_cols Full key shared between row_scores and rec.
#' @param base_size Plot font size.
#' @return ggplot. Domain scores are not relabelled as mechanistic properties.
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic scale_colour_viridis_c
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
#'                       id = sim$id)
#'   plot_rri_state_space(res, group_cols = c("plot", "depth", "plant_id"))
#' }
#' @export
plot_rri_state_space <- function(res, rec = NULL,
  x_property = c("Physio", "Connectivity", "Soil", "Micro"),
  y_property = c("Soil", "Micro", "Physio", "Kinetics"),
  colour_by = c("RRI", "Memory", "trajectory_class"),
  group_cols = c("plot", "depth", "plant_id"), base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Requires ggplot2.")
  x_property <- match.arg(x_property); y_property <- match.arg(y_property)
  colour_by <- match.arg(colour_by); rs <- as.data.frame(res$row_scores)
  match_group <- function(d, column) {
    if (!column %in% names(d)) stop("Missing diagnostic column: ", column)
    a <- .rri_key(rs, group_cols); b <- .rri_key(d, group_cols)
    if (anyDuplicated(b)) stop("Diagnostics must have one row per full group key.")
    ix <- match(a, b)
    if (anyNA(ix)) warning("Some plotted trajectories have no matching diagnostic.")
    d[[column]][ix]
  }
  if (x_property == "Connectivity") {
    cdat <- rri_connectivity_score(res, per_group = TRUE, group_cols = group_cols)$connectivity_score
    rs$.x <- match_group(cdat, "connectivity_score")
    x_label <- "Within-group association magnitude"
  } else {
    if (!x_property %in% names(rs)) stop("Requested x score absent.")
    rs$.x <- rs[[x_property]]; x_label <- paste(x_property, "score")
  }
  if (y_property == "Kinetics") {
    if (is.null(rec)) stop("Provide rec for Kinetics.")
    if (!"kinetics_score" %in% names(rec)) rec <- rri_kinetics_score(rec)
    rs$.y <- match_group(rec, "kinetics_score"); y_label <- "Recovery speed descriptor"
  } else {
    if (!y_property %in% names(rs)) stop("Requested y score absent.")
    rs$.y <- rs[[y_property]]; y_label <- paste(y_property, "score")
  }
  if (colour_by == "RRI") rs$.colour <- rs$RRI else {
    if (is.null(rec)) stop("Provide rec for trajectory-level colouring.")
    if (colour_by == "Memory") {
      if (!"memory_index" %in% names(rec)) rec <- rri_memory_index(rec)
      rs$.colour <- match_group(rec, "memory_index")
    } else rs$.colour <- as.factor(match_group(rec, "trajectory_class"))
  }
  p <- ggplot2::ggplot(rs, ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$.colour)) +
    ggplot2::geom_point(size = 2.2, alpha = 0.7) +
    ggplot2::labs(title = "Domain-score space", x = x_label, y = y_label,
      colour = colour_by, subtitle = "One point per sample; trajectory diagnostics matched by identifiers") +
    ggplot2::theme_classic(base_size = base_size)
  if (colour_by == "trajectory_class") p else p + ggplot2::scale_colour_viridis_c()
}
