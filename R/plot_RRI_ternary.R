#' @title Ternary Plot of Relative Domain Scores
#'
#' @description
#' Creates a ternary diagram of the relative magnitudes of plant, soil and
#' microbial domain scores after closure to a unit sum. These coordinates are
#' display quantities, not fractions of causal buffering capacity. Points are
#' filled according to the corresponding composite \code{RRI} value.
#'
#' @param ternary_df A data frame containing compositional columns
#'   \code{Physio}, \code{Soil}, \code{Micro}, and \code{RRI}.
#' @param point_size Numeric; size of ternary points.
#' @param point_alpha Numeric between 0 and 1 controlling point transparency.
#' @param palette Character; viridis palette option.
#' @param show_subtitle Logical; display system-level RRI mean in subtitle.
#' @param show_centroid Logical; add compositional centroid marker.
#' @param centroid_shape Numeric; ggplot2 shape for centroid marker.
#' @param centroid_size Numeric multiplier for centroid size.
#' @param tolerance Numeric; tolerance used for compositional closure checks.
#' @param renormalize Logical; if TRUE, renormalises rows that do not sum to 1.
#' @param centroid_method Character; one of \code{"auto"},
#'   \code{"simplex_mean"}, or \code{"aitchison_mean"}.
#'
#' @details
#' Closure removes absolute score magnitude: rows with proportional domain scores
#' occupy the same position even when their composite scores differ. Do not infer
#' mechanistic allocation, causal contribution or substitution from this plot.
#' If clr-transformed coordinates are attached as an attribute (\code{"clr"}),
#' the centroid can be computed using the Aitchison mean. Otherwise, a simplex
#' arithmetic mean is used.
#'
#' @return A \code{ggtern} object.
#'
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' # Optional packages ggtern and viridis are required.
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2,
#'   n_depth = 2,
#'   n_plant = 2,
#'   n_time = 8,
#'   p_micro = 6,
#'   seed = 1234
#' )
#'
#' # ---- Compute RedoxRRI ----
#' res <- rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   reducer = "per_domain",
#'   scaling = "pnorm"
#' )
#'
#' # ---- Extract compositional scores ----
#' ternary_df <- res$row_scores_comp
#'
#' # ---- Plot ternary allocation ----
#' p <- plot_RRI_ternary(
#'   ternary_df,
#'   point_size = 3,
#'   show_centroid = TRUE
#' )
#' }
#'
#' @export
plot_RRI_ternary <- function(
  ternary_df,
  point_size = 5,
  point_alpha = 0.9,
  palette = "plasma",
  show_subtitle = TRUE,
  show_centroid = TRUE,
  centroid_shape = 23,
  centroid_size = 1.4,
  tolerance = 1e-6,
  renormalize = FALSE,
  centroid_method = c("auto", "simplex_mean", "aitchison_mean")
) {
  centroid_method <- match.arg(centroid_method)

  # ---- dependency checks ----
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {ggtern} package.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {ggplot2} package.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {viridis} package.")
  }

  ternary_df <- as.data.frame(ternary_df)

  required <- c("Physio", "Soil", "Micro", "RRI")
  missing <- setdiff(required, names(ternary_df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  ternary_df[, required] <- .rri_numeric_df(ternary_df[, required, drop = FALSE])

  if (any(as.matrix(ternary_df[, c("Physio", "Soil", "Micro")]) < 0, na.rm = TRUE))
    stop("Ternary components cannot be negative.")
  # ---- remove invalid rows ----
  comps <- ternary_df[, c("Physio", "Soil", "Micro")]
  rs <- rowSums(comps)

  keep <- is.finite(rs) &
    rs > tolerance &
    stats::complete.cases(comps) &
    is.finite(ternary_df$RRI)

  ternary_df <- ternary_df[keep, , drop = FALSE]

  if (nrow(ternary_df) == 0) {
    stop("No valid compositional rows available for ternary plotting.")
  }

  # ---- closure validation ----
  comps <- ternary_df[, c("Physio", "Soil", "Micro")]
  rs2 <- rowSums(comps)

  if (any(abs(rs2 - 1) > tolerance)) {
    if (isTRUE(renormalize)) {
      rs2[rs2 == 0] <- NA_real_
      comps <- comps / rs2
      ternary_df[, c("Physio", "Soil", "Micro")] <- comps
    } else {
      stop("Physio, Soil, and Micro must sum to 1 within tolerance; set renormalize = TRUE to override.")
    }
  }

  # ---- subtitle ----
  rri_index <- attr(ternary_df, "RRI_index")
  if (is.null(rri_index) || !is.finite(rri_index)) {
    rri_index <- mean(ternary_df$RRI, na.rm = TRUE)
  }

  subtitle_text <- if (isTRUE(show_subtitle)) {
    sprintf("Mean composite score: %.3f", rri_index)
  } else {
    NULL
  }

  # ---- centroid ----
  centroid <- NULL

  if (isTRUE(show_centroid)) {
    use_aitchison <- FALSE
    if (centroid_method == "aitchison_mean") use_aitchison <- TRUE
    if (centroid_method == "auto" && !is.null(attr(ternary_df, "clr"))) use_aitchison <- TRUE

    if (use_aitchison) {
      X <- as.matrix(ternary_df[, c("Physio", "Soil", "Micro")])
      if (any(X <= 0)) stop("Aitchison centroid requires positive components; declare zero handling first.")
      z <- exp(colMeans(log(X)))
      centroid <- as.data.frame(as.list(z / sum(z)))
    } else {
      centroid <- data.frame(
        Physio = mean(ternary_df$Physio, na.rm = TRUE),
        Soil   = mean(ternary_df$Soil, na.rm = TRUE),
        Micro  = mean(ternary_df$Micro, na.rm = TRUE)
      )
      centroid <- centroid / sum(centroid)
    }
  }

  # ---- base plot ----
  p <- ggtern::ggtern(
    data = ternary_df,
    ggplot2::aes(
      x = .data$Physio,
      y = .data$Soil,
      z = .data$Micro,
      fill = .data$RRI
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      size = point_size,
      alpha = point_alpha,
      colour = "black",
      stroke = 0.5
    ) +
    viridis::scale_fill_viridis(option = palette, direction = -1) +
    ggplot2::labs(
      title = "Relative Domain-Score Geometry",
      subtitle = subtitle_text,
      x = "Plant physiology",
      y = "Soil redox chemistry",
      z = "Microbial score",
      fill = "RRI"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggtern::theme_nomask() +
    ggtern::theme_showarrows() +
    ggtern::theme_clockwise()

  if (isTRUE(show_centroid) && !is.null(centroid)) {
    p <- p +
      ggplot2::geom_point(
        data = centroid,
        ggplot2::aes(x = .data$Physio, y = .data$Soil, z = .data$Micro),
        inherit.aes = FALSE,
        shape = centroid_shape,
        size = point_size * centroid_size,
        fill = "white",
        colour = "black",
        stroke = 1.2
      )
  }

  return(p)
}
