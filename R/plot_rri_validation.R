#' Scatter of HRRI score against a simulator-defined target
#'
#' @description Plots mean RRI versus a declared synthetic target per aggregate row, annotated
#'   with descriptive Pearson r and direct score-target RMSE, not LOO error.
#' No model is trained or held out by this plotting function.
#'
#' @param pool_agg  Seed-level aggregate data frame with columns \code{RRI} and
#'   \code{truth} (one row per seed).
#' @param rri_col   Name of the HRRI column in pool_agg (default "RRI").
#' @param truth_col Name of the latent truth column (default "truth").
#' @param colour_col Optional column for point colour (e.g., "n_cycles").
#' @param label_col  Optional column for point labels.
#' @param base_size  Base font size.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_abline geom_line geom_point geom_text annotate labs coord_equal scale_fill_viridis_c scale_fill_viridis_d
#' @importFrom stats predict lm
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#'   agg <- data.frame(RRI = res$row_scores$RRI, truth = sim$latent_truth)
#'   plot_rri_validation(agg)
#' }
#' @export
plot_rri_validation <- function(
  pool_agg,
  rri_col    = "RRI",
  truth_col  = "truth",
  colour_col = NULL,
  label_col  = NULL,
  base_size  = 9
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Requires ggplot2.")

  df <- as.data.frame(pool_agg)
  for (col in c(rri_col, truth_col))
    if (!col %in% names(df)) stop("Missing column: ", col)

  if (!is.numeric(df[[truth_col]]) || !is.numeric(df[[rri_col]]))
    stop("Scores and targets must be numeric.")
  df$.x <- df[[truth_col]]
  df$.y <- df[[rri_col]]
  ok    <- is.finite(df$.x) & is.finite(df$.y)
  df_ok <- df[ok, , drop = FALSE]
  if (!nrow(df_ok)) stop("No finite score-target pairs.")

  r    <- .safe_cor(df_ok$.x, df_ok$.y)
  rmse <- .safe_rmse(df_ok$.y, df_ok$.x)
  ann <- sprintf("Descriptive r = %.3f | Direct RMSE = %.3f", r, rmse)
  lm_df <- NULL
  if (is.finite(r)) {
    xr <- range(df_ok$.x)
    lm_df <- data.frame(.x = seq(xr[1], xr[2], length.out = 100))
    lm_df$.y <- stats::predict(stats::lm(.y ~ .x, data = df_ok), lm_df)
  }

  p <- ggplot2::ggplot(df_ok, ggplot2::aes(x = .data$.x, y = .data$.y))

  # 1:1 reference
  p <- p + ggplot2::geom_abline(intercept = 0, slope = 1,
                                 colour = "#AAAAAA", linewidth = 0.4,
                                 linetype = "longdash")

  # OLS regression line
  if (!is.null(lm_df)) p <- p + ggplot2::geom_line(data = lm_df,
                               ggplot2::aes(x = .data$.x, y = .data$.y),
                               colour = .OI$blue, linewidth = 0.8, inherit.aes = FALSE)

  # Points
  if (!is.null(colour_col) && colour_col %in% names(df_ok)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = .data[[colour_col]]),
      shape = 21, size = 2.5, colour = "white", stroke = 0.4
    )
    p <- p + if (is.numeric(df_ok[[colour_col]]))
      ggplot2::scale_fill_viridis_c(option = "plasma", name = colour_col) else
      ggplot2::scale_fill_viridis_d(option = "plasma", name = colour_col)
  } else {
    p <- p + ggplot2::geom_point(
      shape = 21, size = 2.5, fill = .OI$blue,
      colour = "white", stroke = 0.4, alpha = 0.85
    )
  }

  # Optional labels
  if (!is.null(label_col) && label_col %in% names(df_ok)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[[label_col]]),
      size = base_size * 0.25, vjust = -0.8, colour = "#444444"
    )
  }

  p <- p +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = ann,
                      hjust = -0.05, vjust = 1.4,
                      size = base_size * 0.28, colour = "#333333") +
    ggplot2::labs(
      title    = "HRRI diagnostic agreement",
      subtitle = sprintf("n = %d aggregate rows; dashed grey = 1:1 line; blue = OLS", nrow(df_ok)),
      x        = "Simulator-defined target",
      y        = "Observation-derived HRRI score"
    ) +
    ggplot2::coord_equal() +
    .theme_hrri(base_size)

  p
}
