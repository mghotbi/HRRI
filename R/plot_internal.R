# Internal shared plotting styles; not evidence or scores.
.OI <- list(                        # Okabe-Ito accessible palette
  blue   = "#0072B2",
  red    = "#D55E00",
  green  = "#009E73",
  purple = "#CC79A7",
  orange = "#E69F00",
  sky    = "#56B4E9",
  yellow = "#F0E442",
  black  = "#000000",
  gray   = "#555555"
)
.OI_vec <- unlist(.OI[c("blue","red","green","purple","orange","sky","yellow")])

.theme_hrri <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1,
                                             hjust = 0, margin = ggplot2::margin(b = 3)),
    plot.subtitle    = ggplot2::element_text(size = base_size - 1, colour = "#555555",
                                             hjust = 0, margin = ggplot2::margin(b = 4)),
    axis.title       = ggplot2::element_text(size = base_size),
    axis.text        = ggplot2::element_text(size = base_size - 1),
    legend.title     = ggplot2::element_text(size = base_size - 1, face = "bold"),
    legend.text      = ggplot2::element_text(size = base_size - 1),
    legend.key.size  = grid::unit(3, "mm"),
    panel.grid.major = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.3),
    plot.margin      = ggplot2::margin(4, 4, 4, 4, "mm"),
    strip.background = ggplot2::element_blank(),
    strip.text       = ggplot2::element_text(face = "bold", size = base_size - 1)
  )
}

.safe_cor <- function(x, y) .rri_cor(x, y)

.safe_rmse <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 1) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}
