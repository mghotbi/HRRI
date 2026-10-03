#' @title EMS plotting theme
#'
#' @description A simple ggplot theme used for HRRI visualizations.
#'
#' @param base_size Base font size
#' @return A ggplot2 theme object
#' @importFrom ggplot2 theme_classic theme element_text element_blank
#' @examples
#' \donttest{
#'   library(ggplot2)
#'   ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
#'     geom_point() + theme_ems()
#' }
#' @export
theme_ems <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#233441"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1,
        margin = ggplot2::margin(b = 9)),
      plot.subtitle = ggplot2::element_text(size = base_size - 1,
        colour = "#45545e", margin = ggplot2::margin(b = 7)),
      axis.title = ggplot2::element_text(face = "plain", colour = "#233441"),
      axis.text = ggplot2::element_text(colour = "#45545e", size = base_size - 1),
      axis.line = ggplot2::element_line(colour = "#67737c", linewidth = .35),
      axis.ticks = ggplot2::element_line(colour = "#67737c", linewidth = .3),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#edf0f2", linewidth = .3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = base_size - 2, hjust = 0),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"))
}
