#' @title HRRI: Diagnostics for Soil-Plant-Microbial Redox Recovery
#' @keywords internal
#' @aliases HRRI-package
#' @description
#' HRRI keeps two kinds of quantity apart. Mechanistic quantities --
#' event-window accessible electron capacity and stoichiometric oxygen demand --
#' are computed from inventories, accessibility and exchange rates that the user
#' declares, and are never estimated from the scores. Observational quantities --
#' domain scores and recovery descriptors built from aligned soil, plant and
#' microbial measurements -- carry their orientation, coverage and estimability
#' with them, so an unsupported value is reported as missing rather than
#' defaulted.
#'
#' @section Where to start:
#' \itemize{
#'   \item \code{vignette("HRRI_workflow", package = "HRRI")} for an end-to-end run.
#'   \item \code{vignette("HRRI_gallery", package = "HRRI")} for how to read the figures.
#'   \item \code{vignette("HRRI_paper_figures", package = "HRRI")} to reproduce the
#'         six manuscript figures.
#' }
#'
#' @section What the package does not claim:
#' Simulated demonstrations verify accounting and analytical behaviour. They are
#' not environmental validation, and a statistically identifiable
#' \eqn{\alpha} or \eqn{k} is not a mechanistically identified one.
"_PACKAGE"

