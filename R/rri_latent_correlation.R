#' @title Correlation with a Simulator-Defined Target
#'
#' @description
#' Computes a descriptive correlation between per-sample RRI and a prescribed
#' simulator target. This is an internal simulation benchmark and is not
#' predictive accuracy, empirical validation or recovery of a true latent state.
#'
#' @param res An object returned by \code{rri_pipeline_st()}.
#' @param latent_truth Numeric simulator-target vector. The legacy argument name
#' is retained for compatibility.
#' @param method Correlation method. One of \code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"}.
#'
#' @return A single numeric correlation coefficient.
#'
#' @details
#' This function is designed for simulation benchmarking. In empirical
#' datasets, no known latent state exists and this metric should not be used.
#'
#' @importFrom stats cor
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline_st(sim$ROS_flux, sim$Eh_stability)
#'   rri_latent_correlation(res, sim$latent_truth)
#' }
#' @export
rri_latent_correlation <- function(res,
                                   latent_truth,
                                   method = c("pearson", "spearman", "kendall")) {
  method <- match.arg(method)

  if (is.null(res$row_scores$RRI)) {
    stop("res must contain a row_scores$RRI column.", call. = FALSE)
  }

  rri <- res$row_scores$RRI
  if (!is.numeric(rri) || !is.numeric(latent_truth))
    stop("RRI and latent_truth must be numeric vectors.", call. = FALSE)

  if (length(rri) != length(latent_truth)) {
    stop("RRI and latent_truth must have equal length.", call. = FALSE)
  }

  .rri_cor(rri, latent_truth, method = method)
}
