#' Sensitivity to domain aggregation weights
#' @param res RRI result.
#' @param weight_grid Plant weights in (0,1), or a data frame/matrix with
#' named Physio, Soil, Micro columns specifying complete alternative weights.
#' @return Alternative normalized weights, finite-pair count and Spearman correlation.
#' This conditions on the already computed features, reductions and missingness.
#' @importFrom stats cor
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#'   rri_sensitivity(res)
#' }
#' @export
rri_sensitivity <- function(res, weight_grid = seq(0.2, 0.6, by = 0.1)) {
  rs <- res$row_scores; domains <- c("Physio", "Soil", "Micro")
  if (!all(c(domains, "RRI") %in% names(rs))) stop("Required score columns absent.")
  if (is.numeric(weight_grid) && is.null(dim(weight_grid))) {
    if (!length(weight_grid) || any(!is.finite(weight_grid)) || any(weight_grid <= 0 | weight_grid >= 1))
      stop("Vector weights must lie in (0,1).")
    weights <- cbind(Physio = weight_grid, Soil = (1 - weight_grid)/2, Micro = (1 - weight_grid)/2)
  } else {
    weights <- as.matrix(weight_grid)
    if (!is.numeric(weights) || nrow(weights) < 1L)
      stop("Weight alternatives must be a nonempty numeric matrix or data frame.")
    if (!all(domains %in% colnames(weights))) stop("Use Physio, Soil, Micro weight columns.")
    weights <- weights[, domains, drop = FALSE]
  }
  rows <- lapply(seq_len(nrow(weights)), function(i) {
    w <- .rri_weights(weights[i, ])
    alt <- .rri_weighted(rs[, domains], w)$score
    data.frame(weight_physio = w[1], weight_soil = w[2], weight_micro = w[3],
      n_pairs = sum(is.finite(alt) & is.finite(rs$RRI)),
      spearman_rank_correlation = .rri_cor(rs$RRI, alt, "spearman"))
  })
  do.call(rbind, rows)
}
