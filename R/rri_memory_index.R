#' Persistent-displacement and loop-area diagnostic
#'
#' @description A descriptive composite, not proof of ecological memory. Slow
#' relaxation, baseline drift and continuing forcing can also produce displacement.
#' @param rec Recovery metric table.
#' @param H_weight,I_weight,lag_weight Non-negative component weights.
#' @param normalise_inputs FALSE uses bounded dimensionless H and I fractions;
#' TRUE requests cohort-relative min-max scaling. Lag requires TRUE.
#' @return Input with memory_index, memory_coverage and heuristic memory_class.
#' @examples
#' rec <- data.frame(H_hysteresis = c(0.1, 0.3, 0.05),
#'                   incomplete_return_frac = c(0.2, 0.4, 0.1))
#' rri_memory_index(rec)$memory_index
#' @export
rri_memory_index <- function(rec, H_weight = 0.50, I_weight = 0.50,
  lag_weight = 0, normalise_inputs = FALSE) {
  rec <- as.data.frame(rec); n <- nrow(rec)
  w <- .rri_weights(c(H_weight, I_weight, lag_weight))
  if (lag_weight > 0 && !normalise_inputs) stop("Lag has time units; use normalise_inputs=TRUE or zero lag weight.")
  resolve <- function(candidates) {
    nm <- candidates[candidates %in% names(rec)]
    if (!length(nm)) return(rep(NA_real_, n))
    x <- rec[[nm[1]]]
    if (!is.numeric(x)) stop("Memory inputs must be numeric.")
    x[!is.finite(x)] <- NA_real_; abs(x)
  }
  h <- resolve(c("H_hysteresis", "H", "H_abs", "H_norm"))
  if ("H_axis" %in% names(rec)) h[is.na(rec$H_axis) | rec$H_axis != "forcing"] <- NA_real_
  i <- resolve(c("incomplete_return_frac", "I_norm"))
  if (all(is.na(i)) && normalise_inputs) i <- resolve(c("incomplete_return", "I"))
  lag <- resolve("tau_lag")
  mat <- cbind(h, i, lag)
  if (normalise_inputs) {
    for (j in seq_len(ncol(mat))) mat[, j] <- .rri_scale(mat[, j])
  } else {
    mat[] <- pmin(1, mat)
  }
  z <- .rri_weighted(mat, w)
  rec$memory_index <- z$score; rec$memory_coverage <- z$coverage
  rec$memory_class <- cut(z$score, c(-Inf, 1/3, 2/3, Inf),
    c("low_memory", "moderate_memory", "high_memory"))
  attr(rec, "memory_interpretation") <- "Descriptive displacement; class cutoffs are not evidence of distinct ecological states."
  rec
}
