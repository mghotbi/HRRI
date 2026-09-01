#' Construct an explicitly weighted microbial guild contrast
#' @description Summarises supplied guild measurements or proxies. The weights
#' define a contrast, not a universal ordering of microbial resilience.
#' Denitrification, sulfate reduction and methanogenesis may be beneficial or
#' detrimental depending on the specified ecosystem function and disturbance.
#' @param micro_traits Numeric data frame with comparable, justified guild scales.
#' Gene abundance or expression does not by itself measure process rate.
#' @param weights Named finite signed weights. Default legacy weights are
#' illustrative only and trigger a warning; supply scientifically justified weights.
#' @param scale If TRUE, scale finite contrasts within the supplied cohort to `[0, 1]`. Constants map to 0.5 and wholly unobserved rows remain NA.
#' @return Numeric vector with coverage and raw_contrast attributes. The raw
#' contrast is a signed weighted sum divided by the available absolute weight.
#' Missing guilds are not zeros; changing availability changes the estimand.
#' @seealso rri_reference_scores
#' @examples
#' x <- data.frame(EET_reduction = c(0.2, 0.4, NA),
#'                 methanogenesis = c(0.3, 0.1, NA))
#' rri_micro_functional_score(x, weights = c(EET_reduction = 1, methanogenesis = -1))
#' @export
rri_micro_functional_score <- function(micro_traits, weights = NULL, scale = TRUE) {
  # Assisted-by: OpenAI Codex. Distinguish missing observations from measured zeros.
  mt <- .rri_numeric_df(micro_traits, "micro_traits")
  if (is.null(weights)) {
    warning("Default guild weights are illustrative, not validated resilience weights.",
            call. = FALSE)
    weights <- c(EET_reduction = 0.25, Fe_oxidation = 0.20, Mn_oxidation = 0.15,
      denitrification = 0.15, nitrification = 0.10, DNRA = 0.05,
      methane_oxidation = 0.05, sulfate_reduction = -0.10, methanogenesis = -0.20)
  }
  if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) ||
      anyNA(names(weights)) || any(!nzchar(names(weights))) ||
      anyDuplicated(names(weights)) || any(!is.finite(weights)) || sum(abs(weights)) == 0)
    stop("weights must be uniquely named, finite and not all zero.")
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) stop("scale must be TRUE or FALSE.")
  x <- matrix(NA_real_, nrow(mt), length(weights), dimnames = list(NULL, names(weights)))
  present <- intersect(names(weights), names(mt))
  x[, present] <- as.matrix(mt[, present, drop = FALSE])
  w <- abs(weights) / sum(abs(weights))
  signed <- sweep(x, 2L, sign(weights), "*")
  z <- .rri_weighted(signed, w)
  score <- if (scale) .rri_scale(z$score) else z$score
  attr(score, "coverage") <- z$coverage
  attr(score, "raw_contrast") <- z$score
  score
}
