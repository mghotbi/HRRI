#' Score departures from an explicitly defined reference
#'
#' @description An optional, transparent alternative to latent-axis scoring.
#' Each feature score is max(0, 1 - abs(value - target) / tolerance).
#' The result measures proximity to the declared reference, not validated
#' ecosystem functioning or a universal resilience scale.
#' @param data Numeric feature data frame; rows must be aligned with id.
#' @param reference Data frame with feature, domain, target, tolerance, weight.
#' Domains are Physio, Soil or Micro. Tolerance is a positive distance from
#' target at which the feature score reaches zero. Reference rows for unmeasured
#' features may be retained to report coverage against a common panel.
#' @param id Optional aligned identifiers.
#' @param domain_weights Named, non-negative domain weights.
#' @param min_coverage Minimum weighted within-domain feature coverage.
#' @param na_policy Use available domains or require every positive-weight domain.
#' @return An RRI object with fixed-reference domain scores, feature scores,
#' coverage, effective row-specific domain weights, and reference metadata.
#' @importFrom stats setNames complete.cases
#' @examples
#' dat <- data.frame(Eh = c(100, 200, 150), pH = c(5.5, 6.0, 5.8))
#' ref <- data.frame(feature = c("Eh", "pH"), domain = c("Soil", "Soil"),
#'                   target = c(150, 5.8), tolerance = c(100, 0.5), weight = 1)
#' rri_reference_scores(dat, ref)$row_scores
#' @export
rri_reference_scores <- function(data, reference, id = NULL,
  domain_weights = c(Physio = 0.4, Soil = 0.35, Micro = 0.25),
  min_coverage = 0.5, na_policy = c("available", "complete")) {
  na_policy <- match.arg(na_policy)
  data <- .rri_numeric_df(data, "data")
  reference <- as.data.frame(reference)
  domains <- c("Physio", "Soil", "Micro")
  needed <- c("feature", "domain", "target", "tolerance", "weight")
  if (!nrow(data) || !all(needed %in% names(reference)) || !nrow(reference))
    stop("Provide nonempty data and a complete reference table.", call. = FALSE)
  if (anyNA(reference[, needed]) || anyDuplicated(reference$feature) ||
      any(!reference$domain %in% domains) || any(!nzchar(reference$feature)))
    stop("Reference features must be unique and domains valid.", call. = FALSE)
  if (!is.numeric(reference$target) || any(!is.finite(reference$target)) ||
      !is.numeric(reference$tolerance) || any(!is.finite(reference$tolerance)) ||
      any(reference$tolerance <= 0)) stop("Invalid reference targets or tolerances.", call. = FALSE)
  .rri_weights(reference$weight)
  if (length(min_coverage) != 1L || !is.finite(min_coverage) ||
      min_coverage < 0 || min_coverage > 1) stop("min_coverage must be in [0,1].")
  if (is.null(names(domain_weights)) || anyDuplicated(names(domain_weights)) ||
      any(!names(domain_weights) %in% domains)) stop("Use named domain weights.")
  .rri_weights(domain_weights)
  dw <- stats::setNames(rep(0, 3), domains)
  dw[names(domain_weights)] <- domain_weights
  if (!is.null(id) && nrow(as.data.frame(id)) != nrow(data)) stop("id row count differs.")
  n <- nrow(data)
  features <- matrix(NA_real_, n, nrow(reference), dimnames = list(NULL, reference$feature))
  for (j in seq_len(nrow(reference))) {
    f <- reference$feature[j]
    if (!f %in% names(data)) next
    x <- data[[f]]
    if (!is.numeric(x)) stop("Feature must be numeric: ", f)
    x[!is.finite(x)] <- NA_real_
    features[, j] <- pmax(0, 1 - abs(x - reference$target[j]) / reference$tolerance[j])
  }
  scores <- coverage <- matrix(NA_real_, n, 3, dimnames = list(NULL, domains))
  for (d in domains) {
    j <- which(reference$domain == d & reference$weight > 0)
    if (!length(j)) next
    z <- .rri_weighted(features[, j, drop = FALSE], reference$weight[j])
    scores[, d] <- z$score
    scores[z$coverage < min_coverage, d] <- NA_real_
    coverage[, d] <- z$coverage
  }
  combined <- .rri_weighted(scores, dw)
  if (na_policy == "complete") combined$score[combined$coverage < 1 - 1e-12] <- NA_real_
  eff <- sweep(is.finite(scores), 2, dw, "*")
  denom <- rowSums(eff)
  eff[denom > 0, ] <- eff[denom > 0, , drop = FALSE] / denom[denom > 0]
  eff[denom == 0, ] <- NA_real_
  rs <- data.frame(scores, RRI = combined$score, domain_coverage = combined$coverage,
                   n_domains = combined$n_observed)
  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (any(names(id) %in% names(rs))) stop("id names conflict with score columns.")
    rs <- cbind(id, rs)
  }
  comp <- scores / rowSums(scores)
  comp[!stats::complete.cases(scores) | rowSums(scores, na.rm = TRUE) == 0, ] <- NA_real_
  structure(list(row_scores = rs, row_scores_comp = data.frame(comp, RRI = rs$RRI),
    feature_scores = features, coverage = data.frame(coverage), effective_weights = eff,
    dyn_scores = NULL, meta = list(method = "reference_proximity", reference = reference,
      domain_weights = dw, min_coverage = min_coverage, na_policy = na_policy,
      interpretation = "Proximity to a declared reference; domain subsets change the estimand.")),
    class = "RRI")
}
