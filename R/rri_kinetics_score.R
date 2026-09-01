#' Descriptive recovery speed score
#'
#' @description Combines response lag and one rate descriptor. By default k and
#' log(2)/k are not counted as separate evidence. No causal exchange rate is inferred.
#' @param rec Recovery metric data frame.
#' @param forcing_window Positive duration in the same units as recovery time.
#' With a duration, speed is k*T/(1+k*T), and lag score is 1/(1+lag/T).
#' Without one, scores are cohort-relative min-max descriptions.
#' @param lag_weight,rate_weight,halflife_weight Non-negative weights.
#' @param invert_slow TRUE scores faster recovery higher; FALSE reverses all components.
#' @return Input with kinetics_score, component coverage, heuristic class, and lag ratio.
#' @examples
#' rec <- data.frame(tau_lag = c(2, 4, 1), k_recovery = c(0.3, 0.1, 0.5))
#' rri_kinetics_score(rec)$kinetics_score
#' @export
rri_kinetics_score <- function(rec, forcing_window = NULL, lag_weight = 0.30,
  rate_weight = 0.70, halflife_weight = 0, invert_slow = TRUE) {
  rec <- as.data.frame(rec)
  w <- .rri_weights(c(lag_weight, rate_weight, halflife_weight))
  if (!is.null(forcing_window) && (length(forcing_window) != 1L ||
      !is.finite(forcing_window) || forcing_window <= 0)) stop("forcing_window must be positive.")
  if (rate_weight > 0 && halflife_weight > 0)
    warning("k and t_half usually describe the same fitted rate; these weights duplicate evidence.")
  n <- nrow(rec); mat <- matrix(NA_real_, n, 3)
  cols <- c("tau_lag", if ("k_recovery" %in% names(rec)) "k_recovery" else "k", "t_half")
  if (!any(cols %in% names(rec))) stop("No recovery timing columns found.")
  for (j in seq_along(cols)) {
    if (!cols[j] %in% names(rec)) next
    x <- rec[[cols[j]]]
    if (!is.numeric(x)) stop("Recovery metrics must be numeric.")
    x[!is.finite(x) | x < 0] <- NA_real_
    if (is.null(forcing_window)) {
      mat[, j] <- .rri_scale(x)
      if (j != 2) mat[, j] <- 1 - mat[, j]
    } else {
      mat[, j] <- if (j == 2) 1 - 1 / (1 + x * forcing_window) else 1 / (1 + x / forcing_window)
    }
  }
  if (!isTRUE(invert_slow)) mat <- 1 - mat
  z <- .rri_weighted(mat, w)
  rec$kinetics_score <- z$score; rec$kinetics_coverage <- z$coverage
  rec$kinetics_class <- cut(z$score, c(-Inf, 1/3, 2/3, Inf), c("slow", "moderate", "fast"))
  rec$response_to_forcing_ratio <- if (!is.null(forcing_window) && "tau_lag" %in% names(rec))
    rec$tau_lag / forcing_window else rep(NA_real_, n)
  attr(rec, "kinetics_interpretation") <- "Descriptive score; class cutoffs are not empirically calibrated."
  rec
}
