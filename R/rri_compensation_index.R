#' Cross-domain asynchrony diagnostic
#'
#' @description Association or variance cancellation, not evidence of causal
#' buffering. Variance_ratio is 1 - var(rowSums(X))/sum(var(X_j)); positive
#' values indicate cancellation. mean_neg_cor is the negative mean correlation.
#' @param res RRI result.
#' @param per_group Compute by trajectory/group.
#' @param group_cols Group identifiers.
#' @param id Optional aligned identifiers.
#' @param method mean_neg_cor or variance_ratio.
#' @param scale_output For correlations maps `[-1, 1]` to `[0, 1]`; for variance
#' cancellation truncates negative values to zero. A correlation score of 0.5
#' means zero mean correlation, not moderate biological compensation.
#' @return Diagnostic score, correlations and interpretation.
#' @importFrom stats var
#' @importFrom utils combn
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#'   rri_compensation_index(res)
#' }
#' @export
rri_compensation_index <- function(res, per_group = FALSE, group_cols = NULL,
  id = NULL, method = c("mean_neg_cor", "variance_ratio"), scale_output = TRUE) {
  method <- match.arg(method); df <- as.data.frame(res$row_scores)
  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (nrow(id) != nrow(df)) stop("id row count differs.")
    for (nm in intersect(names(id), names(df)))
      if (!identical(as.character(id[[nm]]), as.character(df[[nm]])))
        stop("id and row_scores conflict in column: ", nm)
    df <- cbind(df, id[, setdiff(names(id), names(df)), drop = FALSE])
  }
  calc <- function(d) {
    names_used <- intersect(c("Physio", "Soil", "Micro"), names(d))
    names_used <- names_used[vapply(d[, names_used, drop = FALSE], function(x) any(is.finite(x)), logical(1))]
    if (length(names_used) < 2L) return(list(score = NA_real_, correlations = NULL))
    x <- as.matrix(d[, names_used, drop = FALSE]); x <- x[apply(is.finite(x), 1, all), , drop = FALSE]
    raw <- NA_real_; corrs <- NULL
    if (method == "mean_neg_cor" && nrow(x) >= 4L) {
      cmb <- utils::combn(seq_len(ncol(x)), 2)
      corrs <- apply(cmb, 2, function(p) .rri_cor(x[, p[1]], x[, p[2]], min_n = 4L))
      names(corrs) <- apply(cmb, 2, function(p) paste(names_used[p], collapse="_"))
      raw <- -.rri_mean(corrs)
    }
    if (method == "variance_ratio" && nrow(x) >= 2L) {
      denominator <- sum(apply(x, 2, stats::var))
      if (denominator > 0) raw <- 1 - stats::var(rowSums(x)) / denominator
    }
    score <- if (!scale_output) raw else if (method == "mean_neg_cor") (raw + 1)/2 else max(0, raw)
    list(score = score, correlations = corrs)
  }
  note <- "Descriptive asynchrony only; does not establish functional compensation."
  if (!per_group) {
    z <- calc(df)
    return(list(compensation_index = z$score, pairwise_cor = z$correlations,
      method_used = method, interpretation = note))
  }
  if (!length(group_cols)) stop("Provide group_cols.")
  key <- .rri_key(df, group_cols)
  rows <- lapply(unique(key), function(k) {
    d <- df[key == k, , drop = FALSE]; z <- calc(d)
    cbind(d[1, group_cols, drop = FALSE], compensation_index = z$score, interpretation = note)
  })
  list(compensation_index = do.call(rbind, rows), pairwise_cor = NULL,
    method_used = paste0(method, "_per_group"), interpretation = note)
}
