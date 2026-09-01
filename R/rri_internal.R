# Internal utilities shared by the reviewed modules.
.rri_scale <- function(x, method = "minmax", range = NULL) {
  if (!is.numeric(x)) stop("Score inputs must be numeric.", call. = FALSE)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  if (method == "reference") {
    if (!is.numeric(range) || length(range) != 2L ||
        any(!is.finite(range)) || range[2] <= range[1]) {
      stop("Each reference range must contain two increasing finite values.", call. = FALSE)
    }
    out[ok] <- (x[ok] - range[1]) / diff(range)
  } else if (method == "pnorm") {
    s <- stats::sd(x[ok])
    out[ok] <- if (!is.finite(s) || s == 0) 0.5 else
      stats::pnorm((x[ok] - mean(x[ok])) / s)
  } else {
    rr <- base::range(x[ok])
    out[ok] <- if (diff(rr) == 0) 0.5 else (x[ok] - rr[1]) / diff(rr)
  }
  pmax(0, pmin(1, out))
}

.rri_weights <- function(w) {
  if (!is.numeric(w) || !length(w) || any(!is.finite(w)) ||
      any(w < 0) || sum(w) <= 0) {
    stop("Weights must be finite, non-negative, and sum to > 0.", call. = FALSE)
  }
  w / sum(w)
}

.rri_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

.rri_weighted <- function(x, w) {
  x <- as.matrix(x)
  if (!is.numeric(x) || ncol(x) != length(w))
    stop("Numeric score columns and weights must have matching lengths.", call. = FALSE)
  w <- .rri_weights(w)
  ok <- is.finite(x)
  x[!ok] <- 0
  den <- as.numeric(ok %*% w)
  out <- as.numeric(x %*% w)
  out[den > 0] <- out[den > 0] / den[den > 0]
  out[den == 0] <- NA_real_
  list(score = out, coverage = den, n_observed = rowSums(ok & rep(w > 0, each = nrow(x))))
}

.rri_key <- function(df, cols) {
  if (!length(cols)) return(rep("all", nrow(df)))
  if (any(!cols %in% names(df))) stop("Grouping columns are missing.", call. = FALSE)
  if (anyNA(df[, cols, drop = FALSE])) stop("Grouping identifiers cannot be missing.", call. = FALSE)
  # Length-prefix encoding avoids collisions when identifiers contain separators.
  parts <- lapply(df[, cols, drop = FALSE], function(z) {
    z <- as.character(z); paste0(nchar(z, type = "bytes"), ":", z)
  })
  do.call(paste0, parts)
}

.rri_cor <- function(x, y, method = "pearson", min_n = 3L) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y))
    stop("Correlation inputs must be equally sized numeric vectors.", call. = FALSE)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_n || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  stats::cor(x[ok], y[ok], method = method)
}

# Assisted-by: OpenAI Codex. Explicit schema and missingness checks.
.rri_numeric_df <- function(x, label = "data") {
  x <- as.data.frame(x)
  if (anyDuplicated(names(x))) stop(label, " has duplicate column names.", call. = FALSE)
  bad <- names(x)[!vapply(x, is.numeric, logical(1))]
  if (length(bad)) stop(label, " must be numeric; invalid columns: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  x[] <- lapply(x, function(v) { v[!is.finite(v)] <- NA_real_; v })
  x
}

.rri_observed_block <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x)
  hidden <- c("alpha_accept", "alpha_donate", "k_accept", "k_donate",
    "k_accept_h", "k_donate_h", "q_accept", "q_donate", "cacc_raw",
    "cacc_eac", "cacc_edc", "cacc_total", "cacc_fraction",
    "net_oxidative_balance", "latent_truth", "latent_physio", "latent_soil",
    "latent_micro", "memory", "memory_state", "microbial_memory",
    "capacity", "connectivity", "kinetics", "truth", "crystallinity",
    "redox_position")
  keep <- !tolower(names(x)) %in% hidden
  if (any(!keep)) warning("Excluding simulator-derived hidden columns from scoring: ",
                          paste(names(x)[!keep], collapse = ", "), call. = FALSE)
  x <- x[, keep, drop = FALSE]
  if (!ncol(x)) NULL else .rri_numeric_df(x, "Domain data")
}

.rri_align_scores <- function(scores, id, keys) {
  scores <- as.data.frame(scores); id <- as.data.frame(id)
  if (anyDuplicated(names(scores)) || anyDuplicated(names(id)))
    stop("Identifier and score columns must be uniquely named.")
  if (nrow(scores) != nrow(id)) stop("Identifiers and scores have different row counts.")
  if (all(keys %in% names(scores)) && length(keys)) {
    a <- .rri_key(id, keys); b <- .rri_key(scores, keys)
    if (anyDuplicated(a) || anyDuplicated(b)) stop("Observation keys must be unique.")
    idx <- match(a, b)
  } else if ("row_id" %in% names(id) && "row_id" %in% names(scores)) {
    a <- .rri_key(id, "row_id"); b <- .rri_key(scores, "row_id")
    if (anyDuplicated(a) || anyDuplicated(b)) stop("row_id must be unique.")
    idx <- match(a, b)
  } else {
    # Scores-only output requires the documented input-order contract.
    idx <- seq_len(nrow(id))
  }
  if (anyNA(idx)) stop("Some identifiers have no matching score row.")
  scores <- scores[idx, , drop = FALSE]
  for (nm in intersect(names(id), names(scores))) {
    if (!identical(as.character(id[[nm]]), as.character(scores[[nm]])))
      stop("Conflicting score identifiers: ", nm)
  }
  out <- cbind(id, scores[, setdiff(names(scores), names(id)), drop = FALSE])
  rownames(out) <- NULL
  out
}
