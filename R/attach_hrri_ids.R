#' Attach design identifiers to a score table with explicit alignment checks
#'
#' @description Joins experimental design identifiers onto a pipeline score
#' table. Alignment is established by a shared unique \code{row_id}, or by a
#' complete shared observation key, or - only as a last resort and with an
#' explicit warning - by preserved input row order. Matching row counts alone
#' do not establish alignment, so the method actually used is recorded in the
#' \code{id_alignment} attribute of the returned data frame.
#'
#' @param scores Data frame of row-level scores, typically
#'   \code{res$row_scores} from \code{\link{rri_pipeline}} or
#'   \code{\link{rri_pipeline_st}}.
#' @param id Data frame of design identifiers, typically \code{sim$id} from
#'   \code{\link{simulate_redox_holobiont}}.
#' @param key Optional character vector naming the observation key columns to
#'   join on. \code{NULL} (default) selects a key automatically: \code{row_id}
#'   when present and unique in both inputs, otherwise the intersection of
#'   \code{c("plot", "depth", "plant_id", "time")} present in both inputs.
#'
#' @return \code{scores} with the non-conflicting columns of \code{id} attached.
#' The \code{id_alignment} attribute is a list giving the \code{method}
#' (\code{"row_id"}, \code{"observation_key"} or \code{"row_order"}), the
#' \code{key} columns used, and the number of rows matched. Shared identifier
#' columns already present in \code{scores} are checked for conflicts rather
#' than silently overwritten.
#'
#' @seealso \code{\link{rri_pipeline}}, \code{\link{rri_pipeline_st}}
#'
#' @examples
#' scores <- data.frame(row_id = 1:4, RRI = c(0.4, 0.6, 0.5, 0.7))
#' ids <- data.frame(row_id = 1:4,
#'                   plot = c("P1", "P1", "P2", "P2"),
#'                   time = c(1, 2, 1, 2))
#' out <- attach_hrri_ids(scores, ids)
#' attr(out, "id_alignment")$method
#' head(out)
#'
#' @export
attach_hrri_ids <- function(scores, id, key = NULL) {
  scores <- as.data.frame(scores)
  id <- as.data.frame(id)

  if (!nrow(scores)) stop("`scores` has no rows.")
  if (!nrow(id)) stop("`id` has no rows.")

  # ---- choose an alignment key ------------------------------------------
  auto_key <- function() {
    if ("row_id" %in% names(scores) && "row_id" %in% names(id) &&
        !anyDuplicated(scores$row_id) && !anyDuplicated(id$row_id)) {
      return(list(key = "row_id", method = "row_id"))
    }
    cand <- intersect(c("plot", "depth", "plant_id", "time"), names(scores))
    cand <- intersect(cand, names(id))
    if (length(cand) &&
        !anyDuplicated(scores[, cand, drop = FALSE]) &&
        !anyDuplicated(id[, cand, drop = FALSE])) {
      return(list(key = cand, method = "observation_key"))
    }
    list(key = character(0), method = "row_order")
  }

  if (is.null(key)) {
    sel <- auto_key()
  } else {
    key <- as.character(key)
    if (!all(key %in% names(scores)) || !all(key %in% names(id)))
      stop("All `key` columns must be present in both `scores` and `id`.")
    if (anyDuplicated(scores[, key, drop = FALSE]) ||
        anyDuplicated(id[, key, drop = FALSE]))
      stop("`key` must be unique within both `scores` and `id`.")
    sel <- list(key = key,
                method = if (identical(key, "row_id")) "row_id" else "observation_key")
  }

  # ---- conflict check on shared non-key columns -------------------------
  shared <- setdiff(intersect(names(scores), names(id)), sel$key)
  add_cols <- setdiff(names(id), c(sel$key, shared))

  # ---- join --------------------------------------------------------------
  if (sel$method == "row_order") {
    if (nrow(scores) != nrow(id))
      stop("Cannot align: no shared unique key and row counts differ (",
           nrow(scores), " vs ", nrow(id), ").")
    warning("No shared unique key found; aligning by input row order. ",
            "This is valid only if the pipeline preserved row order.",
            call. = FALSE)
    out <- cbind(id[, c(sel$key, add_cols), drop = FALSE], scores)
    n_matched <- nrow(scores)
  } else {
    idx <- match(
      do.call(paste, c(scores[, sel$key, drop = FALSE], sep = "\r")),
      do.call(paste, c(id[, sel$key, drop = FALSE], sep = "\r"))
    )
    n_matched <- sum(!is.na(idx))
    if (n_matched == 0L)
      stop("No rows of `scores` matched `id` on key: ",
           paste(sel$key, collapse = ", "))
    if (n_matched < nrow(scores))
      warning(nrow(scores) - n_matched, " of ", nrow(scores),
              " score rows did not match `id`; identifiers are NA there.",
              call. = FALSE)

    # verify shared columns agree where both are observed
    for (nm in shared) {
      a <- scores[[nm]]
      b <- id[[nm]][idx]
      both <- !is.na(a) & !is.na(b)
      if (any(both) && !isTRUE(all.equal(a[both], b[both], check.attributes = FALSE)))
        stop("`scores` and `id` disagree in shared column: ", nm)
    }

    out <- if (length(add_cols))
      cbind(id[idx, add_cols, drop = FALSE], scores) else scores
  }

  rownames(out) <- NULL
  attr(out, "id_alignment") <- list(
    method      = sel$method,
    key         = sel$key,
    n_rows      = nrow(scores),
    n_matched   = n_matched,
    cols_added  = add_cols,
    cols_shared = shared
  )
  out
}
