#' Oxidative-oriented soil feature composite
#'
#' @description Standardizes selected measured features and averages them with
#' declared weights. EAC is positive and EDC is inverted after scaling. This is
#' an oxidative-oriented descriptor, not accessible capacity, a redox potential,
#' or a universal ranking of resilience. High EDC means greater reducing capacity.
#' @param soil_df Numeric soil measurements.
#' @param eac_col,edc_col Columns of electron-accepting/donating capacity.
#' @param reactive_fe_col,poorly_cryst_fe_col Optional Fe-pool proxy columns.
#' @param humic_col Optional organic redox proxy column.
#' @param w_eac,w_edc,w_fe,w_humic Non-negative weights.
#' @param scaling pnorm, minmax, or reference; pnorm is not a calibrated probability.
#' @param ref_ranges Named increasing finite ranges in original measurement units;
#' mandatory for every selected column when scaling=reference.
#' @return capacity_score, EAC/(EAC+EDC) ratio, contributors and observed coverage.
#' Coverage concerns selected available columns; overlapping pools must not be
#' mistaken for independent evidence. Common assays and reference ranges are
#' necessary but insufficient for cross-study comparability.
#' @examples
#' df <- data.frame(EAC = c(10, 20, 30), EDC = c(5, 8, 12))
#' rri_capacity_index(df)$capacity_score
#' @importFrom stats pnorm
#' @export
rri_capacity_index <- function(
  soil_df,
  eac_col = "EAC",
  edc_col = "EDC",
  reactive_fe_col = NULL,
  poorly_cryst_fe_col = NULL,
  humic_col = NULL,
  w_eac = 0.35,
  w_edc = 0.25,
  w_fe = 0.25,
  w_humic = 0.15,
  scaling = c("pnorm", "minmax", "reference"),
  ref_ranges = NULL
) {
  scaling <- match.arg(scaling)
  soil_df <- as.data.frame(soil_df)

  # ---- resolve which columns are present ----------------------------------------
  col_spec <- list(
    eac   = list(col = eac_col, w = w_eac, label = "EAC"),
    edc   = list(col = edc_col, w = w_edc, label = "EDC"),
    fe    = list(col = reactive_fe_col, w = w_fe, label = "reactive_Fe"),
    humic = list(col = humic_col, w = w_humic, label = "humic_redox"),
    pfe   = list(col = poorly_cryst_fe_col, w = w_fe * 0.6, label = "poorly_cryst_Fe")
  )

  present <- Filter(
    function(x) !is.null(x$col) && x$col %in% names(soil_df),
    col_spec
  )

  if (length(present) == 0) {
    stop(
      "None of the specified columns were found in `soil_df`. ",
      "Check column names: eac_col, edc_col, reactive_fe_col, humic_col.",
      call. = FALSE
    )
  }
  selected_columns <- vapply(present, `[[`, character(1), "col")
  if (anyDuplicated(selected_columns))
    stop("A measured column cannot be counted as more than one capacity component.")
  used_columns <- unique(selected_columns)
  soil_df[, used_columns] <- .rri_numeric_df(
    soil_df[, used_columns, drop = FALSE], "Selected soil features")

  has_eac <- !is.null(eac_col) && eac_col %in% names(soil_df)
  has_edc <- !is.null(edc_col) && edc_col %in% names(soil_df)

  if (!has_eac && !has_edc) {
    warning(
      "Neither EAC nor EDC column found. Capacity index will reflect ",
      "Fe and humic pools only.",
      call. = FALSE
    )
  }

  # ---- scale helper --------------------------------------------------------------
  scale_col <- function(x, col_name = NULL) {
    .rri_scale(x, scaling, if (scaling == "reference") ref_ranges[[col_name]] else NULL)
  }

  # ---- build weighted composite --------------------------------------------------
  declared_weights <- vapply(present, `[[`, numeric(1), "w")
  .rri_weights(declared_weights)
  total_w <- sum(declared_weights)
  scores <- matrix(NA_real_, nrow = nrow(soil_df), ncol = length(present))
  weights <- numeric(length(present))
  labels <- character(length(present))

  for (i in seq_along(present)) {
    item <- present[[i]]
    raw <- soil_df[[item$col]]

    # Scale in the measured units first; then orient EDC toward oxidation.
    # EDC is reducing capacity, not a universally adverse resilience attribute.
    scores[, i] <- scale_col(raw, col_name = item$col)
    if (identical(item$label, "EDC")) scores[, i] <- 1 - scores[, i]
    weights[i] <- item$w / total_w # re-normalise
    labels[i] <- item$label
  }

  combined <- .rri_weighted(scores, weights)
  capacity_score <- combined$score
  capacity_score <- pmax(0, pmin(1, capacity_score))

  # ---- EAC/(EAC+EDC) oxidative buffering ratio -----------------------------------
  eac_edc_ratio <- rep(NA_real_, nrow(soil_df))

  if (has_eac && has_edc) {
    eac_raw <- soil_df[[eac_col]]
    edc_raw <- soil_df[[edc_col]]
    eac_raw[!is.finite(eac_raw) | eac_raw < 0] <- NA_real_
    edc_raw[!is.finite(edc_raw) | edc_raw < 0] <- NA_real_
    denom <- eac_raw + edc_raw
    denom[denom <= 0] <- NA_real_
    eac_edc_ratio <- eac_raw / denom
  }

  # ---- contributor table ---------------------------------------------------------
  contributors <- data.frame(
    column = vapply(present, `[[`, character(1), "col"),
    property_label = labels,
    effective_weight = weights,
    stringsAsFactors = FALSE
  )

  list(
    capacity_score = capacity_score,
    eac_edc_ratio  = eac_edc_ratio,
    contributors   = contributors,
    coverage = combined$coverage,
    interpretation = "Oxidative-oriented feature composite; not accessible capacity or universal resilience"
  )
}
