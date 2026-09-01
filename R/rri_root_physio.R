#' Exploratory root-trait composite
#' @description A weighted standardized trait summary. Trait direction is
#' context-dependent: greater ROL, porosity, aerenchyma or SRL is not universally
#' better plant performance or greater oxygen delivery to every root region.
#' @param plant_df Numeric root-trait measurements.
#' @param biomass_col,length_col Columns for root biomass and length/density.
#' @param rol_col Column for measured radial oxygen loss, with consistent units.
#' @param aerenchyma_col,porosity_col Optional aeration trait columns; avoid
#' double weighting correlated measures of the same anatomical attribute.
#' @param srl_col Column for specific root length.
#' @param w_biomass,w_length,w_rol,w_aerenchyma,w_porosity,w_srl Non-negative weights.
#' @param directions Optional named numeric vector assigning +1 (larger maps to
#' a larger score) or -1 (larger maps to a smaller score) to each selected
#' measurement column. NULL uses +1 for backward compatibility and warns.
#' @param scaling pnorm (normal-CDF scaling, not a calibrated probability) or minmax.
#' @return root_physio_score and trait contributions; unobserved values stay NA.
#' @importFrom stats setNames sd
#' @examples
#' df <- data.frame(ROL = c(0.5, 1.0, 1.5), biomass = c(10, 20, 30))
#' rri_root_physio(df, rol_col = "ROL", biomass_col = "biomass",
#'   directions = c(ROL = 1, biomass = 1))$root_physio_score
#' @export
rri_root_physio <- function(
  plant_df,
  biomass_col = NULL,
  length_col = NULL,
  rol_col = NULL,
  aerenchyma_col = NULL,
  porosity_col = NULL,
  srl_col = NULL,
  w_biomass = 0.20,
  w_length = 0.20,
  w_rol = 0.30,
  w_aerenchyma = 0.20,
  w_porosity = 0.05,
  w_srl = 0.05,
  directions = NULL,
  scaling = c("pnorm", "minmax")
) {
  scaling <- match.arg(scaling)
  plant_df <- as.data.frame(plant_df)

  # -- trait registry ---------------------------------------------------------
  trait_spec <- list(
    biomass    = list(col = biomass_col, w = w_biomass, label = "root_biomass"),
    length     = list(col = length_col, w = w_length, label = "root_length_density"),
    rol        = list(col = rol_col, w = w_rol, label = "ROL"),
    aerenchyma = list(col = aerenchyma_col, w = w_aerenchyma, label = "aerenchyma"),
    porosity   = list(col = porosity_col, w = w_porosity, label = "root_porosity"),
    srl        = list(col = srl_col, w = w_srl, label = "SRL")
  )

  # keep only traits where column is non-NULL and present in data frame
  present <- Filter(
    function(x) !is.null(x$col) && x$col %in% names(plant_df),
    trait_spec
  )

  if (length(present) == 0) {
    stop(
      "None of the specified root trait columns were found in `plant_df`.\n",
      "Provide at least one of: biomass_col, length_col, rol_col, ",
      "aerenchyma_col, porosity_col, srl_col.",
      call. = FALSE
    )
  }
  selected_columns <- vapply(present, `[[`, character(1), "col")
  if (anyDuplicated(selected_columns))
    stop("A measured column cannot be counted as more than one root trait.")
  used_columns <- unique(selected_columns)
  plant_df[, used_columns] <- .rri_numeric_df(
    plant_df[, used_columns, drop = FALSE], "Selected plant traits")

  if (is.null(directions)) {
    warning("No trait directions supplied; using +1 for every selected root trait.",
            call. = FALSE)
    directions <- stats::setNames(rep(1, length(used_columns)), used_columns)
  }
  if (!is.numeric(directions) || is.null(names(directions)) ||
      anyDuplicated(names(directions)) || any(!used_columns %in% names(directions)) ||
      any(!directions[used_columns] %in% c(-1, 1)))
    stop("directions must provide named +1/-1 values for every selected column.")

  # -- scaling helper ---------------------------------------------------------
  scale_col <- function(x) .rri_scale(x, scaling)

  # -- build weighted composite -----------------------------------------------
  n <- nrow(plant_df)
  declared_weights <- vapply(present, `[[`, numeric(1), "w")
  .rri_weights(declared_weights)
  total_w <- sum(declared_weights)

  score_mat <- matrix(NA_real_, nrow = n, ncol = length(present))
  weights <- numeric(length(present))
  labels <- character(length(present))
  cols <- character(length(present))
  mean_vals <- numeric(length(present))

  for (i in seq_along(present)) {
    item <- present[[i]]
    raw <- plant_df[[item$col]]
    sc <- scale_col(raw)
    if (directions[[item$col]] < 0) sc <- 1 - sc
    score_mat[, i] <- sc
    weights[i] <- item$w / total_w
    labels[i] <- item$label
    cols[i] <- item$col
    mean_vals[i] <- .rri_mean(sc)
  }

  combined <- .rri_weighted(score_mat, weights)
  root_physio_score <- combined$score

  root_physio_score <- pmax(0, pmin(1, root_physio_score))

  # -- contributor table ------------------------------------------------------
  contributors <- data.frame(
    column = cols,
    trait_label = labels,
    direction = unname(directions[cols]),
    effective_weight = round(weights, 4),
    mean_scaled_value = round(mean_vals, 4),
    stringsAsFactors = FALSE
  )

  list(
    root_physio_score = root_physio_score,
    contributors      = contributors,
    n_traits_used     = length(present),
    coverage = combined$coverage
  )
}
