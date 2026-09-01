#' Score observed soil, plant and microbial panels
#'
#' @description Exploratory integration of available numeric domain observations.
#' A larger score is not automatically greater resilience: justify feature
#' orientation, the reference function and the observation window.
#' @param dat Optional wide data frame with canonical observation names.
#' @param soil,plant,micro Optional numeric data frames with aligned rows.
#' Supply these instead of dat to use custom measurement names or partial panels.
#' @param id Optional identifier data frame in the same row order.
#' @param domain_weights Named nonnegative weights for Physio, Soil and Micro.
#' Available positive weights are renormalized per row; absent domains stay NA.
#' @param ... Arguments to rri_pipeline_st, excluding its domain inputs,
#' identifiers and w1/w2/w3 (use domain_weights instead).
#' @details Known hidden simulator columns are excluded. This is a safeguard,
#' not an automatic detector of every possible source of target leakage.
#' Cohort-fitted PCA and scaling must not be interpreted as a trained predictor.
#' Use rri_reference_scores for fixed, independently justified reference anchors.
#' A reduced panel changes the estimand; compare panels through sensitivity
#' analysis rather than treating their scores as interchangeable.
#' @return An RRI object with row_scores, a scores alias, effective_weights,
#' per-row domain_coverage and n_domains, and a call_mode field.
#' @seealso rri_pipeline_st, rri_reference_scores, benchmark_hrri
#' @examples
#' x <- data.frame(Eh = c(50, 100, 150, 200), pH = c(5, 5.5, 6, 6.5))
#' z <- rri_pipeline(soil = x, method_soil = "scale",
#'                   direction_anchor_soil = "Eh")
#' z$row_scores
#' @importFrom utils modifyList
#' @export
rri_pipeline <- function(dat = NULL, soil = NULL, plant = NULL, micro = NULL,
                         id = NULL,
                         domain_weights = c(Physio = 0.4, Soil = 0.35, Micro = 0.25),
                         ...) {
  # Assisted-by: OpenAI Codex. Explicit input routing and weight forwarding.
  domains <- c("Physio", "Soil", "Micro")
  if (is.null(names(domain_weights)) || anyDuplicated(names(domain_weights)) ||
      !setequal(names(domain_weights), domains))
    stop("domain_weights must have exactly the names Physio, Soil and Micro.")
  w <- .rri_weights(domain_weights[domains])
  dots <- list(...)
  blocked <- intersect(names(dots), c("ROS_flux", "Eh_stability", "micro_data",
                                      "id", "w1", "w2", "w3"))
  if (length(blocked)) stop("Use the named wrapper arguments, not: ",
                            paste(blocked, collapse = ", "))
  call_mode <- "raw"
  if (!is.null(dat)) {
    if (!is.data.frame(dat)) stop("dat must be a wide data frame, not a simulator list.")
    if (!is.null(soil) || !is.null(plant) || !is.null(micro))
      stop("Supply either dat or separate domain blocks, not both.")
    if (anyDuplicated(names(dat))) stop("dat has duplicate column names.")
    soil_names <- c("Eh", "pH", "FeII_mmol_kg", "FeIII_mmol_kg",
      "FeII_aq_mmol_kg", "FeIII_reactive_mmol_kg", "FeIII_crystalline_mmol_kg",
      "Fe2_aq", "Fe2.Fe3", "EAC", "EDC", "DOC", "water_content",
      "air_filled_porosity", "pore_connectivity", "aqueous_connectivity",
      "oxygen_availability", "MnII_mmol_kg", "MnIII_mmol_kg", "MnIV_mmol_kg",
      "NO3_mmol_kg", "NH4_mmol_kg", "SO4_mmol_kg", "sulfide_mmol_kg")
    plant_names <- c("SPAD", "FvFm", "PhiPSII", "NPQ", "ROL", "root_exudates",
      "organic_acids", "phenolics", "exudate_redox_activity", "aerenchyma",
      "ROL_barrier", "root_oxidative_stress", "root_redox_buffering", "Fe_plaque_proxy")
    micro_names <- c("CH4", "Denitrification", "EET", "aerobic_respiration",
      "denitrification", "Fe_Mn_reduction", "EET_potential", "EET_reduction",
      "sulfate_reduction", "methanogenesis", "microbial_redox_flexibility")
    pick <- function(nms) {
      keep <- intersect(nms, names(dat))
      if (length(keep)) dat[, keep, drop = FALSE] else NULL
    }
    plant <- pick(plant_names); soil <- pick(soil_names); micro <- pick(micro_names)
    if (is.null(id)) {
      keys <- intersect(c("seed", "plot", "depth", "plant_id", "time", "row_id",
                          "unit_id", "history_pair", "history", "scenario",
                          "rescue", "cycle", "phase"), names(dat))
      if (length(keys)) id <- dat[, keys, drop = FALSE]
    }
    call_mode <- "data.frame"
  }
  if (is.null(soil) && is.null(plant) && is.null(micro))
    stop("No observation block found. For custom names supply soil, plant or micro explicitly.")
  res <- do.call(rri_pipeline_st, c(list(ROS_flux = plant, Eh_stability = soil,
    micro_data = micro, id = id, w1 = unname(w[1L]), w2 = unname(w[2L]),
    w3 = unname(w[3L])), dots))
  res$scores <- res$row_scores
  res$call_mode <- call_mode
  res
}
