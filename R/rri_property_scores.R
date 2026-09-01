#' Summarise supported diagnostics without fabricating missing properties
#'
#' @param res RRI result.
#' @param rec Optional recovery table.
#' @param soil_df Optional soil capacity measurements.
#' @param eac_col,edc_col,humic_col Capacity-related columns.
#' @param connectivity_method Association or network summary.
#' @param H_weight,I_weight Memory-diagnostic weights.
#' @param forcing_window Optional timescale for the recovery speed score.
#' @return Scores and a provenance table. Unavailable properties stay NA.
#' Capacity here is an oxidative-oriented feature composite, Connectivity an
#' association/topology descriptor, Kinetics a recovery-speed descriptor, and
#' Memory a persistent-displacement descriptor. None proves the named mechanism.
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#'   rri_property_scores(res)
#' }
#' @export
rri_property_scores <- function(res, rec = NULL, soil_df = NULL,
  eac_col = "EAC", edc_col = "EDC", humic_col = NULL,
  connectivity_method = "cross_domain_magnitude", H_weight = 0.50,
  I_weight = 0.50, forcing_window = NULL) {
  scores <- stats::setNames(rep(NA_real_, 4), c("Capacity", "Connectivity", "Kinetics", "Memory"))
  methods <- stats::setNames(rep("unavailable", 4), names(scores))
  if (!is.null(soil_df)) {
    cap <- rri_capacity_index(soil_df, eac_col = eac_col, edc_col = edc_col, humic_col = humic_col)
    scores["Capacity"] <- .rri_mean(cap$capacity_score)
    methods["Capacity"] <- "Oxidative-oriented feature composite; not Cacc"
  }
  conn <- rri_connectivity_score(res, method = connectivity_method)
  scores["Connectivity"] <- .rri_mean(conn$connectivity_score)
  methods["Connectivity"] <- conn$method_used
  if (!is.null(rec)) {
    k <- rri_kinetics_score(rec, forcing_window = forcing_window)
    m <- rri_memory_index(rec, H_weight = H_weight, I_weight = I_weight)
    scores["Kinetics"] <- .rri_mean(k$kinetics_score)
    scores["Memory"] <- .rri_mean(m$memory_index)
    methods["Kinetics"] <- if (is.null(forcing_window)) "Cohort-relative recovery speed" else "Timescale-relative recovery speed"
    methods["Memory"] <- "Loop-area/persistent-displacement diagnostic"
  }
  list(property_scores = scores,
    property_table = data.frame(property = names(scores), score = unname(scores),
      method = unname(methods), available = is.finite(scores), row.names = NULL),
    rri_summary = .rri_mean(res$row_scores$RRI))
}
