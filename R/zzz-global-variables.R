#' @importFrom utils globalVariables
#' @keywords internal
NULL

utils::globalVariables(
  c(
    # domain scores
    "Physio",
    "Soil",
    "Micro",
    "RRI",
    # dynamic mode scores
    "Physio_dyn",
    "Soil_dyn",
    "Micro_dyn",
    "RRI_dyn",
    # property scores (Capacity--Connectivity--Kinetics--Memory)
    "Capacity",
    "Connectivity",
    "Kinetics",
    "Memory",
    "memory_index",
    "memory_class",
    "capacity_score",
    "connectivity_score",
    # plot helpers
    ".x",
    ".y",
    ".colour",
    ".trajectory",
    ".order_value",
    # benchmark columns (ggplot2 aes bare names)
    "truth",
    "RMSE",
    # recovery metric columns
    "A_norm",
    "O_norm",
    "I_norm",
    "tau_lag",
    "tau_r",
    "t_half",
    "H",
    "trajectory_class",
    "metric_label",
    "value_scaled",
    "value"
  )
)
# Package-wide audit and corrective patches in version 0.99.1 were
# Assisted-by: OpenAI Codex. The maintainers remain responsible for scientific
# interpretation, redistribution rights, runtime verification and maintenance.
