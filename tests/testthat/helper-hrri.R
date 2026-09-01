## tests/testthat/helper-hrri.R
## Shared helpers — automatically sourced by testthat before any test file.

#' Small reproducible simulation (2 plots x 2 depths x 2 plants x 10 time)
small_sim <- function() {
  simulate_redox_holobiont(
    n_plot = 2, n_depth = 2, n_plant = 2, n_time = 10,
    p_micro = 5, seed = 42,
    scenario = "flood_drain", n_cycles = 1
  )
}

#' Reservoir spec matching simulate_redox_holobiont() column names
sim_reservoirs <- function() {
  list(
    reactive_FeIII = list(
      Q_col = "FeIII_poor_crystalline_mmol_kg",
      alpha = "alpha_accept", k = "k_accept_h", type = "EAC"
    ),
    crystalline_FeIII = list(
      Q_col = "FeIII_crystalline_mmol_kg",
      alpha = 0.20, k = 0.008, type = "EAC"
    ),
    FeII_pool = list(
      Q_col = "FeII_mmol_kg",
      alpha = "alpha_donate", k = "k_donate_h", type = "EDC"
    )
  )
}
