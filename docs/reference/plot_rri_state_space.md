# Plot domain-score space with correctly matched trajectory diagnostics

Plot domain-score space with correctly matched trajectory diagnostics

## Usage

``` r
plot_rri_state_space(
  res,
  rec = NULL,
  x_property = c("Physio", "Connectivity", "Soil", "Micro"),
  y_property = c("Soil", "Micro", "Physio", "Kinetics"),
  colour_by = c("RRI", "Memory", "trajectory_class"),
  group_cols = c("plot", "depth", "plant_id"),
  base_size = 12
)
```

## Arguments

- res:

  RRI result with aligned identifiers in row_scores.

- rec:

  Optional one-row-per-trajectory diagnostic table.

- x_property, y_property:

  Domain scores, or group association/kinetics.

- colour_by:

  RRI, Memory, or an explicitly supplied trajectory_class.

- group_cols:

  Full key shared between row_scores and rec.

- base_size:

  Plot font size.

## Value

ggplot. Domain scores are not relabelled as mechanistic properties.

## Examples

``` r
# \donttest{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
                      id = sim$id)
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
  plot_rri_state_space(res, group_cols = c("plot", "depth", "plant_id"))

# }
```
