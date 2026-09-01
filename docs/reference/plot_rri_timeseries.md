# Aligned time series with separate physical units

Shows Eh, EAC and the observation-derived score in separate panels
sharing time. It does not place unlike units on a common axis.

## Usage

``` r
plot_rri_timeseries(
  sim,
  res,
  plot_id = "P1",
  depth_id = "D1",
  plant_id = "Plant1",
  perturb_start = NULL,
  perturb_end = NULL,
  base_size = 9
)
```

## Arguments

- sim:

  Simulator output containing id and soil_data.

- res:

  Pipeline output containing row_scores.

- plot_id, depth_id, plant_id:

  Identifiers for a single trajectory.

- perturb_start, perturb_end:

  Optional disturbance interval in input time units.

- base_size:

  Base font size.

## Value

A ggplot with three vertically aligned panels.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
                      id = sim$id)
  plot_rri_timeseries(sim, res)
} # }
```
