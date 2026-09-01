# Sensitivity to domain aggregation weights

Sensitivity to domain aggregation weights

## Usage

``` r
rri_sensitivity(res, weight_grid = seq(0.2, 0.6, by = 0.1))
```

## Arguments

- res:

  RRI result.

- weight_grid:

  Plant weights in (0,1), or a data frame/matrix with named Physio,
  Soil, Micro columns specifying complete alternative weights.

## Value

Alternative normalized weights, finite-pair count and Spearman
correlation. This conditions on the already computed features,
reductions and missingness.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
  rri_sensitivity(res)
} # }
```
