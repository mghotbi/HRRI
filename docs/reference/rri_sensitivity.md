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
# \donttest{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
  rri_sensitivity(res)
#>         weight_physio weight_soil weight_micro n_pairs
#> Physio            0.2        0.40         0.40    1440
#> Physio1           0.3        0.35         0.35    1440
#> Physio2           0.4        0.30         0.30    1440
#> Physio3           0.5        0.25         0.25    1440
#> Physio4           0.6        0.20         0.20    1440
#>         spearman_rank_correlation
#> Physio                  0.9812053
#> Physio1                 0.9966331
#> Physio2                 0.9987503
#> Physio3                 0.9844539
#> Physio4                 0.9571665
# }
```
