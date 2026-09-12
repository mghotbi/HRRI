# Cross-domain association or graph-topology summary

Retains the legacy function name. Correlation magnitude is association,
not electron-transfer encounter probability or measured alpha. Network
summaries concern unweighted topology, not biochemical connectivity.

## Usage

``` r
rri_connectivity_score(
  res,
  method = c("cross_domain_magnitude", "network"),
  per_group = FALSE,
  group_cols = NULL
)
```

## Arguments

- res:

  RRI result; graph method uses meta\$graph.

- method:

  cross_domain_magnitude or network.

- per_group:

  Compute association by group.

- group_cols:

  Required grouping columns when per_group=TRUE.

## Value

Score, association coefficients and method/provenance information.

## Examples

``` r
# \donttest{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
  rri_connectivity_score(res)
#> $connectivity_score
#> [1] 0.4272939
#> 
#> $method_used
#> [1] "cross_domain_magnitude"
#> 
#> $pairwise_abs_cor
#>  Physio_Soil Physio_Micro   Soil_Micro 
#>    0.4272939           NA           NA 
#> 
#> $n_pairs_used
#> [1] 1
#> 
#> $network_metrics
#> NULL
#> 
#> $interpretation
#> [1] "Descriptive association; shared forcing and repeated measures can explain correlation"
#> 
# }
```
