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
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
  rri_connectivity_score(res)
} # }
```
