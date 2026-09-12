# Summarise supported diagnostics without fabricating missing properties

Summarise supported diagnostics without fabricating missing properties

## Usage

``` r
rri_property_scores(
  res,
  rec = NULL,
  soil_df = NULL,
  eac_col = "EAC",
  edc_col = "EDC",
  humic_col = NULL,
  connectivity_method = "cross_domain_magnitude",
  H_weight = 0.5,
  I_weight = 0.5,
  forcing_window = NULL
)
```

## Arguments

- res:

  RRI result.

- rec:

  Optional recovery table.

- soil_df:

  Optional soil capacity measurements.

- eac_col, edc_col, humic_col:

  Capacity-related columns.

- connectivity_method:

  Association or network summary.

- H_weight, I_weight:

  Memory-diagnostic weights.

- forcing_window:

  Optional timescale for the recovery speed score.

## Value

Scores and a provenance table. Unavailable properties stay NA. Capacity
here is an oxidative-oriented feature composite, Connectivity an
association/topology descriptor, Kinetics a recovery-speed descriptor,
and Memory a persistent-displacement descriptor. None proves the named
mechanism.

## Which inputs each property needs

Only Connectivity is derived from `res` alone. The other three require
an additional argument, and are returned as `NA` with method
`"unavailable"` when it is absent:

- **Capacity** needs `soil_df` containing the columns named by `eac_col`
  and `edc_col`. Passing only `res` and `rec` is the usual reason
  Capacity comes back `NA`.

- **Kinetics** and **Memory** need `rec`, the table returned by
  [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md).

A message names any missing input. `NA` here means "not supplied", never
"measured and found to be zero".

## Examples

``` r
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 2, n_time = 40, p_micro = 10,
  seed = 2026
)
res <- suppressWarnings(rri_pipeline_st(
  ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data, id = sim$id
))
rec <- rri_recovery_metrics(
  res = res, id = sim$id, time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = 12, perturb_end = 22
)

## All four properties available: soil_df supplies Capacity, rec supplies
## Kinetics and Memory.
full <- rri_property_scores(res, rec = rec, soil_df = sim$soil_data)
full$property_table
#>       property     score                                         method
#> 1     Capacity 0.5188017 Oxidative-oriented feature composite; not Cacc
#> 2 Connectivity 0.5555565                         cross_domain_magnitude
#> 3     Kinetics 0.5000000                 Cohort-relative recovery speed
#> 4       Memory 0.1657948   Loop-area/persistent-displacement diagnostic
#>   available
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
#> 4      TRUE

## Omitting soil_df leaves Capacity unavailable, and says so.
partial <- rri_property_scores(res, rec = rec)
#> Unavailable property returned as NA.
#>   Capacity: supply `soil_df` with EAC/EDC columns (see `eac_col`, `edc_col`)
partial$property_table
#>       property     score                                       method available
#> 1     Capacity        NA                                  unavailable     FALSE
#> 2 Connectivity 0.5555565                       cross_domain_magnitude      TRUE
#> 3     Kinetics 0.5000000               Cohort-relative recovery speed      TRUE
#> 4       Memory 0.1657948 Loop-area/persistent-displacement diagnostic      TRUE
```
