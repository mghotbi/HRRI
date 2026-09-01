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

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
  rri_property_scores(res)
} # }
```
