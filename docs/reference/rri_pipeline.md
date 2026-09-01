# Score observed soil, plant and microbial panels

Exploratory integration of available numeric domain observations. A
larger score is not automatically greater resilience: justify feature
orientation, the reference function and the observation window.

## Usage

``` r
rri_pipeline(
  dat = NULL,
  soil = NULL,
  plant = NULL,
  micro = NULL,
  id = NULL,
  domain_weights = c(Physio = 0.4, Soil = 0.35, Micro = 0.25),
  ...
)
```

## Arguments

- dat:

  Optional wide data frame with canonical observation names.

- soil, plant, micro:

  Optional numeric data frames with aligned rows. Supply these instead
  of dat to use custom measurement names or partial panels.

- id:

  Optional identifier data frame in the same row order.

- domain_weights:

  Named nonnegative weights for Physio, Soil and Micro. Available
  positive weights are renormalized per row; absent domains stay NA.

- ...:

  Arguments to rri_pipeline_st, excluding its domain inputs, identifiers
  and w1/w2/w3 (use domain_weights instead).

## Value

An RRI object with row_scores, a scores alias, effective_weights,
per-row domain_coverage and n_domains, and a call_mode field.

## Details

Known hidden simulator columns are excluded. This is a safeguard, not an
automatic detector of every possible source of target leakage.
Cohort-fitted PCA and scaling must not be interpreted as a trained
predictor. Use rri_reference_scores for fixed, independently justified
reference anchors. A reduced panel changes the estimand; compare panels
through sensitivity analysis rather than treating their scores as
interchangeable.

## See also

rri_pipeline_st, rri_reference_scores, benchmark_hrri

## Examples

``` r
x <- data.frame(Eh = c(50, 100, 150, 200), pH = c(5, 5.5, 6, 6.5))
z <- rri_pipeline(soil = x, method_soil = "scale",
                  direction_anchor_soil = "Eh")
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
z$row_scores
#>   Physio      Soil Micro       RRI domain_coverage n_domains Micro_abundance
#> 1     NA 0.0000000    NA 0.0000000            0.35         1              NA
#> 2     NA 0.3333333    NA 0.3333333            0.35         1              NA
#> 3     NA 0.6666667    NA 0.6666667            0.35         1              NA
#> 4     NA 1.0000000    NA 1.0000000            0.35         1              NA
#>   Micro_network Micro_mfa
#> 1            NA        NA
#> 2            NA        NA
#> 3            NA        NA
#> 4            NA        NA
```
