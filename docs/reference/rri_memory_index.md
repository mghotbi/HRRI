# Persistent-displacement and loop-area diagnostic

A descriptive composite, not proof of ecological memory. Slow
relaxation, baseline drift and continuing forcing can also produce
displacement.

## Usage

``` r
rri_memory_index(
  rec,
  H_weight = 0.5,
  I_weight = 0.5,
  lag_weight = 0,
  normalise_inputs = FALSE
)
```

## Arguments

- rec:

  Recovery metric table.

- H_weight, I_weight, lag_weight:

  Non-negative component weights.

- normalise_inputs:

  FALSE uses bounded dimensionless H and I fractions; TRUE requests
  cohort-relative min-max scaling. Lag requires TRUE.

## Value

Input with memory_index, memory_coverage and heuristic memory_class.

## Examples

``` r
rec <- data.frame(H_hysteresis = c(0.1, 0.3, 0.05),
                  incomplete_return_frac = c(0.2, 0.4, 0.1))
rri_memory_index(rec)$memory_index
#> [1] 0.150 0.350 0.075
```
