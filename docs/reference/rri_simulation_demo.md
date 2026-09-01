# Reproducible observable-only HRRI demonstration

Generates both forcing scenarios with the package simulator, calibrates
illustrative targets using a separate baseline simulation, and compares
full, 4-soil/3-plant/2-microbe, and single-domain panels with unchanged
feature targets and tolerances. No hidden capacity, alpha, k, memory or
latent_truth column enters the observed-feature index. This is a
software demonstration, not a validation of latent-state recovery or
ecological prediction.

## Usage

``` r
rri_simulation_demo(seed = 20260830L)
```

## Arguments

- seed:

  Non-negative integer seed.

## Value

Simulation objects, explicit reference specifications and all plotted
tables. The capacity-horizon table is an internal equation check using
known synthetic parameters; it is deliberately separate from observable
scoring.

## Examples

``` r
if (FALSE) { # \dontrun{
  demo <- rri_simulation_demo(seed = 20260830L)
  head(demo$scores)
} # }
```
