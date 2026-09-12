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
# \donttest{
  demo <- rri_simulation_demo(seed = 20260830L)
  head(demo$scores)
#>                  scenario panel      unit_id time       RRI domain_coverage
#> flood_drain.1 flood_drain  full P1.D1.Plant1    1 0.8976458               1
#> flood_drain.2 flood_drain  full P2.D1.Plant1    1 0.8385890               1
#> flood_drain.3 flood_drain  full P1.D2.Plant1    1 0.8688453               1
#> flood_drain.4 flood_drain  full P2.D2.Plant1    1 0.8344818               1
#> flood_drain.5 flood_drain  full P1.D1.Plant2    1 0.7993227               1
#> flood_drain.6 flood_drain  full P2.D1.Plant2    1 0.8950797               1
#>               n_domains
#> flood_drain.1         3
#> flood_drain.2         3
#> flood_drain.3         3
#> flood_drain.4         3
#> flood_drain.5         3
#> flood_drain.6         3
# }
```
