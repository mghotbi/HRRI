# Draw the reproducible simulation demonstration

Draw the reproducible simulation demonstration

## Usage

``` r
plot_rri_simulation_demo(
  demo,
  figure = c("observations", "coverage", "capacity")
)
```

## Arguments

- demo:

  Result of rri_simulation_demo.

- figure:

  observations, coverage or capacity.

## Value

Invisibly returns the plotted data table; draws on the active device.
Uses base graphics, so no optional plotting package is required.

## Examples

``` r
if (FALSE) { # \dontrun{
  demo <- rri_simulation_demo(seed = 20260830L)
  plot_rri_simulation_demo(demo, figure = "observations")
} # }
```
