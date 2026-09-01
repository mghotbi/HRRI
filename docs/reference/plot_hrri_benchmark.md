# Plot descriptive benchmark agreement

Plot descriptive benchmark agreement

## Usage

``` r
plot_hrri_benchmark(bm, print = TRUE, colour = "#0072B2")
```

## Arguments

- bm:

  An hrri_benchmark object.

- print:

  Whether to display the plot.

- colour:

  Scatter colour.

## Value

Invisibly, a two-panel patchwork object. Requires patchwork.

## Examples

``` r
if (FALSE) { # \dontrun{
  b <- benchmark_hrri(domains = "soil", n = 2, missing = 0.1)
  plot_hrri_benchmark(b)
} # }
```
