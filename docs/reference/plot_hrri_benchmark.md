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
# \donttest{
  b <- benchmark_hrri(domains = "soil", n = 2, missing = 0.1)
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> benchmark_hrri: 1/2 seeds attempted
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> benchmark_hrri: 2/2 seeds attempted
  plot_hrri_benchmark(b)

# }
```
