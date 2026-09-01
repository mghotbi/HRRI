# Correlation with a Simulator-Defined Target

Computes a descriptive correlation between per-sample RRI and a
prescribed simulator target. This is an internal simulation benchmark
and is not predictive accuracy, empirical validation or recovery of a
true latent state.

## Usage

``` r
rri_latent_correlation(
  res,
  latent_truth,
  method = c("pearson", "spearman", "kendall")
)
```

## Arguments

- res:

  An object returned by
  [`rri_pipeline_st()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md).

- latent_truth:

  Numeric simulator-target vector. The legacy argument name is retained
  for compatibility.

- method:

  Correlation method. One of `"pearson"`, `"spearman"`, or `"kendall"`.

## Value

A single numeric correlation coefficient.

## Details

This function is designed for simulation benchmarking. In empirical
datasets, no known latent state exists and this metric should not be
used.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline_st(sim$ROS_flux, sim$Eh_stability)
  rri_latent_correlation(res, sim$latent_truth)
} # }
```
