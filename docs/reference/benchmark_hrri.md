# Benchmark diagnostic agreement with a simulator-defined target

Runs independent seeds and reports descriptive agreement, not held-out
prediction, empirical validation or parameter identification.

## Usage

``` r
benchmark_hrri(
  domains = c("soil", "plant", "micro"),
  missing = 0,
  noise = 0,
  n = 50L,
  seed_start = 1L,
  sim_args = NULL,
  verbose = TRUE,
  pipeline_args = list()
)

# S3 method for class 'hrri_benchmark'
print(x, ...)
```

## Arguments

- domains:

  Nonempty subset of soil, plant and micro.

- missing:

  Fraction in `[0, 1)` of uniformly sampled cells removed (MCAR). This
  does not implement informative or MNAR missingness.

- noise:

  Nonnegative Gaussian noise SD on each column's original scale. A
  common SD has different relative effects on differently scaled
  variables.

- n:

  Positive integer number of independent simulations.

- seed_start:

  First integer seed. The caller's RNG state is restored.

- sim_args:

  Named simulator arguments, excluding seed.

- verbose:

  Print progress messages.

- pipeline_args:

  Named additional rri_pipeline arguments, excluding dat, soil, plant,
  micro and id; use this to justify orientations and weights.

- x:

  An hrri_benchmark object.

- ...:

  Unused method arguments.

## Value

An hrri_benchmark list with summary, seed_metrics, row_data, failures
(seed and error), and settings including package/R versions. Summary
RMSE and Bias compare the score directly with the chosen target. r_truth
and rank_truth are pooled descriptive correlations; spread_association
compares within-seed score and target SDs. None is interval coverage or
predictive uncertainty. n_rows counts finite matched rows across all
seeds.

## Details

Soil and plant synthetic observations and log1p gene abundances feed the
pipeline. Latent architecture and microbial activity states are not
scoring inputs. The latent target is still a prescribed simulator
composite; it is not an independently measured recovery outcome.
Replicated rows within seeds are dependent. Compare methods using
held-out seeds and independent process outcomes in a separate validation
design.

## See also

rri_pipeline, simulate_redox_holobiont, plot_hrri_benchmark

## Examples

``` r
if (FALSE) { # \dontrun{
b <- benchmark_hrri(domains = "soil", n = 2, missing = 0.1)
print(b)
} # }
```
