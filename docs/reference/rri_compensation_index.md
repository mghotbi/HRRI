# Cross-domain asynchrony diagnostic

Association or variance cancellation, not evidence of causal buffering.
Variance_ratio is 1 - var(rowSums(X))/sum(var(X_j)); positive values
indicate cancellation. mean_neg_cor is the negative mean correlation.

## Usage

``` r
rri_compensation_index(
  res,
  per_group = FALSE,
  group_cols = NULL,
  id = NULL,
  method = c("mean_neg_cor", "variance_ratio"),
  scale_output = TRUE
)
```

## Arguments

- res:

  RRI result.

- per_group:

  Compute by trajectory/group.

- group_cols:

  Group identifiers.

- id:

  Optional aligned identifiers.

- method:

  mean_neg_cor or variance_ratio.

- scale_output:

  For correlations maps `[-1, 1]` to `[0, 1]`; for variance cancellation
  truncates negative values to zero. A correlation score of 0.5 means
  zero mean correlation, not moderate biological compensation.

## Value

Diagnostic score, correlations and interpretation.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
  rri_compensation_index(res)
} # }
```
