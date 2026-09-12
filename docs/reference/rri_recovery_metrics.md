# Descriptive recovery metrics for a single disturbance

Summarises decline and return of a higher-is-better score. A score
decline does not identify pathway truncation; a displaced plateau does
not establish alternative electron routing. Hysteresis is only reported
for a sufficiently closed, reversing forcing-response path. Temporal
deficit asymmetry is a separate diagnostic. Analyse repeated events
separately.

## Usage

``` r
rri_recovery_metrics(
  res,
  id = NULL,
  time_col = "time",
  group_cols = NULL,
  perturb_start,
  perturb_end,
  rri_col = "RRI",
  forcing_col = NULL,
  min_pts = 3L,
  lag_threshold = 0.05,
  plateau_window = 3L,
  plateau_tol = 0.1
)
```

## Arguments

- res:

  RRI object or data frame.

- id:

  Optional aligned identifiers; common columns must agree.

- time_col:

  Numeric time column; time must be unique within each group.

- group_cols:

  Columns identifying one longitudinal experimental unit.

- perturb_start, perturb_end:

  Finite start and end of one disturbance.

- rri_col:

  Numeric score column.

- forcing_col:

  Optional measured external forcing column, not a response proxy.

- min_pts:

  Minimum finite baseline and recovery observations.

- lag_threshold:

  Fraction of observed decline defining recovery onset.

- plateau_window:

  Number of final observations for plateau assessment.

- plateau_tol:

  Fractional terminal displacement defining a plateau flag.

## Value

One row per group, including diagnostic fit status and observation
counts. k is a log-linear fit of positive baseline deficits after the
observed minimum. It is a conditional trajectory descriptor, not a
mechanistic exchange rate. Legacy alt_routing fields are retained as NA;
use displaced_plateau_flag.

## Why `k_recovery` and `t_half` are often `NA`

The rate fit is deliberately conservative. It uses only recovery
observations that lie after the observed trough, strictly below the
pre-event baseline, and before the first crossing back through that
baseline, and it requires at least `min_pts` such points. Short recovery
windows routinely leave fewer, in which case `fit_status` is
`"insufficient_positive_deficits"` and both `k_recovery` and `t_half`
are returned as `NA` rather than being fitted to two points.

This is missingness by design, not failure. A rate estimated from a
handful of points spanning less than one recovery time constant is not
informative, and reporting it would invite over-interpretation. Always
read `k_recovery` together with `fit_status`, `n_fit` and
`fit_r_squared`.

If most trajectories return `NA`, extend the observation window rather
than lowering `min_pts`: as a rule of thumb the record should continue
for at least two or three times the expected recovery time after
`perturb_end`.

## Examples

``` r
## A window long enough for the rate fit to succeed. The event occupies
## days 12-22 of a 40-day record, leaving 18 recovery observations.
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 2, n_time = 40,
  p_micro = 10, seed = 2026
)

res <- suppressWarnings(rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
))

rec <- rri_recovery_metrics(
  res = res, id = sim$id, time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = 12, perturb_end = 22
)

## Inspect estimability before using any rate.
table(rec$fit_status)
#> 
#> insufficient_positive_deficits          no_resolvable_decline 
#>                              3                              5 
rec[, c("plot", "depth", "plant_id", "depth_min_frac",
        "k_recovery", "n_fit", "fit_status")]
#>   plot depth plant_id depth_min_frac k_recovery n_fit
#> 1   P1    D1   Plant1     0.00000000         NA     0
#> 2   P2    D1   Plant1     0.00000000         NA     0
#> 3   P1    D2   Plant1     0.15479833         NA     0
#> 4   P2    D2   Plant1     0.08137064         NA     0
#> 5   P1    D1   Plant2     0.00000000         NA     0
#> 6   P2    D1   Plant2     0.00000000         NA     0
#> 7   P1    D2   Plant2     0.07117166         NA     0
#> 8   P2    D2   Plant2     0.00000000         NA     0
#>                       fit_status
#> 1          no_resolvable_decline
#> 2          no_resolvable_decline
#> 3 insufficient_positive_deficits
#> 4 insufficient_positive_deficits
#> 5          no_resolvable_decline
#> 6          no_resolvable_decline
#> 7 insufficient_positive_deficits
#> 8          no_resolvable_decline

## Contrast: a short window leaves too few usable points and the rate
## is correctly withheld.
short <- simulate_redox_holobiont(
  n_plot = 1, n_depth = 1, n_plant = 2, n_time = 12,
  p_micro = 5, seed = 1
)
res_s <- suppressWarnings(rri_pipeline_st(
  ROS_flux = short$ROS_flux, Eh_stability = short$Eh_stability,
  micro_data = short$micro_data, id = short$id
))
rec_s <- rri_recovery_metrics(
  res = res_s, id = short$id, time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = 5, perturb_end = 7
)
table(rec_s$fit_status)
#> 
#> insufficient_positive_deficits 
#>                              2 
```
