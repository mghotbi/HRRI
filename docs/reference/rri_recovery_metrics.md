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

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux,
                      id = sim$id)
  rri_recovery_metrics(res, time_col = "time",
    group_cols = c("plot", "depth", "plant_id"),
    perturb_start = 20, perturb_end = 30)
} # }
```
