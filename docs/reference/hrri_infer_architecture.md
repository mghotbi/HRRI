# Fit a conditional capacity-recovery curve (legacy function name)

Fits y(t) = B \* (1 - A \* exp(-r \* t)) after a specified disturbance
peak. B is fixed from the observed baseline. A is a fractional recovery
deficit and r is a trajectory recovery rate. They are NOT the
accessibility alpha and reservoir exchange k in the accessible-capacity
model.

Tabulates estimates and independently supplied targets. There is no
default mapping from recovery-curve amplitude to accessibility or from
recovery rate to reservoir exchange kinetics.

## Usage

``` r
hrri_infer_architecture(
  Eh_stability,
  id = NULL,
  perturb_time = NULL,
  tau_unit = c("day", "hour", "week"),
  group_cols = c("plot", "depth"),
  fit_edc = FALSE,
  min_points = 3L,
  verbose = TRUE,
  control = list(),
  baseline_end = NULL
)

# S3 method for class 'hrri_arch'
print(x, ...)

# S3 method for class 'hrri_arch'
summary(object, ...)

# S3 method for class 'hrri_arch'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

validate_architecture(arch, params = NULL)

# S3 method for class 'hrri_arch_validation'
print(x, ...)
```

## Arguments

- Eh_stability:

  Data frame containing EAC and optionally EDC.

- id:

  Aligned identifiers containing time and grouping columns if absent.

- perturb_time:

  Disturbance peak time; NULL detects the EAC minimum. Detection is
  descriptive and can select noise or the last observation.

- tau_unit:

  Time unit: day, hour or week.

- group_cols:

  Columns defining one trajectory. Include plant_id when plant-level
  trajectories are separate; duplicate times are rejected.

- fit_edc:

  Also fit increasing EDC recovery toward its own baseline. This is
  inappropriate for a decreasing EDC trajectory; default FALSE.

- min_points:

  Minimum distinct post-peak times (at least three). This is a
  computational threshold, not proof of identifiability.

- verbose:

  Print fit status.

- control:

  Named list passed to stats::nls.control.

- baseline_end:

  Explicit end of the baseline window, strictly before the peak. NULL
  uses only the first observation time as a stated assumption.

- x:

  An hrri_arch object.

- ...:

  Additional method arguments.

- object:

  An hrri_arch object.

- row.names, optional:

  Standard data-frame method arguments.

- arch:

  An hrri_arch object with a keyed ground_truth data frame.

- params:

  Explicit strings of the form estimated_column:true_column.

## Value

An hrri_arch object containing estimates, conditional fit status, fit
objects and optional summaries of supplied latent columns. Legacy
columns Q_eac, alpha_eac, k_eac and M_eac remain for compatibility.
Prefer the aliases baseline_eac, deficit_fraction_eac, recovery_rate_eac
and terminal_ratio_eac. M_eac is a terminal-to-baseline ratio, not a
measure of causal memory. k_eac_h converts the trajectory rate to
reciprocal hours, not exchange kinetics.

Data frame with finite_pair flags and descriptive aggregate statistics.
The caller must establish that paired columns have the same definition
and units.

## Details

A nonlinear fit converging does not establish structural or practical
identifiability. Baseline uncertainty is not propagated into coefficient
SEs. With Cacc(t) = Q \* alpha \* (1 - exp(-k\*t)), Q and alpha cannot
be separated using that curve alone. This function does not solve that
inverse problem. Fit only a single recovery window. Nonmonotonic and
multi-event trajectories require a different model and residual
diagnostics.

## See also

rri_recovery_metrics, rri_accessible_capacity

## Examples

``` r
tt <- 0:12
yy <- ifelse(tt <= 2, 10, 10 * (1 - 0.7 * exp(-0.3 * (tt - 3))))
dat <- data.frame(plot = "P1", depth = "D1", time = tt, EAC = yy)
a <- hrri_infer_architecture(dat, perturb_time = 3, baseline_end = 2,
                             verbose = FALSE)
a$estimates
#>   plot depth t_perturb Q_eac alpha_eac k_eac k_eac_h alpha_eac_se     k_eac_se
#> 1   P1    D1         3    10       0.7   0.3  0.0125 3.654351e-12 1.964802e-12
#>       M_eac n_pre n_rec method_eac          ident_eac Q_edc alpha_edc k_edc
#> 1 0.9447268     3     9   nls_port fitted_conditional    NA        NA    NA
#>   k_edc_h M_edc ident_edc baseline_eac deficit_fraction_eac recovery_rate_eac
#> 1      NA    NA      <NA>           10                  0.7               0.3
#>   terminal_ratio_eac
#> 1          0.9447268
```
