# Descriptive recovery speed score

Combines response lag and one rate descriptor. By default k and log(2)/k
are not counted as separate evidence. No causal exchange rate is
inferred.

## Usage

``` r
rri_kinetics_score(
  rec,
  forcing_window = NULL,
  lag_weight = 0.3,
  rate_weight = 0.7,
  halflife_weight = 0,
  invert_slow = TRUE
)
```

## Arguments

- rec:

  Recovery metric data frame.

- forcing_window:

  Positive duration in the same units as recovery time. With a duration,
  speed is k*T/(1+k*T), and lag score is 1/(1+lag/T). Without one,
  scores are cohort-relative min-max descriptions.

- lag_weight, rate_weight, halflife_weight:

  Non-negative weights.

- invert_slow:

  TRUE scores faster recovery higher; FALSE reverses all components.

## Value

Input with kinetics_score, component coverage, heuristic class, and lag
ratio.

## Examples

``` r
rec <- data.frame(tau_lag = c(2, 4, 1), k_recovery = c(0.3, 0.1, 0.5))
rri_kinetics_score(rec)$kinetics_score
#> [1] 0.55 0.00 1.00
```
