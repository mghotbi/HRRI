# Plot a recovery landscape from RRI perturbation-recovery metrics

Visualises user-level perturbation-recovery metrics from
[`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md).
Each row is one trajectory and each column is a core resilience metric.
Values are scaled within metric columns to make amplitudes, overshoot,
incomplete recovery, and recovery times comparable.

## Usage

``` r
plot_rri_recovery_landscape(
  rec,
  group_cols = c("plot", "depth", "plant_id"),
  metrics = c("A_norm", "O_norm", "I_norm", "k", "tau_r", "t_half"),
  order_by = "I_norm",
  base_size = 12
)
```

## Arguments

- rec:

  A data frame returned by
  [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md).

- group_cols:

  Character vector identifying trajectory labels.

- metrics:

  Character vector of recovery metric columns to plot.

- order_by:

  Character scalar. Metric used to order trajectories.

- base_size:

  Numeric. Base font size.

## Value

A `ggplot` object.

## Examples

``` r
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 3,
  n_plant = 2,
  n_time = 12,
  p_micro = 20,
  seed = 109
)

res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h

rec <- rri_recovery_metrics(
  res = res,
  id = sim$id,
  time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = 5,
  perturb_end = 7
)

head(rec)
#>   plot depth plant_id baseline_rri   min_rri depth_min depth_min_frac tau_lag
#> 1   P1    D1   Plant1    0.3143189 0.1583739 0.1559451      0.4961364       2
#> 2   P2    D1   Plant1    0.5070936 0.2082781 0.2988155      0.5892709       1
#> 3   P1    D2   Plant1    0.5105035 0.2354041 0.2750993      0.5388785      NA
#> 4   P2    D2   Plant1    0.7260577 0.3609971 0.3650605      0.5027983       2
#> 5   P1    D3   Plant1    0.6753192 0.3943502 0.2809690      0.4160537       1
#> 6   P2    D3   Plant1    0.8497697 0.4842407 0.3655291      0.4301507       1
#>   k_recovery   t_half  overshoot overshoot_frac H_hysteresis      H_axis
#> 1         NA       NA 0.00000000     0.00000000           NA unavailable
#> 2 0.03074131 22.54774 0.00000000     0.00000000           NA unavailable
#> 3         NA       NA 0.00000000     0.00000000           NA unavailable
#> 4         NA       NA 0.00000000     0.00000000           NA unavailable
#> 5         NA       NA 0.00000000     0.00000000           NA unavailable
#> 6         NA       NA 0.03657464     0.04304065           NA unavailable
#>   temporal_asymmetry incomplete_return incomplete_return_frac
#> 1         -0.3378744        -0.1382381             -0.4398019
#> 2         -0.3671227        -0.2916808             -0.5752012
#> 3         -0.3441857        -0.2555410             -0.5005667
#> 4         -0.2780719        -0.2854339             -0.3931284
#> 5         -0.2472709        -0.2029588             -0.3005375
#> 6         -0.1513147        -0.2205439             -0.2595337
#>   displaced_plateau_flag displaced_plateau_level alt_routing_flag
#> 1                   TRUE               0.1760809               NA
#> 2                   TRUE               0.2154128               NA
#> 3                  FALSE                      NA               NA
#> 4                  FALSE                      NA               NA
#> 5                  FALSE                      NA               NA
#> 6                  FALSE                      NA               NA
#>   alt_routing_level n_pre n_perturb n_recovery n_missing n_fit
#> 1                NA     4         3          5         0     2
#> 2                NA     4         3          5         0     3
#> 3                NA     4         3          5         0     1
#> 4                NA     4         3          5         0     2
#> 5                NA     4         3          5         0     1
#> 6                NA     4         3          5         0     1
#>                       fit_status fit_r_squared fit_start_time
#> 1 insufficient_positive_deficits            NA             11
#> 2 fitted_conditional_exponential     0.8905163             10
#> 3 insufficient_positive_deficits            NA             12
#> 4 insufficient_positive_deficits            NA             11
#> 5 insufficient_positive_deficits            NA             12
#> 6 insufficient_positive_deficits            NA             12
#>   final_observation_time hysteresis_status          k  H         I H_abs H_norm
#> 1                     12     not_evaluated         NA NA 0.1382381    NA     NA
#> 2                     12     not_evaluated 0.03074131 NA 0.2916808    NA     NA
#> 3                     12     not_evaluated         NA NA 0.2555410    NA     NA
#> 4                     12     not_evaluated         NA NA 0.2854339    NA     NA
#> 5                     12     not_evaluated         NA NA 0.2029588    NA     NA
#> 6                     12     not_evaluated         NA NA 0.2205439    NA     NA
#>      I_norm
#> 1 0.4398019
#> 2 0.5752012
#> 3 0.5005667
#> 4 0.3931284
#> 5 0.3005375
#> 6 0.2595337
plot_rri_recovery_landscape(rec)
#> Error: Missing metric columns: A_norm, O_norm, tau_r
```
