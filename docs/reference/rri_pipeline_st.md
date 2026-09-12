# Exploratory domain-score integration (legacy interface)

Integrates plant, soil and microbial latent scores. Direction is not
biologically identifiable without justified anchors. All three domains
may be incomplete. Available positive domain weights are renormalized
per row. Use rri_reference_scores for externally anchored,
fixed-reference comparisons.

## Usage

``` r
rri_pipeline_st(
  ROS_flux = NULL,
  Eh_stability = NULL,
  micro_data = NULL,
  graph = NULL,
  id = NULL,
  time_col = NULL,
  group_cols = NULL,
  mode = c("snapshot", "rolling", "event"),
  window = 3,
  align = c("right", "center", "left"),
  event_col = NULL,
  baseline_label = "pre",
  recovery_labels = "recovery",
  alpha_micro = 0.5,
  method_phys = "pca",
  method_soil = "pca",
  method_micro = "pca",
  direction_phys = c("auto", "higher_is_better", "lower_is_better"),
  direction_soil = c("auto", "higher_is_better", "lower_is_better"),
  direction_micro = c("auto", "higher_is_better", "lower_is_better"),
  direction_anchor_phys = NULL,
  direction_anchor_soil = NULL,
  direction_anchor_micro = NULL,
  scale_by = NULL,
  network_agg = c("equation", "mean"),
  w1 = 0.4,
  w2 = 0.35,
  w3 = 0.25,
  add_coupling = FALSE,
  coupling_weight = 0,
  coupling_fun = c("geometric_mean", "agreement"),
  norm_method = NULL,
  reducer = c("per_domain", "mfa"),
  scaling = c("minmax_legacy", "pnorm"),
  comp_space = c("closure_legacy", "clr"),
  ref_stats = NULL,
  add_compensation = FALSE,
  compensation_weight = 0
)
```

## Arguments

- ROS_flux:

  Data frame of plant physiological variables (rows = samples).

- Eh_stability:

  Data frame of soil redox chemistry variables (rows = samples).

- micro_data:

  Optional data frame of microbial abundance or functional features.

- graph:

  Optional `igraph` object or list of `igraph` objects representing
  microbial network structure.

- id:

  Optional data frame describing experimental design (same number of
  rows as inputs).

- time_col:

  Optional character. Name of time column in `id`.

- group_cols:

  Optional character vector of grouping variables in `id`.

- mode:

  Character. One of `"snapshot"`, `"rolling"`, or `"event"`.

- window:

  Integer \>= 2. Rolling window size (for mode = "rolling").

- align:

  Character. Alignment rule for rolling window: `"right"`, `"center"`,
  or `"left"`.

- event_col:

  Optional character. Column in `id` identifying event phases.

- baseline_label:

  Character. Label identifying baseline phase.

- recovery_labels:

  Character vector identifying recovery phases.

- alpha_micro:

  Numeric between 0 and 1 controlling blending of microbial abundance
  and network components.

- method_phys:

  Character. Reduction method for plant block.

- method_soil:

  Character. Reduction method for soil block.

- method_micro:

  Character. Reduction method for microbial block.

- direction_phys:

  Character. Orientation rule for plant latent dimension.

- direction_soil:

  Character. Orientation rule for soil latent dimension.

- direction_micro:

  Character. Orientation rule for microbial latent dimension.

- direction_anchor_phys:

  Optional character. Anchor variable for plant orientation.

- direction_anchor_soil:

  Optional character. Anchor variable for soil orientation.

- direction_anchor_micro:

  Optional character. Anchor variable for microbial orientation.

- scale_by:

  Optional character vector of grouping variables used for scaling.

- network_agg:

  Character. Network aggregation method: `"equation"` or `"mean"`.

- w1:

  Numeric weight for plant domain.

- w2:

  Numeric weight for soil domain.

- w3:

  Numeric weight for microbial domain. Must sum with w1 and w2 to 1.

- add_coupling:

  Logical. If TRUE, adds cross-domain coherence term.

- coupling_weight:

  Numeric between 0 and 1 controlling weight of coupling term.

- coupling_fun:

  Character. Coupling function: `"geometric_mean"` or `"agreement"`.

- norm_method:

  Optional character. If provided, overrides block-specific methods.

- reducer:

  Character. Reduction strategy: `"per_domain"` or `"mfa"`.

- scaling:

  Character. Scaling rule: `"minmax_legacy"` or `"pnorm"`.

- comp_space:

  Character. Compositional projection method: `"closure_legacy"` or
  `"clr"`.

- ref_stats:

  Optional list of reference statistics used for scaling.

- add_compensation:

  Logical. If TRUE, includes covariance-based compensation term.

- compensation_weight:

  Numeric between 0 and 1 controlling compensation weight.

## Value

RRI object; identifiers accompany scores and rolling output retains
original input order. Stochastic and advanced reducers need separate
validation.

## Details

MFA is disabled pending a validated implementation. Scaling statistics
do not freeze PCA/FA loadings, so ref_stats is not a trained prediction
model. The CLR round trip changes display coordinates only: inversion
returns closure. Grouping does not imply within-group scaling; request
scale_by explicitly. Missing data are median-imputed for exploratory
reduction, not corrected for MNAR. Event scores are descriptive products
of resistance and reference proximity; baseline/recovery label defaults
must be matched to the supplied data.

## Examples

``` r
# \donttest{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline_st(sim$ROS_flux, sim$Eh_stability, id = sim$id)
#> Warning: Unanchored latent axes have arbitrary signs; RRI is exploratory, not directionally validated resilience.
#> Warning: Excluding simulator-derived hidden columns from scoring: Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction, net_oxidative_balance, alpha_accept, alpha_donate, k_accept_h, k_donate_h
  head(res$row_scores)
#>   plot depth plant_id time row_id      unit_id history_pair history    scenario
#> 1   P1    D1   Plant1    1      1 P1.D1.Plant1  P1_D1_Pair1   naive flood_drain
#> 2   P2    D1   Plant1    1      2 P2.D1.Plant1  P2_D1_Pair1   naive flood_drain
#> 3   P3    D1   Plant1    1      3 P3.D1.Plant1  P3_D1_Pair1   naive flood_drain
#> 4   P4    D1   Plant1    1      4 P4.D1.Plant1  P4_D1_Pair1   naive flood_drain
#> 5   P1    D2   Plant1    1      5 P1.D2.Plant1  P1_D2_Pair1   naive flood_drain
#> 6   P2    D2   Plant1    1      6 P2.D2.Plant1  P2_D2_Pair1   naive flood_drain
#>   rescue cycle    phase event_intensity      WFPS water_table_cm    Physio
#> 1   none     1 baseline      0.06572853 0.5780600      -5.464198 0.3104766
#> 2   none     1 baseline      0.06572853 0.6064134      -7.448936 0.3368434
#> 3   none     1 baseline      0.06572853 0.5707389      -4.951720 0.3281934
#> 4   none     1 baseline      0.06572853 0.6558207     -10.907448 0.3670551
#> 5   none     1 baseline      0.06572853 0.7280600      -7.964198 0.2724980
#> 6   none     1 baseline      0.06572853 0.7564134      -9.948936 0.3439660
#>        Soil Micro       RRI domain_coverage n_domains Micro_abundance
#> 1 0.1043981    NA 0.2143066            0.75         2              NA
#> 2 0.1005309    NA 0.2265642            0.75         2              NA
#> 3 0.1566488    NA 0.2481392            0.75         2              NA
#> 4 0.0661552    NA 0.2266351            0.75         2              NA
#> 5 0.6788061    NA 0.4621085            0.75         2              NA
#> 6 0.5062489    NA 0.4196980            0.75         2              NA
#>   Micro_network Micro_mfa
#> 1            NA        NA
#> 2            NA        NA
#> 3            NA        NA
#> 4            NA        NA
#> 5            NA        NA
#> 6            NA        NA
# }
```
