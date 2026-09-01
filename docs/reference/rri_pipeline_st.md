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
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(seed = 1)
  res <- rri_pipeline_st(sim$ROS_flux, sim$Eh_stability, id = sim$id)
  head(res$row_scores)
} # }
```
