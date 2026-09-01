# HRRI ![](reference/figures/logo.png)

> Diagnostics for soil–plant–microbial redox recovery across
> hydroclimatic disturbance events.

**HRRI** provides transparent, assumption-explicit diagnostics for
longitudinal soil, plant and microbial observations spanning redox
disturbance and recovery. It computes stoichiometric oxygen demand,
event-window accessible electron capacity, six recovery signatures,
fixed-reference domain scores, and exploratory multiblock composites —
each with coverage diagnostics and documented limits on what may be
inferred.

------------------------------------------------------------------------

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mghotbi/HRRI", build_vignettes = TRUE)
```

## Quick start

``` r

library(HRRI)

## 1 — Generate illustrative trajectories (flood–drain, one cycle)
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 3, n_time = 30,
  scenario = "flood_drain", n_cycles = 1, seed = 42
)

## 2 — Score the three observed domains
res <- rri_pipeline(
  plant = sim$ROS_flux,
  soil  = sim$Eh_stability,
  micro = log1p(sim$micro_gene_abundance),
  id    = sim$id,
  direction_anchor_phys  = "FvFm",   # anchor otherwise-arbitrary PCA signs
  direction_anchor_soil  = "Eh",
  direction_anchor_micro = "mtrA"
)

## 3 — Extract recovery signatures (one row per trajectory)
scored <- attach_hrri_ids(res$row_scores, sim$id)
agg    <- aggregate(RRI ~ plot + depth + time, data = scored, FUN = mean)

rri_recovery_metrics(
  agg, time_col = "time", group_cols = c("plot", "depth"),
  perturb_start = 8, perturb_end = 18
)
```

## The four hidden-state controls

| Property | Symbol | Interpretation |
|:---|:--:|:---|
| Capacity | *Q* | Electron-accepting / donating inventory (mmol e⁻ kg⁻¹) |
| Connectivity | *α* | Fraction of *Q* connected to operative transfer pathways |
| Kinetics | *k* | Exchange rate (h⁻¹); encodes mineralogical crystallinity |
| Memory | *M* | Legacy of prior disturbance, carried in Fe-phase composition |

Accessible capacity over an event window of duration *τ*:

``` math
C_{\mathrm{acc}}(\tau) \;=\; \sum_i Q_i \, \alpha_i \left(1 - e^{-k_i \tau}\right)
```

## Function reference

**Simulation**

- [`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
  — mass-conserved Fe/Mn trajectories with plant and microbial
  observation models
- [`rri_simulation_demo()`](https://mghotbi.github.io/HRRI/reference/rri_simulation_demo.md)
  — reproducible end-to-end demonstration

**Capacity and stoichiometry**

- [`rri_accessible_capacity()`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md)
  — event-window *C*_(acc) for declared reservoirs
- [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  — complete-oxidation O₂ demand from reduced-pool inventories
- [`rri_capacity_index()`](https://mghotbi.github.io/HRRI/reference/rri_capacity_index.md)
  — oxidative-oriented soil feature composite

**Scoring**

- [`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  — convenience wrapper over available observed domains
- [`rri_pipeline_st()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md)
  — full-control multiblock interface
- [`rri_reference_scores()`](https://mghotbi.github.io/HRRI/reference/rri_reference_scores.md)
  — fixed-reference, externally anchored scoring
- [`attach_hrri_ids()`](https://mghotbi.github.io/HRRI/reference/attach_hrri_ids.md)
  — join design identifiers with explicit alignment checks

**Recovery and diagnostics**

- [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md)
  — lag, overshoot, hysteresis, depth, incomplete return, displaced
  plateau
- [`rri_memory_index()`](https://mghotbi.github.io/HRRI/reference/rri_memory_index.md),
  [`rri_kinetics_score()`](https://mghotbi.github.io/HRRI/reference/rri_kinetics_score.md),
  [`rri_connectivity_score()`](https://mghotbi.github.io/HRRI/reference/rri_connectivity_score.md),
  [`rri_compensation_index()`](https://mghotbi.github.io/HRRI/reference/rri_compensation_index.md)
- [`rri_property_scores()`](https://mghotbi.github.io/HRRI/reference/rri_property_scores.md)
  — property summary with provenance
- [`rri_sensitivity()`](https://mghotbi.github.io/HRRI/reference/rri_sensitivity.md)
  — sensitivity to domain aggregation weights

**Visualisation**

- [`plot_rri_timeseries()`](https://mghotbi.github.io/HRRI/reference/plot_rri_timeseries.md),
  [`plot_rri_state_space()`](https://mghotbi.github.io/HRRI/reference/plot_rri_state_space.md),
  [`plot_RRI_ternary()`](https://mghotbi.github.io/HRRI/reference/plot_RRI_ternary.md),
  [`plot_rri_recovery_map()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_map.md),
  [`plot_rri_properties()`](https://mghotbi.github.io/HRRI/reference/plot_rri_properties.md),
  [`plot_rri_validation()`](https://mghotbi.github.io/HRRI/reference/plot_rri_validation.md)

## Scope and limits

HRRI is deliberately conservative about inference. Please note:

- Latent axis directions are **not** biologically identifiable without
  justified anchors — supply `direction_anchor_*` arguments.
- A score decline does not identify pathway truncation; a displaced
  plateau does not establish alternative electron routing.
- Gene abundance indicates potential, not process rate.
- Simulator benchmarks measure agreement with a prescribed synthetic
  target. They are internal consistency checks, not empirical
  validation.
- `pnorm` scaling is a monotone transform, not a calibrated probability.

Each function’s help page states what its output does and does not
support.

## Vignette

``` r

vignette("HRRI_workflow", package = "HRRI")
```

Walks through simulation, accessible-capacity estimation, domain
scoring, recovery signatures, and the mineralogical-ratchet
disturbance-history experiment.

## Citation

``` r

citation("HRRI")
```

Ghotbi, M., Ghotbi, M., Guerreiro, M., Komluski, J., &
Holtgrewe-Stukenbrock, E. H. (2026). HRRI: A direction-aware R framework
for quantifying soil–plant–microbiome redox resilience across
hydroclimatic disturbance events. Manuscript submitted.

## License

MIT © Mitra Ghotbi. See
[LICENSE](https://mghotbi.github.io/HRRI/LICENSE).
