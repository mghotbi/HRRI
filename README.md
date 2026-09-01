# HRRI

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.3](https://img.shields.io/badge/R-%3E%3D%204.3-276DC3.svg)](https://cran.r-project.org/)
<!-- badges: end -->

> Diagnostics for soil–plant–microbial redox recovery across hydroclimatic
> disturbance events.

**HRRI** provides transparent, assumption-explicit diagnostics for longitudinal
soil, plant and microbial observations spanning redox disturbance and recovery.
It computes stoichiometric oxygen demand, event-window accessible electron
capacity, six recovery signatures, fixed-reference domain scores, and
exploratory multiblock composites — each with coverage diagnostics and
documented limits on what may be inferred.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("mghotbi/HRRI", build_vignettes = TRUE)
```

## Quick start

```r
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
|:---|:---:|:---|
| Capacity     | *Q* | Electron-accepting and electron-donating inventory available within the system (mmol e⁻ kg⁻¹) |
| Connectivity | *α* | Fraction of capacity functionally connected to active electron-transfer pathways |
| Kinetics     | *k* | Characteristic rate of electron exchange under physicochemical and biological constraints (h⁻¹) |
| Memory       | *M* | Legacy of prior disturbances retained through persistent biogeochemical, microbial and physiological states that influence future system responses |

Accessible capacity over an event window of duration *τ*:

$$C_{\mathrm{acc}}(\tau) \;=\; \sum_i Q_i \, \alpha_i \left(1 - e^{-k_i \tau}\right)$$

## Function reference

**Simulation**

- `simulate_redox_holobiont()` — mass-conserved Fe/Mn trajectories with plant and microbial observation models
- `rri_simulation_demo()` — reproducible end-to-end demonstration

**Capacity and stoichiometry**

- `rri_accessible_capacity()` — event-window *C*<sub>acc</sub> for declared reservoirs
- `rri_o2_demand()` — complete-oxidation O₂ demand from reduced-pool inventories
- `rri_capacity_index()` — oxidative-oriented soil feature composite

**Scoring**

- `rri_pipeline()` — convenience wrapper over available observed domains
- `rri_pipeline_st()` — full-control multiblock interface
- `rri_reference_scores()` — fixed-reference, externally anchored scoring
- `attach_hrri_ids()` — join design identifiers with explicit alignment checks

**Recovery and diagnostics**

- `rri_recovery_metrics()` — lag, overshoot, hysteresis, depth, incomplete return, displaced plateau
- `rri_memory_index()`, `rri_kinetics_score()`, `rri_connectivity_score()`, `rri_compensation_index()`
- `rri_property_scores()` — property summary with provenance
- `rri_sensitivity()` — sensitivity to domain aggregation weights

**Visualisation**

- `plot_rri_timeseries()`, `plot_rri_state_space()`, `plot_RRI_ternary()`,
  `plot_rri_recovery_map()`, `plot_rri_properties()`, `plot_rri_validation()`

## Scope and limits

HRRI is deliberately conservative about inference. Please note:

- Latent axis directions are **not** biologically identifiable without justified
  anchors — supply `direction_anchor_*` arguments.
- A score decline does not identify pathway truncation; a displaced plateau does
  not establish alternative electron routing.
- Gene abundance indicates potential, not process rate.
- Simulator benchmarks measure agreement with a prescribed synthetic target.
  They are internal consistency checks, not empirical validation.
- `pnorm` scaling is a monotone transform, not a calibrated probability.

Each function's help page states what its output does and does not support.

## Vignette

```r
vignette("HRRI_workflow", package = "HRRI")
```

Walks through simulation, accessible-capacity estimation, domain scoring,
recovery signatures, and the mineralogical-ratchet disturbance-history
experiment.

## Citation

```r
citation("HRRI")
```

Ghotbi, M., Ghotbi, M., Guerreiro, M., Komluski, J., & Holtgrewe-Stukenbrock, E. H. (2026). 
HRRI: A direction-aware R framework for quantifying soil–plant–microbiome redox resilience 
across hydroclimatic disturbance events. Manuscript submitted.


## License

MIT © Mitra Ghotbi. See [LICENSE](LICENSE).
