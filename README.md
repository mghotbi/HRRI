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
vignette("HRRI_gallery", package = "HRRI")
```

Walks through simulation, accessible-capacity estimation, domain scoring,
recovery signatures, and the mineralogical-ratchet disturbance-history
experiment.

## Citation

```r
citation("HRRI")
```

## Companion manuscripts

HRRI is the software component of three manuscripts, none yet published; two
are under review. They are listed here because the package implements what they
describe. Please cite the published version once available, and treat the
entries below as provisional until then.

**Software and diagnostics — the paper this package accompanies**
Ghotbi, M., Ghotbi, M., Komluski, J., & Holtgrewe-Stukenbrock, E. H. HRRI:
direction-aware diagnostics for soil–plant–microbiome redox recovery across
hydroclimatic disturbances. *In preparation.*
→ implemented by `rri_pipeline_st()`, `rri_property_scores()`,
`rri_recovery_metrics()`, `rri_accuracy()`.

**Theory — where the four hidden states come from**
Ghotbi, M., Kolody, B. C., Ghotbi, M., & Holtgrewe-Stukenbrock, E. A Theory of
Hydroclimatic Redox Resilience. *Submitted to Communications Earth &
Environment.*
→ the capacity–connectivity–kinetics–memory decomposition and the
accessible-capacity expression `C_acc(tau) = sum_j Q_j a_j (1 - exp(-k_j tau))`,
implemented by `rri_accessible_capacity()`.

**Mechanistic review — the biology the simulator encodes**
Ghotbi, M., Ghotbi, M., Mühling, K. H., & Stukenbrock, E. H. Rhizosphere redox
recovery after hydrological disturbances: mechanisms across the
soil–plant–microbiome continuum. *Submitted to Soil Biology & Biochemistry.*
→ the plant, microbial and mineralogical legacy terms in
`simulate_redox_holobiont()`.

None of these is required to use the package, and none is cited in `DESCRIPTION`:
CRAN asks that the `Description` field carry only references a reader can
actually retrieve, so it lists the published methods sources instead.

## Published methods this package builds on

| Reference | What HRRI takes from it |
|---|---|
| Sander, Hofstetter & Gorski (2015) *Environ. Sci. Technol.* 49:5862 [doi:10.1021/acs.est.5b00006](https://doi.org/10.1021/acs.est.5b00006) | Mediated electrochemical measurement of EAC and EDC — the capacity the index scores |
| Klüpfel, Piepenbrock, Kappler & Sander (2014) *Nat. Geosci.* 7:195 [doi:10.1038/ngeo2084](https://doi.org/10.1038/ngeo2084) | Electron-accepting capacity regenerates across repeated anoxic periods — the basis for treating memory as a state, not a trend |
| Thompson, Chadwick, Rancourt & Chorover (2006) *Geochim. Cosmochim. Acta* 70:1710 [doi:10.1016/j.gca.2005.12.005](https://doi.org/10.1016/j.gca.2005.12.005) | Iron-oxide crystallinity increases under redox oscillation — the mineralogical ratchet |
| Keiluweit, Wanzek, Kleber, Nico & Fendorf (2017) *Nat. Commun.* 8:1771 [doi:10.1038/s41467-017-01406-6](https://doi.org/10.1038/s41467-017-01406-6) | Anaerobic microsites persist in aerobic soil — why connectivity is separated from capacity |
| Lin (1989) *Biometrics* 45:255 [doi:10.2307/2532051](https://doi.org/10.2307/2532051) | Concordance correlation coefficient, reported by `rri_accuracy()` |
| Kobayashi & Salam (2000) *Agron. J.* 92:345 [doi:10.2134/agronj2000.922345x](https://doi.org/10.2134/agronj2000.922345x) | MSE partition into bias, variance mismatch and lack of correlation |


## License

MIT © Mitra Ghotbi. See [LICENSE](LICENSE).
