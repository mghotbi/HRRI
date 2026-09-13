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
vignette(package = "HRRI")     
browseVignettes("HRRI")

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
are under review. They are listed because the package implements what they
describe. Please cite the published version once available.

**Software and framework — the paper this package accompanies**
Ghotbi, M., Ghotbi, M., Komluski, J., & Holtgrewe-Stukenbrock, E. H. HRRI: a
framework for diagnosing redox recovery in soil–plant–microbiome systems.
*In preparation.*
→ implemented by `rri_pipeline_st()`, `rri_property_scores()`,
`rri_recovery_metrics()`, `rri_accuracy()`.

**Theory — where the four hidden states come from**
Ghotbi, M., Kolody, B. C., Ghotbi, M., & Holtgrewe-Stukenbrock, E. A Theory of
Hydroclimatic Redox Resilience. *Submitted to Communications Earth &
Environment.*
→ the capacity–connectivity–kinetics–memory decomposition and the
accessible-capacity expression, implemented by `rri_accessible_capacity()`.

**Mechanistic review — the biology the simulator encodes**
Ghotbi, M., Ghotbi, M., Mühling, K. H., & Stukenbrock, E. H. Rhizosphere redox
recovery after hydrological disturbances: mechanisms across the
soil–plant–microbiome continuum. *Submitted to Soil Biology & Biochemistry.*
→ the plant, microbial and mineralogical legacy terms in
`simulate_redox_holobiont()`.

None is required to use the package, and none is cited in `DESCRIPTION`: CRAN
asks that the `Description` field carry only references a reader can retrieve,
so it lists the published methods sources instead.

## Published methods this package builds on

Every DOI below was resolved against Crossref before being listed.

**Measuring the four quantities**

| Reference | What HRRI takes from it |
|---|---|
| Sander, Hofstetter & Gorski (2015) *Environ. Sci. Technol.* 49:5862 [doi:10.1021/acs.est.5b00006](https://doi.org/10.1021/acs.est.5b00006) | Mediated electrochemical measurement of EAC and EDC in electron equivalents — the capacity *Q* the index scores, rather than an elemental concentration |
| Dorau et al. (2022) *Eur. J. Soil Sci.* 73:e13165 [doi:10.1111/ejss.13165](https://doi.org/10.1111/ejss.13165) | *Connected* air-filled porosity, not total air content, governs the shift toward oxidising conditions — the measurement behind α |
| Peiffer et al. (2021) *Nat. Geosci.* 14:264–272 [doi:10.1038/s41561-021-00742-z](https://doi.org/10.1038/s41561-021-00742-z) | Framework coupling redox-active compound pools to hydrological forcing — why inventory and event timescale must be carried separately |

**Why bulk state variables are not enough**

| Reference | What HRRI takes from it |
|---|---|
| Rooney et al. (2024) *Commun. Earth Environ.* 5 [doi:10.1038/s43247-024-01927-1](https://doi.org/10.1038/s43247-024-01927-1) | Redox processes decouple from soil saturation — moisture recovery does not imply redox recovery |
| Keiluweit et al. (2017) *Nat. Commun.* 8:1771 [doi:10.1038/s41467-017-01406-6](https://doi.org/10.1038/s41467-017-01406-6) | Anaerobic microsites persist in aerobic soil — why connectivity is separated from capacity rather than folded into it |
| Angle et al. (2017) *Nat. Commun.* 8:1567 [doi:10.1038/s41467-017-01753-4](https://doi.org/10.1038/s41467-017-01753-4) | Methanogenesis in oxygenated soils — reducing metabolism where a bulk measurement would not predict it |

**Memory as a state, not a trend**

| Reference | What HRRI takes from it |
|---|---|
| Thompson et al. (2006) *Geochim. Cosmochim. Acta* 70:1710–1727 [doi:10.1016/j.gca.2005.12.005](https://doi.org/10.1016/j.gca.2005.12.005) | Iron-oxide crystallinity increases under redox oscillation — the mineralogical ratchet |
| Aeppli et al. (2019) *Environ. Sci. Technol.* 53:3568–3578 [doi:10.1021/acs.est.8b07190](https://doi.org/10.1021/acs.est.8b07190) | Reducibility falls as ferrihydrite transforms abiotically to goethite and magnetite — why the ratchet lowers the ceiling |
| Aeppli et al. (2019) *Environ. Sci. Technol.* 53:8736–8746 [doi:10.1021/acs.est.9b01299](https://doi.org/10.1021/acs.est.9b01299) | The same loss of reducibility under *microbial* reductive dissolution — the ratchet is not solely abiotic |
| Klüpfel et al. (2014) *Nat. Geosci.* 7:195–200 [doi:10.1038/ngeo2084](https://doi.org/10.1038/ngeo2084) | Humic electron-accepting capacity is fully regenerable across repeated anoxic periods — a cycle, not a ratchet, so the two components cannot share one decay term |
| Meisner et al. (2021) *ISME J.* 15:1207–1221 [doi:10.1038/s41396-020-00844-3](https://doi.org/10.1038/s41396-020-00844-3) | Microbial legacies differ by disturbance type — why soil, plant and microbial legacies are returned separately |

**Limits on what may be inferred**

| Reference | What HRRI takes from it |
|---|---|
| Louca et al. (2018) *Nat. Ecol. Evol.* 2:936–943 [doi:10.1038/s41559-018-0519-1](https://doi.org/10.1038/s41559-018-0519-1) | Functional redundancy decouples taxonomy from function — gene abundance indicates potential, not process rate |
| Gloor et al. (2017) *Front. Microbiol.* 8:2224 [doi:10.3389/fmicb.2017.02224](https://doi.org/10.3389/fmicb.2017.02224) | Microbiome data are compositional — why a log-ratio workflow must be declared rather than assumed |

**Statistics**

| Reference | What HRRI takes from it |
|---|---|
| Lin (1989) *Biometrics* 45:255–268 [doi:10.2307/2532051](https://doi.org/10.2307/2532051) | Concordance correlation coefficient — agreement, not merely association, reported by `rri_accuracy()` |
| Kobayashi & Salam (2000) *Agron. J.* 92:345–352 [doi:10.2134/agronj2000.922345x](https://doi.org/10.2134/agronj2000.922345x) | Exact partition of mean squared error into bias, variance mismatch and lack of correlation |

**Rate context for the O₂ demand calculation**

| Reference | What HRRI takes from it |
|---|---|
| Stumm & Lee (1961) *Ind. Eng. Chem.* 53:143–146 [doi:10.1021/ie50614a030](https://doi.org/10.1021/ie50614a030) | Fe(II) oxygenation kinetics — stoichiometric demand is pH-independent, the rate is not |
| Millero, Sotolongo & Izaguirre (1987) *Geochim. Cosmochim. Acta* 51:793–801 [doi:10.1016/0016-7037(87)90093-7](https://doi.org/10.1016/0016-7037(87)90093-7) | The ~100-fold rate increase per unit pH that separates a ceiling from a realised consumption |

## License

MIT © Mitra Ghotbi. See [LICENSE](LICENSE).
