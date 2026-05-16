# RRI

**RRI** quantifies holobiont redox resilience by integrating plant physiology,
soil redox chemistry, and microbial resilience into a unified, interpretable
Redox Resilience Index.

## Overview

RRI is an R package for computing, visualizing, and interpreting a
**Holobiont Redox Resilience Index (RRI)**. The framework treats resilience as
an emergent property of coordinated redox buffering across three interacting
domains:

- **Physio**: plant physiological buffering and stress response.
- **Soil**: soil redox chemistry and electron-acceptor stability.
- **Micro**: microbial abundance, function, traits, or network resilience.

The package supports spatial structure, temporal dynamics, directionally
identified latent scores, ternary visualization, sensitivity analysis, and
perturbation-recovery metrics.



<p align="center">
  <img src="https://github.com/user-attachments/assets/c64a2965-0413-41dc-952d-a7d0ce4017a4"
       alt="Conceptual overview of the RRI framework"
       width="900">
</p>

## Key features

- Integrates plant, soil, and microbial redox-resilience indicators.
- Computes per-sample RRI scores scaled to `[0, 1]`.
- Produces ternary-ready Physio-Soil-Micro compositional scores.
- Supports PCA, factor analysis, NMF, UMAP, WGCNA, and MFA-style workflows.
- Uses direction anchoring so higher scores consistently indicate stronger resilience.
- Supports microbial abundance, functional genes, MetaT expression, upregulation summaries, traits, and optional networks.
- Computes perturbation-recovery metrics including amplitude, lag, overshoot, incomplete recovery, recovery rate, hysteresis, and trajectory class.
- Provides diagnostic and sensitivity utilities for validation and robustness checks.

  
## Brief theory: redox resilience as buffering, connectivity, memory, and recovery

Redox resilience reflects the capacity of a soil-plant-microbiome system to
absorb, transmit, and recover from perturbations that alter electron flow.
In this framework, resilience is not treated as a single trait, but as an
emergent property of several interacting redox processes:

- **Capacity**: the size of electron-accepting and electron-donating pools,
  such as Fe/Mn oxides, organic matter, clays, humic substances, and other
  redox-active phases. High capacity allows the system to buffer electron
  imbalance before crossing critical redox thresholds.

- **Connectivity**: the physical and biochemical linkage among electron
  donors, electron acceptors, roots, microbes, minerals, and water-filled pore
  networks. Connectivity controls whether electrons, oxygen, substrates, and
  microbial processes remain spatially coupled or become locally isolated.

- **Kinetics**: the rate at which redox reactions proceed relative to
  hydrological or environmental forcing. Fast forcing combined with slow
  microbial, mineral, or diffusion-mediated responses can generate lag,
  transient disequilibrium, and redox pulses.

- **Memory**: the persistence of altered mineral, organic, microbial, or plant
  states after the original perturbation ends. Memory can arise from residual
  reduced minerals, redox-active organic matter, shifted microbial communities,
  enzyme pools, or plant acclimation.

- **Hysteresis**: path dependence in which wetting and drying, reduction and
  oxidation, or disturbance and recovery do not follow the same trajectory.
  Hysteresis indicates that the current state depends not only on present
  conditions, but also on the sequence of previous states.

- **Recovery**: the trajectory after perturbation, including return to the
  original state, transient overshoot, incomplete recovery, or transition to a
  new equilibrium.

RRI summarizes these interacting processes by integrating plant physiological
buffering, soil redox stability, and microbial resilience into a single
directionally aligned index. Perturbation-recovery functions then quantify how
the RRI trajectory changes through time, including amplitude, lag, overshoot,
incomplete recovery, recovery rate, hysteresis, and trajectory class.


## Installation

```r
install.packages("remotes")

remotes::install_github("mghotbi/RRI")
```

**To install with vignettes:**

```r
install.packages("remotes")
install.packages("BiocManager")

BiocManager::install(c("BiocStyle", "rmarkdown", "knitr"))

remotes::install_github(
  "mghotbi/RRI",
  build_vignettes = TRUE,
  dependencies = TRUE)

```

**Basic example**

```r
library(RRI)

sim <- simulate_redox_holobiont(seed = 1)

res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm",
  direction_phys = "auto",
  direction_anchor_phys = "FvFm",
  direction_soil = "auto",
  direction_anchor_soil = "Eh",
  direction_micro = "higher_is_better")

head(res$row_scores)
head(res$row_scores_comp)

```


**Adding microbial functional information
**
Users are not limited to microbial abundance data. Functional microbial
information can be added to the microbial domain when available.

```r
micro_features <- sim$micro_data

if (!is.null(sim$micro_traits)) {
  micro_features <- cbind(micro_features, sim$micro_traits)
}

if (!is.null(sim$gene_abundance)) {
  micro_features <- cbind(micro_features, sim$gene_abundance)
}

if (!is.null(sim$gene_log2fc)) {
  micro_features <- cbind(micro_features, sim$gene_log2fc)
}

res_full <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = micro_features,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm",
  direction_phys = "auto",
  direction_anchor_phys = "FvFm",
  direction_soil = "auto",
  direction_anchor_soil = "Eh",
  direction_micro = "higher_is_better")

```

The minimal and expanded workflows are both valid. The expanded workflow simply
uses a richer definition of microbial resilience.

Ternary visualization

```r
plot_RRI_ternary(res$row_scores_comp)

```

<img src="https://github.com/user-attachments/assets/abf56f94-582a-49d3-a28d-83087a01dd61" alt="Ternary visualization of Physio, Soil, and Micro contributions to RRI" width="900">

Perturbation-recovery metrics

RRI trajectories can be summarized using direction-aware perturbation-recovery
metrics.

```r
rec <- rri_recovery_metrics(
  res = res,
  id = sim$id,
  time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = 5,
  perturb_end = 7)

head(rec)

rri_metric_table()

```

Recovery metrics include perturbation amplitude, lag time, overshoot,
incomplete recovery, recovery rate, recovery fit diagnostics, recovery
timescale, half-recovery time, hysteresis, and trajectory class.


**Output structure**

The primary output from rri_pipeline_st() contains:

row_scores: absolute Physio, Soil, Micro, and RRI scores.
row_scores_comp: compositional Physio-Soil-Micro allocation scores.
meta: model settings and system-level RRI metadata.

The compositional output is directly compatible with plot_RRI_ternary().

**Additional utilities**

RRI includes helper functions for diagnostics, sensitivity analysis, recovery
metrics, and visualization:


```r
rri_latent_correlation(res, sim$latent_truth)
rri_compensation_index(res)
rri_sensitivity(res)

plot_RRI_ternary(res$row_scores_comp)
plot_rri_recovery_map(rec)
plot_rri_recovery_landscape(rec)


```

**Citation**

```r
citation("RRI")
```

Ghotbi, M. et al. RRI: A framework for holobiont-level redox resilience (manuscript in preparation).

**License**   GPL-3 © 2025 Mitra Ghotbi

**Contact**

Maintainer: Mitra Ghotbi
Email: mitra.ghotbi@gmail.com
ORCID: 0000-0001-9185-9993



