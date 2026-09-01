# HRRI: Holobiont Redox Resilience Index — End-to-End Workflow

## Introduction

The **HRRI** package implements exploratory, multi-domain diagnostics
for describing how soil–plant–microbiome systems buffer and recover from
hydroclimatic redox disturbances (Ghotbi *et al.*, 2026). The
theoretical framing distinguishes four properties. They are not all
identifiable from a single observation curve:

| Property | Symbol | Interpretation |
|----|----|----|
| Capacity | $`Q`$ | Total electron-accepting/donating inventory (mmol e⁻ kg⁻¹) |
| Connectivity | $`\alpha`$ | Fraction of $`Q`$ electrochemically accessible to porewater |
| Kinetics | $`k`$ | Exchange rate (h⁻¹); encodes mineralogical crystallinity |
| Memory | $`M`$ | Legacy of prior disturbances; encoded in Fe-phase composition |

These combine through the accessible-capacity formula:

\$\$C\_{\rm acc} = \sum_i Q_i \cdot \alpha_i \cdot \left(1 - e^{-k_i
\tau}\right)\$\$

and aggregate into the Holobiont Redox Resilience Index:

``` math
\mathrm{RRI}_{it} = w_P \cdot P_{it} + w_S \cdot S_{it} + w_M \cdot M_{it}
```

where $`P`$ (Physiology), $`S`$ (Soil), and $`M`$ (Microbial) are domain
scores.

This vignette walks through the full workflow from in-silico data
generation to RRI computation, accessible-capacity estimation, and
recovery-signature metrics — all driven by
[`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md).

## In-silico Data Generation

[`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
is the package’s master simulator. It generates synthetic longitudinal
observations across all three holobiont domains.

Only Fe and Mn inventories carry closed-balance checks. Carbon,
nitrogen, sulfur and oxygen budgets are **not** closed, and all rate
parameters are illustrative rather than field-calibrated. Simulated
output supports software demonstration and falsifiable model checks, not
empirical ecological inference.

``` r

library(HRRI)
packageVersion("HRRI")
#> [1] '0.99.1'

## Compatibility shim -----------------------------------------------------
## rri_pipeline() is the convenience wrapper around rri_pipeline_st().
## If the *installed* HRRI predates the wrapper, define an equivalent local
## version so this vignette knits against either release. Reinstall the
## package (see README) to use the exported function directly.
if (!exists("rri_pipeline", mode = "function")) {
  message("Installed HRRI has no rri_pipeline(); using a vignette-local wrapper.")
  rri_pipeline <- function(dat = NULL, soil = NULL, plant = NULL,
                           micro = NULL, id = NULL,
                           domain_weights = c(Physio = 0.4, Soil = 0.35,
                                              Micro = 0.25), ...) {
    stopifnot(setequal(names(domain_weights),
                       c("Physio", "Soil", "Micro")))
    w <- domain_weights[c("Physio", "Soil", "Micro")]
    w <- w / sum(w)
    res <- rri_pipeline_st(
      ROS_flux = plant, Eh_stability = soil, micro_data = micro, id = id,
      w1 = unname(w[1]), w2 = unname(w[2]), w3 = unname(w[3]), ...
    )
    res$scores <- res$row_scores
    res
  }
}

## Reproducible 1-cycle flood-drain experiment
## n_plot=2, n_depth=2, n_plant=3, n_time=30 -> 360 rows
sim <- simulate_redox_holobiont(
  n_plot               = 2,
  n_depth              = 2,
  n_plant              = 3,
  n_time               = 30,
  p_micro              = 20,
  seed                 = 42,
  scenario             = "flood_drain",
  n_cycles             = 1,
  disturbance_strength = 0.70,
  history_strength     = 0.55,
  decoupling           = 0.20
)

## Top-level structure
names(sim)
#>  [1] "id"                   "forcing"              "latent_state"        
#>  [4] "soil_data"            "plant_data"           "micro_gene_abundance"
#>  [7] "micro_metat_counts"   "micro_metat_metadata" "micro_traits"        
#> [10] "fluxes"               "conservation_checks"  "ROS_flux"            
#> [13] "Eh_stability"         "micro_data"           "latent_truth"        
#> [16] "graph"                "metadata"
nrow(sim$id)                   # one row per plot × depth × plant × time
#> [1] 360
```

### Design identifiers

``` r

head(sim$id[, c("plot","depth","plant_id","time","cycle","phase","WFPS")])
#>   plot depth plant_id time cycle    phase      WFPS
#> 1   P1    D1   Plant1    1     1 baseline 0.6489965
#> 2   P2    D1   Plant1    1     1 baseline 0.5812486
#> 3   P1    D2   Plant1    1     1 baseline 0.7989965
#> 4   P2    D2   Plant1    1     1 baseline 0.7312486
#> 5   P1    D1   Plant2    1     1 baseline 0.6489965
#> 6   P2    D1   Plant2    1     1 baseline 0.5812486
```

### Soil geochemical outputs

``` r

head(sim$soil_data[, c("EAC","EDC","Cacc_EAC","Cacc_total","Cacc_fraction",
                        "FeIII_poor_crystalline_mmol_kg",
                        "FeII_mmol_kg","Eh","pH")])
#>        EAC      EDC Cacc_EAC Cacc_total Cacc_fraction
#> 1 266.3053 112.0382 110.4102   147.6319     0.3902060
#> 2 298.1276 119.3425 125.0644   162.0105     0.3880769
#> 3 337.3575 141.5196 107.5981   155.6755     0.3250845
#> 4 356.4329 146.4849 135.8578   188.5886     0.3749889
#> 5 279.5978 104.1723 113.5038   146.9287     0.3828560
#> 6 286.3268 115.3987 116.9400   152.2080     0.3788855
#>   FeIII_poor_crystalline_mmol_kg FeII_mmol_kg        Eh       pH
#> 1                       67.13652     17.26132 46.285671 6.657346
#> 2                       70.93895     17.91741 79.808173 6.544344
#> 3                       93.28308     24.82342 -4.972318 6.653453
#> 4                       99.62313     26.34155 19.816007 6.475302
#> 5                       67.21437     17.18353 56.597319 6.725813
#> 6                       70.95279     17.89324        NA 6.800352
```

### Fe mass-balance verification

``` r

## Maximum absolute error should be < 0.01 mmol kg-1
sim$conservation_checks
#>                                    check        value
#> 1 maximum_absolute_Fe_mass_balance_error 1.136868e-13
#> 2 maximum_absolute_Mn_mass_balance_error 1.421085e-14
#> 3                 minimum_simulated_pool 1.875881e-04
#> 4                 nonnegative_pool_check 1.000000e+00
```

### Plant physiology

``` r

head(sim$plant_data[, c("SPAD","FvFm","ROL","ROS_load","aerenchyma")])
#>       SPAD      FvFm       ROL  ROS_load aerenchyma
#> 1 41.64823 0.8149366 0.3165926 0.1962837  0.1169200
#> 2 41.74653 0.8053264 0.3388381 0.1817459  0.1169200
#> 3 43.48615 0.7906729 0.3122744 0.2424776  0.1313152
#> 4 43.51788 0.8002860 0.3490689 0.1989678  0.1169200
#> 5 43.79311 0.8167222 0.3719281 0.1988047  0.2164600
#> 6 42.65734 0.7955340 0.4006285 0.1607862  0.2164600
```

### Microbial functional genes

``` r

## 18 genes spanning Fe-cycling, denitrification, nitrification,
## methanogenesis, and sulfur cycling
colnames(sim$micro_gene_abundance)
#>  [1] "mtrA"     "omcS"     "cyc2"     "mnxG"     "narG"     "napA"    
#>  [7] "nirK"     "nirS"     "norB"     "nosZ"     "nrfA"     "amoA_AOA"
#> [13] "amoA_AOB" "nxrB"     "dsrA"     "dsrB"     "mcrA"     "pmoA"
summary(sim$micro_gene_abundance[, "mcrA"])   # methanogenesis gene
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   56717  140868  185855  195049  233312  578788
```

## Accessible-Capacity Estimation

[`rri_accessible_capacity()`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md)
computes \$C\_{\rm acc}\$ for arbitrary mineralogical reservoirs with
explicitly supplied accessibility and exchange-rate parameters. The
values below are illustrative model inputs, not estimates or validated
literature defaults.

``` r

## Subset one plot-depth unit for illustration
idx <- sim$id$plot == "P1" & sim$id$depth == "D1" & sim$id$plant_id == "Plant1"
sdf <- sim$soil_data[idx, ]

## Define reservoir specifications
## (Q_col names must match columns in sdf)
res_spec <- list(
  reactive_FeIII = list(
    Q_col = "FeIII_poor_crystalline_mmol_kg",
    alpha = "alpha_accept",   # column name: per-row connectivity
    k     = "k_accept_h",     # column name: per-row kinetics
    type  = "EAC"
  ),
  crystalline_FeIII = list(
    Q_col = "FeIII_crystalline_mmol_kg",
    alpha = 0.20,             # attenuated connectivity for crystalline phases
    k     = 0.008,            # h-1: slow exchange (goethite/hematite)
    type  = "EAC"
  ),
  FeII_pool = list(
    Q_col = "FeII_mmol_kg",
    alpha = "alpha_donate",
    k     = "k_donate_h",
    type  = "EDC"
  )
)

## tau = 24 h (diurnal event timescale)
cap <- rri_accessible_capacity(sdf, res_spec, tau = 24,
                                normalise = FALSE, return_components = TRUE)

## Per-component summary (returned because return_components = TRUE)
cap$components
#>           reservoir type  alpha k_per_tau_unit saturation_fraction mean_Q
#> 1    reactive_FeIII  EAC 0.5061         0.0672              0.7931 62.820
#> 2 crystalline_FeIII  EAC 0.2000         0.0080              0.1747 54.059
#> 3         FeII_pool  EDC 0.4807         0.0453              0.6602 20.946
#>   mean_contribution
#> 1            25.339
#> 2             1.889
#> 3             6.623

## Mean accessible vs. total inventory
cat("Mean Cacc_raw:", mean(cap$cacc_raw, na.rm=TRUE), "mmol e- kg-1\n")
#> Mean Cacc_raw: 33.85065 mmol e- kg-1
cat("Mean fraction :", mean(cap$cacc_fraction, na.rm=TRUE), "\n")
#> Mean fraction : 0.2455965
if ("ck_limited" %in% names(cap) && length(cap$ck_limited)) {
  cat("CK-limited rows:", sum(cap$ck_limited, na.rm=TRUE),
      "/", sum(!is.na(cap$ck_limited)), "classified rows\n")
} else {
  cat("CK-limited classification is not returned by this HRRI version.\n")
}
#> CK-limited rows: 0 / 0 classified rows
```

### Effect of event timescale τ

``` r

tau_vals <- c(1, 6, 24, 72, 168, 720)   # 1 h to 30 d
cacc_tau <- sapply(tau_vals, function(tt) {
  r <- rri_accessible_capacity(sdf, res_spec, tau = tt, normalise = FALSE)
  mean(r$cacc_raw, na.rm = TRUE)
})
data.frame(tau_h = tau_vals, Cacc_mean = round(cacc_tau, 2))
#>   tau_h Cacc_mean
#> 1     1      2.61
#> 2     6     13.47
#> 3    24     33.85
#> 4    72     45.84
#> 5   168     49.84
#> 6   720     52.63
```

## RRI Pipeline

[`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
integrates available observed domains into an exploratory composite. The
simulator’s hidden architecture columns are excluded. We use explicit
measured anchors to orient otherwise arbitrary PCA axes.

``` r

rri_out <- rri_pipeline(
  plant        = sim$ROS_flux,
  soil         = sim$Eh_stability,
  micro        = log1p(sim$micro_gene_abundance),
  id           = sim$id,
  mode         = "snapshot",
  scaling      = "pnorm",
  direction_anchor_phys = "FvFm",
  direction_anchor_soil = "Eh",
  direction_anchor_micro = "mtrA",
  domain_weights = c(Physio=0.35, Soil=0.40, Micro=0.25)
)

## Align once and reuse this identifier-complete table downstream.
rri_scored <- attach_hrri_ids(rri_out$row_scores, sim$id)
attr(rri_scored, "id_alignment")
#> [1] "observation keys"
summary(rri_scored$RRI)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2229  0.3861  0.5100  0.5124  0.6413  0.8111
head(rri_scored[, c("plot", "depth", "plant_id", "time",
                    "RRI", "Physio", "Soil", "Micro")])
#>   plot depth plant_id time       RRI    Physio      Soil     Micro
#> 1   P1    D1   Plant1    1 0.5325589 0.4003217 0.8663227 0.1836689
#> 2   P2    D1   Plant1    1 0.6679582 0.4420314 0.7785811 0.8072591
#> 3   P1    D2   Plant1    1 0.4450769 0.4103019 0.2531689 0.8008149
#> 4   P2    D2   Plant1    1 0.4132135 0.5230932 0.1604719 0.6637687
#> 5   P1    D1   Plant2    1 0.6320338 0.4016891 0.8255425 0.6449023
#> 6   P2    D1   Plant2    1 0.5843764 0.4715212 0.8202560 0.3649661
```

Alignment uses complete observation keys, or a shared unique `row_id`
when available. If neither is returned the pipeline must preserve input
row order — matching row counts alone do **not** establish alignment.
Shared identifier columns are checked for conflicts rather than silently
overwritten.

### Correlation with latent truth

``` r

## rri_scored is aligned to sim$id, and hence to its latent_truth vector.
truth <- sim$latent_truth
if (!is.numeric(truth) || length(truth) != nrow(rri_scored)) {
  stop("latent_truth must be a numeric vector with one value per sim$id row.")
}
ok <- is.finite(rri_scored$RRI) & is.finite(truth)
if (sum(ok) >= 3L && stats::sd(rri_scored$RRI[ok]) > 0 &&
    stats::sd(truth[ok]) > 0) {
  r_val <- stats::cor(rri_scored$RRI[ok], truth[ok])
} else {
  r_val <- NA_real_
  cat("Correlation unavailable: too few finite pairs or a constant vector.\n")
}
cat("r(RRI, latent truth):", round(r_val, 3), "\n")
#> r(RRI, latent truth): 0.452
```

This is a within-simulation association, not independent predictive
validation or evidence that the individual latent parameters have been
identified.

## Recovery Signatures

[`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md)
summarizes the RRI trajectory around a specified disturbance window. The
output fields depend on the installed package version; the complete
returned table is displayed below without assuming legacy names. A score
minimum or persistent departure alone does not establish biochemical
pathway truncation or alternative routing. Those interpretations require
independent process evidence.

The function requires **one row per time point per group**. Because the
pipeline returns one row per plant, we must first average RRI across
plants within each plot × depth × time cell before calling recovery
metrics.

``` r

## Step 1 — use the aligned score table created in the pipeline chunk.
## Step 2 — aggregate to one row per plot × depth × time (mean over plants).
## The data-frame method avoids formula-level complete-case filtering.
## An all-missing group remains NA; it is not replaced by zero.
rri_agg <- stats::aggregate(
  x = rri_scored["RRI"],
  by = rri_scored[c("plot", "depth", "time")],
  FUN = finite_mean
)
rri_agg <- rri_agg[order(rri_agg$plot, rri_agg$depth, rri_agg$time), ,
                   drop = FALSE]
rownames(rri_agg) <- NULL
stopifnot(!anyDuplicated(rri_agg[c("plot", "depth", "time")]))

## Step 3 — extract recovery signatures.
## Pass rri_agg directly (group columns are already inside it; no id= needed).
metrics <- rri_recovery_metrics(
  res           = rri_agg,
  time_col      = "time",
  group_cols    = c("plot","depth"),
  perturb_start = 8,
  perturb_end   = 18,
  rri_col       = "RRI",
  forcing_col   = NULL   # no measured forcing supplied
)

if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
  stop("rri_recovery_metrics() returned no nonempty recovery data frame.")
}
names(metrics)
#>  [1] "plot"                    "depth"                  
#>  [3] "baseline_rri"            "min_rri"                
#>  [5] "depth_min"               "depth_min_frac"         
#>  [7] "tau_lag"                 "k_recovery"             
#>  [9] "t_half"                  "overshoot"              
#> [11] "overshoot_frac"          "H_hysteresis"           
#> [13] "H_axis"                  "temporal_asymmetry"     
#> [15] "incomplete_return"       "incomplete_return_frac" 
#> [17] "displaced_plateau_flag"  "displaced_plateau_level"
#> [19] "alt_routing_flag"        "alt_routing_level"      
#> [21] "n_pre"                   "n_perturb"              
#> [23] "n_recovery"              "n_missing"              
#> [25] "n_fit"                   "fit_status"             
#> [27] "fit_r_squared"           "fit_start_time"         
#> [29] "final_observation_time"  "hysteresis_status"      
#> [31] "k"                       "H"                      
#> [33] "I"                       "H_abs"                  
#> [35] "H_norm"                  "I_norm"
metrics
#>   plot depth baseline_rri   min_rri  depth_min depth_min_frac tau_lag
#> 1   P1    D1    0.5587637 0.5045502 0.05421348     0.09702398       2
#> 2   P1    D2    0.3290230 0.2886692 0.04035387     0.12264754       2
#> 3   P2    D1    0.5893014 0.6032840 0.00000000     0.00000000      NA
#> 4   P2    D2    0.3430854 0.2762618 0.06682359     0.19477249       2
#>   k_recovery t_half overshoot overshoot_frac H_hysteresis      H_axis
#> 1         NA     NA 0.1742331      0.3118190           NA unavailable
#> 2         NA     NA 0.1379898      0.4193924           NA unavailable
#> 3         NA     NA 0.1621001      0.2750716           NA unavailable
#> 4         NA     NA 0.1794376      0.5230114           NA unavailable
#>   temporal_asymmetry incomplete_return incomplete_return_frac
#> 1                  1         0.1579226              0.2826286
#> 2                  1         0.0954700              0.2901621
#> 3                 NA         0.1470789              0.2495817
#> 4                  1         0.1206503              0.3516627
#>   displaced_plateau_flag displaced_plateau_level alt_routing_flag
#> 1                   TRUE               0.7166863               NA
#> 2                  FALSE                      NA               NA
#> 3                   TRUE               0.7363803               NA
#> 4                   TRUE               0.4637357               NA
#>   alt_routing_level n_pre n_perturb n_recovery n_missing n_fit
#> 1                NA     7        11         12         0     0
#> 2                NA     7        11         12         0     0
#> 3                NA     7        11         12         0     0
#> 4                NA     7        11         12         0     0
#>                       fit_status fit_r_squared fit_start_time
#> 1 insufficient_positive_deficits            NA             18
#> 2 insufficient_positive_deficits            NA             18
#> 3          no_resolvable_decline            NA             NA
#> 4 insufficient_positive_deficits            NA             18
#>   final_observation_time hysteresis_status  k  H         I H_abs H_norm
#> 1                     30     not_evaluated NA NA 0.1579226    NA     NA
#> 2                     30     not_evaluated NA NA 0.0954700    NA     NA
#> 3                     30     not_evaluated NA NA 0.1470789    NA     NA
#> 4                     30     not_evaluated NA NA 0.1206503    NA     NA
#>      I_norm
#> 1 0.2826286
#> 2 0.2901621
#> 3 0.2495817
#> 4 0.3516627
```

With `forcing_col = NULL`, a returned hysteresis-related field must be
interpreted according to the installed function’s documented definition;
it does not demonstrate a measured forcing–response loop. The limits 8
and 18 are example analysis boundaries, in the units of `time`; verify
them against the simulator’s event schedule before interpreting recovery
rates.

## Disturbance-History Sensitivity

Repeated-cycle simulations can test consequences of this simulator’s
stated rules. They do not establish a universal mineralogical ratchet.
In this version, crystallisation is continuous rather than restricted to
reoxidation events, and end-state EAC is not constrained to decline
monotonically with cycle count.

``` r

history <- do.call(rbind, lapply(1:4, function(nc) {
  z <- simulate_redox_holobiont(n_plot=1, n_depth=1, n_plant=2,
    n_time=30, p_micro=5, seed=99, n_cycles=nc,
    disturbance_strength=0.70)
  keep <- z$id$plant_id=="Plant1"
  data.frame(n_cycles=nc,
    EAC_end=tail(z$soil_data$EAC[keep],1),
    memory_end=tail(z$latent_state$memory[keep],1))
}))
history
#>   n_cycles  EAC_end memory_end
#> 1        1 301.2744          0
#> 2        2 296.6573          0
#> 3        3 292.8844          0
#> 4        4 290.9509          0
```

Interpret the direction and magnitude as a model sensitivity result. An
empirical claim about hydrological memory requires independent
observations.

## Session Information

``` r

sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.5.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] HRRI_0.99.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6       jsonlite_2.0.0     dplyr_1.2.1        compiler_4.5.1    
#>  [5] tidyselect_1.2.1   tidyr_1.3.2        jquerylib_0.1.4    systemfonts_1.3.2 
#>  [9] scales_1.4.0       textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0     
#> [13] ggplot2_4.0.3      R6_2.6.1           generics_0.1.4     igraph_2.3.3      
#> [17] knitr_1.51         htmlwidgets_1.6.4  tibble_3.3.1       desc_1.4.3        
#> [21] bslib_0.11.0       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.3.0       
#> [25] cachem_1.1.0       xfun_0.58          fs_2.1.0           sass_0.4.10       
#> [29] S7_0.2.2           otel_0.2.0         cli_3.6.6          pkgdown_2.2.1     
#> [33] magrittr_2.0.5     digest_0.6.39      grid_4.5.1         rstudioapi_0.18.0 
#> [37] lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5     glue_1.8.1        
#> [41] farver_2.1.2       ragg_1.5.2         purrr_1.2.2        rmarkdown_2.31    
#> [45] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.9
```

## References

Ghotbi, M., Kolody, B. C., Ghotbi, M., & Holtgrewe-Stukenbrock, E.
(2026). HRRI: a direction-aware R framework for quantifying
soil–plant–microbiome redox resilience across hydroclimatic disturbance
events *Methods in Ecology and Evolution*. Manuscript submitted.
