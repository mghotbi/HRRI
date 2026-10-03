## ----setup, include=FALSE-----------------------------------------------------
# Select a compatible installed HRRI before any calculations.
# Check the required API, rather than rejecting a compatible 1.0.7 build.
.hrri_ready <- function(ns) {
  required <- c("plot_rri_framework", "plot_rri_identifiability",
                "plot_rri_recovery_diagnostics", "plot_rri_timeseries",
                "plot_rri_accuracy")
  all(required %in% getNamespaceExports(ns)) &&
    all(c("forcing_threshold", "time_label") %in%
          names(formals(getExportedValue("HRRI", "plot_rri_timeseries")))) &&
    "cluster_label" %in%
      names(formals(getExportedValue("HRRI", "plot_rri_accuracy")))
}
.hrri_select <- function() {
  if ("HRRI" %in% loadedNamespaces() && .hrri_ready(asNamespace("HRRI")))
    return(invisible(getNamespaceInfo(asNamespace("HRRI"), "path")))

  candidates <- unique(c(path.expand("~/R/HRRI-1.0.8"),
                         path.expand("~/R/HRRI-repaired"), .libPaths()))
  candidates <- candidates[file.exists(file.path(candidates, "HRRI", "DESCRIPTION"))]
  if ("package:HRRI" %in% search())
    detach("package:HRRI", unload = FALSE, character.only = TRUE)
  if ("HRRI" %in% loadedNamespaces()) unloadNamespace("HRRI")

  for (lib in candidates) {
    ns <- tryCatch(loadNamespace("HRRI", lib.loc = lib), error = function(e) NULL)
    if (!is.null(ns) && .hrri_ready(ns))
      return(invisible(getNamespaceInfo(ns, "path")))
    if ("HRRI" %in% loadedNamespaces()) unloadNamespace("HRRI")
  }
  stop("No installed HRRI provides the functions and arguments used by this vignette. ",
       "Install the supplied HRRI_1.0.8.tar.gz. Libraries checked: ",
       paste(candidates, collapse = "; "), call. = FALSE)
}
.hrri_selected_path <- .hrri_select()
library(HRRI, lib.loc = dirname(.hrri_selected_path))
rm(.hrri_ready, .hrri_select)
knitr::opts_chunk$set(dev = "png", 
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 7,
  fig.height = 4.2,
  fig.align  = "center",
  dpi        = 200,
  out.width  = "100%",
  message    = FALSE,
  warning    = FALSE
)

## Consistent figure styling across the vignette -------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  old_theme <- ggplot2::theme_set(
    ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(linewidth = 0.3,
                                                 colour = "#dde3e1"),
        axis.title       = ggplot2::element_text(colour = "#4a5451"),
        axis.text        = ggplot2::element_text(colour = "#4a5451"),
        strip.text       = ggplot2::element_text(face = "bold",
                                                 colour = "#1c2321"),
        plot.title       = ggplot2::element_text(face = "bold",
                                                 colour = "#1c2321"),
        legend.position  = "bottom"
      )
  )
}

## Vignette-local adapter: retain every required design identifier even when
## a pipeline release returns only scores or a subset of identifier columns.
attach_hrri_ids <- function(scores, id) {
  scores <- as.data.frame(scores)
  id <- as.data.frame(id)
  keys <- c("plot", "depth", "plant_id", "time")

  if (anyDuplicated(names(scores)) || anyDuplicated(names(id))) {
    stop("Score and identifier tables must have unique column names.")
  }
  missing_keys <- setdiff(keys, names(id))
  if (length(missing_keys)) {
    stop("sim$id is missing: ", paste(missing_keys, collapse = ", "))
  }
  if (!nrow(id) || nrow(scores) != nrow(id)) {
    stop("Expected one pipeline score row per sim$id row.")
  }
  if (!"RRI" %in% names(scores) || !is.numeric(scores[["RRI"]])) {
    stop("Pipeline scores must contain a numeric RRI column.")
  }

  make_key <- function(x) {
    if (anyNA(x[keys])) stop("Observation identifiers cannot be missing.")
    values <- lapply(x[keys], as.character)
    if (any(vapply(values, function(v) any(grepl("\034", v, fixed = TRUE)),
                   logical(1)))) {
      stop("Observation identifiers contain the reserved key separator.")
    }
    do.call(paste, c(values, list(sep = "\034")))
  }

  id_key <- make_key(id)
  if (anyDuplicated(id_key)) {
    stop("sim$id must uniquely identify plot x depth x plant x time rows.")
  }

  ## Put scores in the original simulator order whenever usable keys exist.
  if (all(keys %in% names(scores))) {
    score_key <- make_key(scores)
    if (anyDuplicated(score_key)) stop("Duplicate observation keys in scores.")
    idx <- match(id_key, score_key)
    alignment <- "observation keys"
  } else if ("row_id" %in% names(id) && "row_id" %in% names(scores)) {
    id_row <- as.character(id[["row_id"]])
    score_row <- as.character(scores[["row_id"]])
    if (anyNA(id_row) || anyNA(score_row) ||
        anyDuplicated(id_row) || anyDuplicated(score_row)) {
      stop("row_id must be unique and nonmissing for key-based alignment.")
    }
    idx <- match(id_row, score_row)
    alignment <- "row_id"
  } else {
    ## Scores-only releases must preserve their input row order. Row count
    ## equality alone cannot prove this; do not sort either table beforehand.
    idx <- seq_len(nrow(id))
    alignment <- "input row order (pipeline contract)"
  }
  if (anyNA(idx)) stop("Some simulator identifiers have no matching scores.")
  scores <- scores[idx, , drop = FALSE]

  shared <- intersect(c(keys, "row_id"), intersect(names(scores), names(id)))
  for (nm in shared) {
    if (!identical(as.character(scores[[nm]]), as.character(id[[nm]]))) {
      stop("Conflicting or misaligned identifier column: ", nm)
    }
  }

  ## Rebuild identifiers once; never create depth.x/depth.y pairs.
  out <- cbind(id[keys], scores[setdiff(names(scores), keys)])
  rownames(out) <- NULL
  attr(out, "id_alignment") <- alignment
  out
}

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

## ----simulate-----------------------------------------------------------------
library(HRRI)
packageVersion("HRRI")

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
nrow(sim$id)                   # one row per plot × depth × plant × time

## ----design-------------------------------------------------------------------
head(sim$id[, c("plot","depth","plant_id","time","cycle","phase","WFPS")])

## ----soil---------------------------------------------------------------------
head(sim$soil_data[, c("EAC","EDC","Cacc_EAC","Cacc_total","Cacc_fraction",
                        "FeIII_poor_crystalline_mmol_kg",
                        "FeII_mmol_kg","Eh","pH")])

## ----conservation-------------------------------------------------------------
## Maximum absolute error should be < 0.01 mmol kg-1
sim$conservation_checks

## ----plant--------------------------------------------------------------------
head(sim$plant_data[, c("SPAD","FvFm","ROL","ROS_load","aerenchyma")])

## ----genes--------------------------------------------------------------------
## 18 genes spanning Fe-cycling, denitrification, nitrification,
## methanogenesis, and sulfur cycling
colnames(sim$micro_gene_abundance)
summary(sim$micro_gene_abundance[, "mcrA"])   # methanogenesis gene

## ----cacc---------------------------------------------------------------------
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

## Mean accessible vs. total inventory
cat("Mean Cacc_raw:", mean(cap$cacc_raw, na.rm=TRUE), "mmol e- kg-1\n")
cat("Mean fraction :", mean(cap$cacc_fraction, na.rm=TRUE), "\n")
if ("ck_limited" %in% names(cap) && length(cap$ck_limited)) {
  cat("CK-limited rows:", sum(cap$ck_limited, na.rm=TRUE),
      "/", sum(!is.na(cap$ck_limited)), "classified rows\n")
} else {
  cat("CK-limited classification is not returned by this HRRI version.\n")
}

## ----tau_sweep----------------------------------------------------------------
tau_vals <- c(1, 6, 24, 72, 168, 720)   # 1 h to 30 d
cacc_tau <- sapply(tau_vals, function(tt) {
  r <- rri_accessible_capacity(sdf, res_spec, tau = tt, normalise = FALSE)
  mean(r$cacc_raw, na.rm = TRUE)
})
data.frame(tau_h = tau_vals, Cacc_mean = round(cacc_tau, 2))

## ----rri_pipeline-------------------------------------------------------------
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
summary(rri_scored$RRI)
head(rri_scored[, c("plot", "depth", "plant_id", "time",
                    "RRI", "Physio", "Soil", "Micro")])

## ----validation---------------------------------------------------------------
cat("HRRI: starting the 24-plot agreement example (100 bootstrap resamples).\n",
    file = stderr())
acc_sim <- simulate_redox_holobiont(
  n_plot=24, n_depth=2, n_plant=3, n_time=40,
  seed=4096, scenario="flood_drain", disturbance_strength=.70,
  n_cycles=2L, disturbance_center=NULL, disturbance_width=.08)
acc_res <- rri_pipeline_st(
  ROS_flux=acc_sim$plant_data, Eh_stability=acc_sim$Eh_stability,
  micro_data=log1p(acc_sim$micro_gene_abundance), id=acc_sim$id,
  time_col="time", group_cols=c("plot","depth","plant_id"), mode="snapshot",
  direction_anchor_phys="FvFm", direction_anchor_soil="Eh", direction_anchor_micro="mtrA")
acc <- rri_accuracy(
  score=acc_res$row_scores$RRI, target=acc_sim$latent_truth,
  cluster=acc_sim$id$plot, n_boot=100, n_perm=0, seed=20260913)
cat("HRRI: agreement calculations finished; preparing figures.\n", file = stderr())
acc

## ----validation_decomp--------------------------------------------------------
acc$decomposition[, c("component", "percent")]

## Exactness check: the residual is numerical noise, not a rounding allowance.
c(mse      = attr(acc$decomposition, "mse"),
  residual = attr(acc$decomposition, "residual"))

## ----validation_figure, fig.width=7.4, fig.height=6.5, out.width="100%", fig.alt="Agreement with the prescribed target"----
# Use the installed namespace, avoiding a stale function in the workspace.
accuracy_plot <- HRRI::plot_rri_accuracy
plot_formals <- names(formals(accuracy_plot))
plot_args <- list(acc = acc, score_label = "RRI",
                  target_label = "Prescribed target",
                  base_size = 9, show_clusters = FALSE)
# Older HRRI builds do not accept these presentation arguments.
if ("cluster_label" %in% plot_formals) plot_args$cluster_label <- "Plots"
if ("style" %in% plot_formals) plot_args$style <- "paper"
accuracy_figure <- do.call(accuracy_plot, plot_args)
# Without optional patchwork, HRRI returns a named list of plots.
if (inherits(accuracy_figure, c("ggplot", "patchwork"))) {
  print(accuracy_figure)
} else {
  for (panel in accuracy_figure) print(panel)
}

## ----recovery-----------------------------------------------------------------
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
metrics

## ----history_sensitivity------------------------------------------------------
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

## ----restore-theme, include=FALSE---------------------------------------------
## theme_set() changes state that persists for the rest of the session.
## Vignettes build in their own process so nothing outside is affected, but
## restoring is the same courtesy CRAN asks for with par() and options().
if (exists("old_theme")) ggplot2::theme_set(old_theme)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

