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
knitr::opts_chunk$set(dev = "png", collapse=TRUE,comment="#>",fig.width=7.4,
  fig.height=4.8,dpi=200,out.width="100%",message=FALSE,warning=FALSE)
library(HRRI)
# Fast documentation build; export script defaults to all 2000 paper resamples.
paper_n_boot <- 100L

if (!exists("paper_n_boot", inherits = FALSE)) paper_n_boot <- 2000L
library(HRRI)
gal <- simulate_redox_holobiont(n_plot=2,n_depth=2,n_plant=3,n_time=40,
  seed=2026,scenario="flood_drain",disturbance_strength=.72,include_graph=TRUE,
  n_cycles=1L,disturbance_center=17,
  disturbance_width=5.5/sqrt(-2*log(.35))/40)
fit <- function(sim) rri_pipeline_st(ROS_flux=sim$plant_data,Eh_stability=sim$Eh_stability,
  micro_data=log1p(sim$micro_gene_abundance),id=sim$id,time_col="time",
  group_cols=c("plot","depth","plant_id"),mode="snapshot",
  direction_anchor_phys="FvFm",direction_anchor_soil="Eh",direction_anchor_micro="mtrA")
gal_res <- fit(gal)
event_start <- min(gal$id$time[gal$id$phase=="disturbance"])
event_end <- max(gal$id$time[gal$id$phase=="disturbance"])
gal_scores <- attach_hrri_ids(gal_res$row_scores,gal$id)
gal_scores$WFPS <- gal$forcing$WFPS
gal_rec <- rri_recovery_metrics(gal_scores,time_col="time",
  group_cols=c("plot","depth","plant_id"),perturb_start=event_start,
  perturb_end=event_end,forcing_col="WFPS")
props <- rri_property_scores(gal_res,rec=gal_rec,soil_df=gal$soil_data)
acc_sim <- simulate_redox_holobiont(n_plot=24,n_depth=2,n_plant=3,n_time=40,
  seed=4096,scenario="flood_drain",disturbance_strength=.70,n_cycles=2L,
  disturbance_center=NULL,disturbance_width=.08)
acc_res <- fit(acc_sim)
acc <- rri_accuracy(acc_res$row_scores$RRI,acc_sim$latent_truth,
  cluster=acc_sim$id$plot,n_boot=paper_n_boot,n_perm=0,seed=20260913)
paper_figures <- list(
  Fig1_framework=plot_rri_framework(),
  Fig2_identifiability=plot_rri_identifiability(),
  Fig3_trajectories=plot_rri_timeseries(gal,gal_res,
    perturb_start=event_start,perturb_end=event_end,forcing_threshold=.35,time_label="Time (days)"),
  Fig4_recovery_diagnostics=plot_rri_recovery_diagnostics(gal_res,gal$id,gal_rec,
    perturb_start=event_start,perturb_end=event_end),
  Fig5_operational_profile=plot_rri_properties(props,rec=gal_rec),
  Fig6_internal_agreement=plot_rri_accuracy(acc,base_size=9,
    score_label="Observation-derived score",target_label="Prescribed target",
    show_clusters=FALSE,cluster_label="Plots")
)

show_figure <- function(p) {
  if(inherits(p,c("ggplot","patchwork"))) print(p) else for(q in p) print(q)
}

## ----framework, fig.alt="Figure 1: complementary routes"----------------------
show_figure(paper_figures$Fig1_framework)

## ----identifiability, fig.height=6.4, fig.alt="Figure 2: identifiability, sensitivity and accounting"----
show_figure(paper_figures$Fig2_identifiability)

## ----trajectories, fig.height=6, fig.alt="Figure 3: forcing and responses"----
show_figure(paper_figures$Fig3_trajectories)

## ----recovery, fig.width=8, fig.height=4.6, fig.alt="Figure 4: recovery diagnostics"----
show_figure(paper_figures$Fig4_recovery_diagnostics)

## ----properties, fig.height=4.5, fig.alt="Figure 5: operational properties"----
show_figure(paper_figures$Fig5_operational_profile)

## ----agreement, fig.height=6.5, fig.alt="Figure 6: conditional agreement"-----
show_figure(paper_figures$Fig6_internal_agreement)

## ----export, eval=FALSE-------------------------------------------------------
# figure_dir <- file.path(getwd(), "HRRI_paper_figures")
# paper_n_boot <- 2000L
# source(system.file("examples", "export_paper_figures.R", package="HRRI"), local=TRUE)

