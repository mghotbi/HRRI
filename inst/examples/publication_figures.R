
library(HRRI)
library(ggplot2)
figure_dir <- file.path(tempdir(), "HRRI_publication_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
sim <- simulate_redox_holobiont(n_plot = 6, n_depth = 2, n_plant = 1,
                               n_time = 40, p_micro = 20, seed = 2026)
res <- rri_pipeline_st(ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data, id = sim$id, reducer = "per_domain", scaling = "pnorm")
# Illustrative analysis window. In real studies specify endpoints from the
# protocol or an explicit forcing rule, and examine sensitivity to that choice.
event_start <- 12
event_end <- 22
rec <- rri_recovery_metrics(res, id = sim$id, time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  perturb_start = event_start, perturb_end = event_end)
props <- rri_property_scores(res, rec = rec, soil_df = sim$soil_data)
figures <- list(
  profile = plot_rri_properties(props, base_size = 11),
  recovery = plot_rri_recovery_landscape(rec, base_size = 10),
  recovery_map = plot_rri_recovery_map(res, sim$id,
    perturb_start = event_start, perturb_end = event_end),
  timeseries = plot_rri_timeseries(sim, res,
    perturb_start = event_start, perturb_end = event_end)
)
# Separate generic agreement example: cluster is the independent sampling unit.
# Replace BOTH paired vectors and the cluster key with your aligned data.
set.seed(103)
cluster <- rep(seq_len(12), each = 20)
target <- rep(runif(12, .2, .8), each = 20) + rnorm(240, 0, .04)
score <- .1 + .7 * target + rnorm(240, 0, .04)
acc <- rri_accuracy(score, target, cluster = cluster,
                    n_boot = 500, n_perm = 0, seed = 104)
figures$agreement <- plot_rri_accuracy(acc, show_clusters = FALSE)
for (nm in names(figures)) {
  w <- if (nm == "agreement") 11 else 8
  h <- if (nm == "agreement") 8 else 5.4
  # PDF keeps vector geometry and text; PNG is for previews.
  ggsave(file.path(figure_dir, paste0(nm, ".pdf")), figures[[nm]],
         width = w, height = h, units = "in", device = "pdf")
  ggsave(file.path(figure_dir, paste0(nm, ".png")), figures[[nm]],
         width = w, height = h, units = "in", dpi = 300)
}
write.csv(rec, file.path(figure_dir, "recovery_values_and_fit_status.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(figure_dir, "sessionInfo.txt"))
message("Figures and underlying recovery table: ", figure_dir)
