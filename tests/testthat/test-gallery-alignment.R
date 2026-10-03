test_that("gallery scores preserve forcing metadata and align after row reordering", {
  sim <- simulate_redox_holobiont(n_plot=2, n_depth=2, n_plant=3,
    n_time=40, seed=2026, scenario="flood_drain", disturbance_strength=.72,
    include_graph=TRUE, n_cycles=1L, disturbance_center=17,
    disturbance_width=5.5/sqrt(-2*log(.35))/40)
  expect_warning(res <- rri_pipeline_st(ROS_flux=sim$plant_data, Eh_stability=sim$Eh_stability,
    micro_data=log1p(sim$micro_gene_abundance), id=sim$id, time_col="time",
    group_cols=c("plot","depth","plant_id"), mode="snapshot",
    direction_anchor_phys="FvFm", direction_anchor_soil="Eh", direction_anchor_micro="mtrA"), "Excluding simulator-derived hidden columns")
  expect_identical(res$row_scores$WFPS, sim$id$WFPS)
  plot <- plot_rri_timeseries(sim,res,perturb_start=12,perturb_end=22,
    forcing_threshold=.35,time_label="Time (days)")
  expect_silent(ggplot2::ggplot_build(plot))
  reordered <- res
  reordered$row_scores <- res$row_scores[rev(seq_len(nrow(res$row_scores))),]
  again <- plot_rri_timeseries(sim,reordered,perturb_start=12,perturb_end=22,
    forcing_threshold=.35,time_label="Time (days)")
  expect_equal(plot$data, again$data)
  conflicting <- res
  conflicting$row_scores$WFPS <- conflicting$row_scores$WFPS + .001
  expect_error(plot_rri_timeseries(sim,conflicting), "Conflicting score metadata: WFPS")
})
