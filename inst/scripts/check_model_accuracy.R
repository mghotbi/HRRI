## ---------------------------------------------------------------------------
## check_model_accuracy.R
##
## Runs the accuracy assessment against your CURRENTLY INSTALLED HRRI and
## writes a publication figure. It sources the two new files directly, so you
## can see the numbers before rebuilding the package.
##
##   setwd("<this folder>")
##   source("check_model_accuracy.R")
##
## After installing HRRI 0.99.4 the two source() lines below are unnecessary:
## rri_accuracy() and plot_rri_accuracy() are exported.
##
## Replaces the claim "r = 0.470 on 480 observations" with one a reviewer
## cannot take apart: those 480 rows are 12 trajectories x 40 time points.
## ---------------------------------------------------------------------------

library(HRRI)
library(ggplot2)

if (!exists("rri_accuracy", mode = "function"))      source("rri_accuracy.R")
if (!exists("plot_rri_accuracy", mode = "function")) source("plot_rri_accuracy.R")

out_dir <- "accuracy_output"
dir.create(out_dir, showWarnings = FALSE)

set.seed(2026)

## ---- 1. the same experiment as the gallery vignette -----------------------
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 3, n_time = 40,
  p_micro = 25, seed = 2026,
  scenario = "flood_drain", n_cycles = 1,
  disturbance_strength = 0.70, history_strength = 0.55
)

res <- suppressWarnings(rri_pipeline_st(
  ROS_flux     = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data   = sim$micro_data,
  id           = sim$id,
  time_col     = "time",
  group_cols   = c("plot", "depth"),
  reducer      = "per_domain",
  scaling      = "pnorm"
))

scored <- attach_hrri_ids(res$row_scores, sim$id)
score  <- scored$RRI
target <- sim$latent_truth

## One independent experimental unit = one plot x depth x plant trajectory.
traj <- interaction(scored$plot, scored$depth, scored$plant_id, drop = TRUE)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("Design: ", length(score), " rows = ", nlevels(traj),
    " trajectories x ", length(score) / nlevels(traj), " time points\n", sep = "")
cat(strrep("=", 72), "\n\n", sep = "")


## ---- 2. what the naive summary says ---------------------------------------
cat("--- The claim as currently stated ------------------------------------\n")
ok <- is.finite(score) & is.finite(target)
ct <- cor.test(score[ok], target[ok])
cat(sprintf("r = %.3f, 95%% CI [%.3f, %.3f], n = %d\n",
            ct$estimate, ct$conf.int[1], ct$conf.int[2], sum(ok)))
cat("That interval assumes 480 independent observations. There are not 480.\n\n")


## ---- 3. the defensible version --------------------------------------------
acc <- rri_accuracy(
  score = score, target = target, cluster = traj,
  n_boot = 2000, n_perm = 2000, seed = 2026
)
print(acc)


## ---- 4. the figure ---------------------------------------------------------
fig <- plot_rri_accuracy(
  acc,
  score_label  = "RRI",
  target_label = "Prescribed target",
  base_size    = 11
)

ggsave(file.path(out_dir, "Fig_accuracy.pdf"), fig,
       width = 9.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig_accuracy.png"), fig,
       width = 9.5, height = 7.5, dpi = 400)
cat(sprintf("\nFigure written to %s/Fig_accuracy.{pdf,png}\n", out_dir))


## ---- 5. sensitivity: does the conclusion survive re-simulation? -----------
## Repeating with fresh seeds separates "this estimator works" from
## "this seed worked".
cat("\n--- Across 10 independent simulations ---------------------------------\n")
reps <- do.call(rbind, lapply(1:10, function(s) {
  ss <- simulate_redox_holobiont(
    n_plot = 2, n_depth = 2, n_plant = 3, n_time = 40,
    p_micro = 25, seed = 1000 + s,
    scenario = "flood_drain", n_cycles = 1,
    disturbance_strength = 0.70, history_strength = 0.55
  )
  rr <- suppressWarnings(rri_pipeline_st(
    ROS_flux = ss$ROS_flux, Eh_stability = ss$Eh_stability,
    micro_data = ss$micro_data, id = ss$id,
    reducer = "per_domain", scaling = "pnorm"
  ))
  sc <- attach_hrri_ids(rr$row_scores, ss$id)
  tj <- interaction(sc$plot, sc$depth, sc$plant_id, drop = TRUE)
  a  <- rri_accuracy(sc$RRI, ss$latent_truth, cluster = tj,
                     n_boot = 0, n_perm = 0)
  v  <- setNames(a$agreement$row_level, a$agreement$statistic)
  data.frame(seed = 1000 + s, r = v[["pearson_r"]], ccc = v[["lins_ccc"]],
             rmse = v[["rmse"]], bias = v[["bias"]],
             icc = a$dependence$icc, slope = a$calibration$slope)
}))
print(round(reps, 3), row.names = FALSE)

cat(sprintf("\nAcross seeds: r = %.3f (sd %.3f), CCC = %.3f (sd %.3f)\n",
            mean(reps$r), sd(reps$r), mean(reps$ccc), sd(reps$ccc)))
cat("Report the mean and spread across seeds, not a single run.\n")

## Seed-to-seed spread, as a small companion figure.
seed_fig <- ggplot(
  data.frame(stat = rep(c("Pearson r", "Lin's CCC"), each = nrow(reps)),
             value = c(reps$r, reps$ccc)),
  aes(x = stat, y = value)) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = "#eaf0ef",
               colour = "#2f6b6b") +
  geom_jitter(width = 0.08, size = 1.8, alpha = 0.75, colour = "#8c4a2f") +
  labs(title = "Stability across 10 independent simulations",
       subtitle = "Each point is one seed; report the spread, not one run",
       x = NULL, y = NULL) +
  theme_ems(base_size = 11)

ggsave(file.path(out_dir, "Fig_accuracy_seeds.png"), seed_fig,
       width = 5.2, height = 4.2, dpi = 400)
cat(sprintf("Seed-stability figure written to %s/Fig_accuracy_seeds.png\n",
            out_dir))


## ---- 6. which domain carries the score ------------------------------------
cat("\n--- Realised domain influence ------------------------------------------\n")
infl <- rri_domain_influence(res)
print(infl$influence, row.names = FALSE)
cat("\n")
for (nt in infl$notes) {
  cat(paste(strwrap(nt, 76, prefix = "  "), collapse = "\n"), "\n")
}


cat("\n", strrep("=", 72), "\n", sep = "")
cat("Reminder: latent_truth and RRI share a generator. Everything above is\n")
cat("internal consistency of the estimator, not empirical validation.\n")
cat(strrep("=", 72), "\n", sep = "")
