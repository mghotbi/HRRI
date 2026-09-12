# Agreement between a score and a reference target, respecting clustering

Quantifies how closely a score tracks a reference target, using
statistics appropriate to repeated observations of the same experimental
units. A pooled row-wise correlation over a longitudinal panel
overstates precision, because rows within a trajectory are not
independent. This function reports the naive row-level result alongside
the cluster-aware one, so the difference is visible rather than hidden.

## Usage

``` r
rri_accuracy(
  score,
  target,
  cluster = NULL,
  score_col = "RRI",
  n_boot = 1000,
  n_perm = 1000,
  conf = 0.95,
  seed = NULL
)

# S3 method for class 'rri_accuracy'
print(x, ...)
```

## Arguments

- score:

  Numeric vector of scores, or an `RRI` object. If an `RRI` object,
  `score_col` is taken from its `row_scores`.

- target:

  Numeric vector of reference values, the same length as `score`.

- cluster:

  Vector identifying the independent experimental unit for each
  observation, typically one trajectory. Rows sharing a value are
  treated as dependent. If `NULL`, every row is treated as independent
  and the function warns, because that assumption is rarely correct for
  time series.

- score_col:

  Column name used when `score` is an `RRI` object.

- n_boot:

  Number of cluster bootstrap resamples for confidence intervals. Set to
  0 to skip.

- n_perm:

  Number of cluster permutations for the null test. Set to 0 to skip.

- conf:

  Confidence level for intervals.

- seed:

  Optional integer seed for reproducible resampling.

- x:

  An `rri_accuracy` object.

- ...:

  Ignored.

## Value

An object of class `rri_accuracy`: a list with elements `agreement`,
`calibration`, `decomposition`, `dependence`, `null_test` and `notes`,
plus `data` (the complete cases actually used), `draws` (the resampled
and permuted statistics) and `conf`. The last three exist so that
[`plot_rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/plot_rri_accuracy.md)
can draw the sampling distributions without repeating the resampling.
See Details.

## Details

**Why correlation is not agreement.** Pearson's \\r\\ is invariant to
location and scale: a score equal to \\2 \times\\ the target plus a
constant correlates perfectly with it while agreeing with it nowhere.
Lin's concordance correlation coefficient

\$\$\rho_c = \frac{2 s\_{xy}}{s_x^2 + s_y^2 + (\bar{x} - \bar{y})^2}\$\$

penalises departure from the 1:1 line and is reported alongside \\r\\. A
large gap between the two means the score is well correlated but
miscalibrated.

**Why clustering matters.** With \\m\\ observations per unit and
intra-cluster correlation \\\rho\\, the design effect is \\1 + (m -
1)\rho\\ and the effective sample size is \\n / \mathrm{deff}\\. For a
40-point trajectory with \\\rho = 0.3\\ this is roughly a twelvefold
reduction, so a confidence interval computed from the row count is far
too narrow. Intervals here come from a cluster bootstrap that resamples
whole trajectories with replacement, which preserves the within-unit
dependence.

**Error decomposition.** Mean squared error is split following Kobayashi
and Salam (2000) into squared bias, a difference in variability, and a
lack of correlation:

\$\$\mathrm{MSE} = (\bar{x} - \bar{y})^2 + (s_x - s_y)^2 + 2 s_x s_y
(1 - r)\$\$

These say different things. Large squared bias means a systematic
offset, correctable by recentring. Large variance mismatch means the
score is too flat or too volatile. Large lack of correlation means the
score does not track the target's pattern, and no rescaling will fix it.

**What the permutation null asks.** Whole trajectories are exchanged, so
each unit keeps its own temporal shape and only the pairing between
score and target is broken. This is deliberately the harder null.
Permuting individual rows would destroy the shared event-driven shape
that every trajectory has, making almost any score look significant;
exchanging blocks retains that shape and asks whether the score tracks
*this* unit's target beyond what the common disturbance imposes on all
of them. A small p-value under this null is therefore informative, and a
large one is not evidence that the score is uninformative about the
disturbance itself.

**What this does not establish.** If the score and the target come from
the same generative model, as they do for `latent_truth` from
[`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md),
high agreement is internal consistency and nothing more. It is not
predictive accuracy, not out-of-sample error, and not evidence of
ecological validity. Those require a target constructed independently of
the score, and replication at the level of independent experimental
units.

## Methods (by generic)

- `print(rri_accuracy)`: Compact console summary of an accuracy
  assessment.

## References

Lin, L.I. (1989) A concordance correlation coefficient to evaluate
reproducibility. *Biometrics*, **45**, 255–268.

Kobayashi, K. & Salam, M.U. (2000) Comparing simulated and measured
values using mean squared deviation and its components. *Agronomy
Journal*, **92**, 345–352.

## See also

[`plot_rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/plot_rri_accuracy.md)
for the four-panel diagnostic figure;
[`benchmark_hrri()`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
for repeated-seed benchmarking;
[`rri_domain_influence()`](https://mghotbi.github.io/HRRI/reference/rri_domain_influence.md)
for which domain drives the score.

## Examples

``` r
## Synthetic panel: 8 trajectories, 20 time points each.
## The score is deliberately miscalibrated so the r-vs-CCC gap is visible.
set.seed(1)
k <- 8; m <- 20
unit   <- rnorm(k, 0, 0.30)
target <- unlist(lapply(unit, function(u) u + 0.5 + rnorm(m, 0, 0.05)))
score  <- 0.75 * target + 0.10 + rnorm(k * m, 0, 0.06)
traj   <- rep(seq_len(k), each = m)

acc <- rri_accuracy(score, target, cluster = traj,
                    n_boot = 200, n_perm = 200, seed = 1)
acc
#> Agreement with reference target
#> -------------------------------------------------------------- 
#>  statistic row_level cluster_mean_level ci_lower ci_upper
#>  pearson_r    0.9482             0.9974   0.8676   0.9665
#>   lins_ccc    0.9135             0.9464   0.7711   0.9281
#>       rmse    0.0903             0.0683   0.0674   0.1187
#>        mae    0.0709             0.0548   0.0542   0.0971
#>       bias   -0.0297            -0.0297  -0.0802   0.0067
#>         r2    0.8583             0.9157   0.5766   0.8843
#> 
#> Dependence structure
#>   160 observations in 8 clusters (mean size 20.0)
#>   ICC 0.871, design effect 17.6, effective n 9
#> 
#> Error decomposition (percent of MSE)
#>   squared_bias          10.9%
#>   variance_mismatch     31.5%
#>   lack_of_correlation   57.6%
#> 
#> Calibration: target = -0.073 + 1.202 x score  (ideal 0 and 1)
#> Cluster permutation test: p = 0.0050 (200 permutations)
#> 
#> Notes
#>   Rows are strongly clustered (ICC 0.87, design effect 17.6). The 160
#>   observations carry roughly the information of 9 independent ones; quote
#>   the cluster bootstrap interval, not one based on n = 160. 
#>   Ignoring clustering would give a 95% interval for r of width 0.033;
#>   resampling whole trajectories gives width 0.099, 3.0 times wider. Report
#>   the latter. 
#>   Calibration slope is 1.20 rather than 1: the score compresses or
#>   exaggerates the target's range. 
#>   Error is dominated by lack_of_correlation (58% of MSE). 
#>   If score and target derive from the same generator, this is internal
#>   consistency, not validation. 

## Correlation is high, concordance is not: the score compresses the target.
acc$agreement
#>   statistic   row_level cluster_mean_level    ci_lower    ci_upper
#> 1 pearson_r  0.94824740         0.99743125  0.86759471 0.966528380
#> 2  lins_ccc  0.91348052         0.94642199  0.77105411 0.928109217
#> 3      rmse  0.09025121         0.06834325  0.06738456 0.118737771
#> 4       mae  0.07085095         0.05483952  0.05424316 0.097145734
#> 5      bias -0.02974388        -0.02974388 -0.08018564 0.006706731
#> 6        r2  0.85833006         0.91574196  0.57661875 0.884330780
acc$calibration
#>   intercept    slope intercept_lower intercept_upper slope_lower slope_upper
#> 1 -0.073473 1.202293      -0.1053678      0.02958112    1.009207    1.270628

## The row count is not the sample size.
acc$dependence
#>   n_observations n_clusters mean_cluster_size       icc design_effect
#> 1            160          8                20 0.8714487      17.55753
#>   effective_n r_ci_naive_width r_ci_cluster_width
#> 1    9.112901       0.03280238         0.09893367

# \donttest{
## Against the simulator's own prescribed target.
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 3, n_time = 30,
  p_micro = 20, seed = 2026
)
res <- suppressWarnings(rri_pipeline_st(
  ROS_flux = sim$ROS_flux, Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data, id = sim$id
))
scored <- attach_hrri_ids(res$row_scores, sim$id)
tj <- interaction(scored$plot, scored$depth, scored$plant_id, drop = TRUE)
rri_accuracy(scored$RRI, sim$latent_truth, cluster = tj,
             n_boot = 500, n_perm = 500, seed = 1)
#> Agreement with reference target
#> -------------------------------------------------------------- 
#>  statistic row_level cluster_mean_level  ci_lower ci_upper
#>  pearson_r   -0.7455            -0.7482   -0.8424  -0.6131
#>   lins_ccc   -0.0648            -0.0317   -0.1238  -0.0219
#>       rmse    0.3768             0.3559    0.3110   0.4363
#>        mae    0.3445             0.3191    0.2759   0.4169
#>       bias   -0.3191            -0.3191   -0.4088  -0.2341
#>         r2 -127.6569          -384.5064 -351.6890 -64.2928
#> 
#> Dependence structure
#>   360 observations in 12 clusters (mean size 30.0)
#>   ICC 0.685, design effect 20.9, effective n 17
#> 
#> Error decomposition (percent of MSE)
#>   squared_bias          71.7%
#>   variance_mismatch     14.0%
#>   lack_of_correlation   14.2%
#> 
#> Calibration: target = 0.671 + -0.142 x score  (ideal 0 and 1)
#> Cluster permutation test: p = 0.0020 (500 permutations)
#> 
#> Notes
#>   Rows are strongly clustered (ICC 0.68, design effect 20.9). The 360
#>   observations carry roughly the information of 17 independent ones; quote
#>   the cluster bootstrap interval, not one based on n = 360. 
#>   Ignoring clustering would give a 95% interval for r of width 0.095;
#>   resampling whole trajectories gives width 0.229, 2.4 times wider. Report
#>   the latter. 
#>   Calibration slope is -0.14 rather than 1: the score compresses or
#>   exaggerates the target's range. 
#>   R2 is negative: as an absolute predictor the score does worse than the
#>   target's own mean. It may still rank correctly; check pearson_r. 
#>   Error is dominated by squared_bias (72% of MSE). 
#>   If score and target derive from the same generator, this is internal
#>   consistency, not validation. 
# }
```
