#' Benchmark diagnostic agreement with a simulator-defined target
#' @description Runs independent seeds and reports descriptive agreement, not
#' held-out prediction, empirical validation or parameter identification.
#' @param domains Nonempty subset of soil, plant and micro.
#' @param missing Fraction in `[0, 1)` of uniformly sampled cells removed (MCAR).
#' This does not implement informative or MNAR missingness.
#' @param noise Nonnegative Gaussian noise SD on each column's original scale.
#' A common SD has different relative effects on differently scaled variables.
#' @param n Positive integer number of independent simulations.
#' @param seed_start First integer seed. The caller's RNG state is restored.
#' @param sim_args Named simulator arguments, excluding seed.
#' @param verbose Print progress messages.
#' @param pipeline_args Named additional rri_pipeline arguments, excluding
#' dat, soil, plant, micro and id; use this to justify orientations and weights.
#' @return An hrri_benchmark list with summary, seed_metrics, row_data,
#' failures (seed and error), and settings including package/R versions.
#' Summary RMSE and Bias compare the score directly with the chosen target.
#' r_truth and rank_truth are pooled descriptive correlations; spread_association
#' compares within-seed score and target SDs. None is interval coverage or
#' predictive uncertainty. n_rows counts finite matched rows across all seeds.
#' @details Soil and plant synthetic observations and log1p gene abundances
#' feed the pipeline. Latent architecture and microbial activity states are not
#' scoring inputs. The latent target is still a prescribed simulator composite;
#' it is not an independently measured recovery outcome. Replicated rows within
#' seeds are dependent. Compare methods using held-out seeds and independent
#' process outcomes in a separate validation design.
#' @seealso rri_pipeline, simulate_redox_holobiont, plot_hrri_benchmark
#' @examples
#' \dontrun{
#' b <- benchmark_hrri(domains = "soil", n = 2, missing = 0.1)
#' print(b)
#' }
#' @importFrom stats rnorm sd
#' @importFrom utils modifyList packageVersion
#' @export
benchmark_hrri <- function(domains = c("soil", "plant", "micro"),
                            missing = 0, noise = 0, n = 50L,
                            seed_start = 1L, sim_args = NULL, verbose = TRUE,
                            pipeline_args = list()) {
  # Assisted-by: OpenAI Codex. Seed-local corruption, explicit failures and metrics.
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is.character(domains) || !length(domains) ||
      anyNA(domains) || any(!domains %in% c("soil", "plant", "micro")))
    stop("domains must be a nonempty subset of soil, plant and micro.")
  domains <- unique(domains)
  if (!scalar(missing) || missing < 0 || missing >= 1) stop("missing must be in `[0, 1)`.")
  if (!scalar(noise) || noise < 0) stop("noise must be nonnegative and finite.")
  if (!scalar(n) || n < 1 || n != floor(n)) stop("n must be a positive integer.")
  if (!scalar(seed_start) || seed_start < 0 || seed_start != floor(seed_start) ||
      seed_start + n - 1 > .Machine$integer.max)
    stop("Seeds must be nonnegative integers within R's integer range.")
  named_list <- function(x, what) {
    if (!is.list(x) || (length(x) && (is.null(names(x)) ||
        any(!nzchar(names(x))) || anyDuplicated(names(x)))))
      stop(what, " must be a uniquely named list.")
  }
  if (is.null(sim_args)) sim_args <- list()
  named_list(sim_args, "sim_args"); named_list(pipeline_args, "pipeline_args")
  if ("seed" %in% names(sim_args)) stop("Use seed_start, not sim_args$seed.")
  if (any(names(pipeline_args) %in% c("dat", "soil", "plant", "micro", "id")))
    stop("pipeline_args must not replace the benchmark observation inputs.")
  sim_args <- utils::modifyList(list(n_plot = 1, n_depth = 1, n_plant = 2,
                                    n_time = 30, p_micro = 20), sim_args)
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  corrupt <- function(x) {
    x <- .rri_observed_block(x)
    if (is.null(x)) return(NULL)
    m <- as.matrix(x)
    if (noise > 0) m <- m + matrix(stats::rnorm(length(m), sd = noise), nrow(m))
    if (missing > 0) m[sample.int(length(m), floor(length(m) * missing))] <- NA_real_
    as.data.frame(m)
  }
  seeds <- seed_start + seq_len(n) - 1L
  results <- vector("list", n); errors <- vector("list", n)
  for (i in seq_along(seeds)) {
    s <- seeds[i]
    results[i] <- list(tryCatch({
      sim <- do.call(simulate_redox_holobiont, c(sim_args, list(seed = s)))
      truth <- sim$latent_truth
      if (!is.numeric(truth) || length(truth) != nrow(sim$id))
        stop("Simulator truth must be numeric and aligned to id.")
      set.seed(s)
      plant <- if ("plant" %in% domains) corrupt(sim$ROS_flux) else NULL
      soil <- if ("soil" %in% domains) corrupt(sim$Eh_stability) else NULL
      micro <- if ("micro" %in% domains) {
        if (is.null(sim$micro_gene_abundance))
          stop("Simulator did not supply observable micro_gene_abundance.")
        corrupt(log1p(.rri_numeric_df(sim$micro_gene_abundance)))
      } else NULL
      res <- do.call(rri_pipeline, c(list(soil = soil, plant = plant,
                            micro = micro, id = sim$id), pipeline_args))
      keys <- intersect(c("plot", "depth", "plant_id", "time"), names(sim$id))
      rs <- .rri_align_scores(res$row_scores, sim$id, keys)
      if (!"RRI" %in% names(rs) || !is.numeric(rs$RRI)) stop("Pipeline returned no numeric RRI.")
      ok <- is.finite(rs$RRI) & is.finite(truth)
      if (sum(ok) < 3L) stop("Fewer than three finite matched score/target pairs.")
      rs$seed <- s
      rs$truth <- truth
      rs$residual <- rs$RRI - truth
      rs
    }, error = function(e) {
      errors[i] <<- list(data.frame(seed = s, error = conditionMessage(e)))
      NULL
    }))
    if (isTRUE(verbose)) message("benchmark_hrri: ", i, "/", n, " seeds attempted")
  }
  failures <- do.call(rbind, Filter(Negate(is.null), errors))
  if (is.null(failures)) failures <- data.frame(seed = numeric(), error = character())
  row_data <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(row_data)) stop("All seeds failed: ",
    paste(paste(failures$seed, failures$error, sep = ": "), collapse = "; "))
  if (nrow(failures)) warning(nrow(failures), " seeds failed; inspect $failures.", call. = FALSE)
  rownames(row_data) <- NULL
  metric <- function(d) {
    d <- d[is.finite(d$RRI) & is.finite(d$truth), , drop = FALSE]
    data.frame(RMSE = sqrt(mean(d$residual^2)), Bias = mean(d$residual),
      r_truth = .rri_cor(d$RRI, d$truth),
      rank_truth = .rri_cor(d$RRI, d$truth, method = "spearman"),
      sd_RRI = stats::sd(d$RRI), sd_truth = stats::sd(d$truth), n_rows = nrow(d))
  }
  seed_metrics <- do.call(rbind, lapply(unique(row_data$seed), function(s)
    cbind(seed = s, metric(row_data[row_data$seed == s, , drop = FALSE]))))
  summary <- metric(row_data)[, c("RMSE", "Bias", "r_truth", "rank_truth", "n_rows")]
  summary$spread_association <- .rri_cor(seed_metrics$sd_RRI, seed_metrics$sd_truth,
                                        method = "spearman")
  summary$mean_within_seed_sd <- .rri_mean(seed_metrics$sd_RRI)
  summary$n_requested <- n
  summary$n_seeds <- nrow(seed_metrics)
  summary$n_failed <- nrow(failures)
  structure(list(summary = summary, seed_metrics = seed_metrics,
    row_data = row_data, failures = failures,
    settings = list(domains = domains, missing = missing, missingness = "MCAR",
      noise = noise, n = n, seed_start = seed_start, sim_args = sim_args,
      pipeline_args = pipeline_args, R_version = R.version.string,
      HRRI_version = tryCatch(as.character(utils::packageVersion("HRRI")),
                              error = function(e) NA_character_),
      target = "Simulator composite: descriptive agreement only")),
    class = "hrri_benchmark")
}

#' @rdname benchmark_hrri
#' @param x An hrri_benchmark object.
#' @param ... Unused method arguments.
#' @export
print.hrri_benchmark <- function(x, ...) {
  cat("<hrri_benchmark> descriptive agreement, not predictive validation\n")
  print(x$summary, row.names = FALSE)
  if (nrow(x$failures)) cat("Inspect $failures for unsuccessful seeds.\n")
  invisible(x)
}

#' Plot descriptive benchmark agreement
#' @param bm An hrri_benchmark object.
#' @param print Whether to display the plot.
#' @param colour Scatter colour.
#' @return Invisibly, a two-panel patchwork object. Requires patchwork.
#' @importFrom ggplot2 ggplot aes geom_abline geom_point geom_histogram labs coord_fixed theme_classic
#' @examples
#' \dontrun{
#'   b <- benchmark_hrri(domains = "soil", n = 2, missing = 0.1)
#'   plot_hrri_benchmark(b)
#' }
#' @export
plot_hrri_benchmark <- function(bm, print = TRUE, colour = "#0072B2") {
  if (!inherits(bm, "hrri_benchmark")) stop("bm must be an hrri_benchmark.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Plotting requires patchwork.")
  d <- bm$row_data[is.finite(bm$row_data$RRI) & is.finite(bm$row_data$truth), , drop = FALSE]
  p1 <- ggplot2::ggplot(d, ggplot2::aes(x = truth, y = RRI)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8, colour = colour) +
    ggplot2::labs(x = "Simulator-defined target", y = "Observation-derived score",
      subtitle = sprintf("Descriptive r = %.3f; RMSE = %.3f",
                         bm$summary$r_truth, bm$summary$RMSE)) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_classic(base_size = 10)
  p2 <- ggplot2::ggplot(bm$seed_metrics, ggplot2::aes(x = RMSE)) +
    ggplot2::geom_histogram(bins = 15, fill = colour, colour = "white") +
    ggplot2::labs(x = "Within-seed RMSE", y = "Number of seeds") +
    ggplot2::theme_classic(base_size = 10)
  p <- patchwork::wrap_plots(p1, p2, ncol = 2)
  if (isTRUE(print)) base::print(p)
  invisible(p)
}
