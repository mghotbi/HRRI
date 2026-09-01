#' Simulate illustrative soil-plant-microbe redox trajectories
#'
#' @description Synthetic daily trajectories with explicit Fe/Mn redistribution,
#' plant indicators, gene abundance and transcript/count observation models.
#' Only Fe and Mn inventories have closed-balance checks. C, N, S and oxygen
#' budgets are not fully balanced. Model parameters are illustrative, not fitted.
#' Gene abundance indicates potential; transcript counts are observations, not flux.
#' @param n_plot Positive integer. Number of plots (spatial replicates).
#' @param n_depth Positive integer. Number of depth strata per plot.
#' @param n_plant Positive integer. Number of plants per plot-depth unit.
#'   Must be \eqn{\geq 1}. Plants contribute to ROL and ROS signals.
#' @param n_time Positive integer \eqn{\geq 4}. Number of daily time steps.
#' @param p_micro Positive integer. Number of ASV-like microbial taxonomic
#'   features generated alongside functional gene data.
#' @param seed Integer or \code{NULL}. Random seed passed to
#'   \code{\link[base]{set.seed}} before simulation; \code{NULL} means
#'   no seeding (non-reproducible). All manuscript figures use explicit seeds.
#' @param scenario Character; one of \code{"flood_drain"} (default) or
#'   \code{"drought_rewet"}. Determines the shape of the hydrological
#'   forcing function and the sign of the dominant redox transition.
#' @param n_cycles Positive integer. Number of forcing pulses.
#' @param disturbance_strength Numeric in \eqn{[0, 1]}. Event severity.
#'   Controls peak WFPS, anaerobic volume fraction, and the amplitude of
#'   soil redox transitions.
#' @param disturbance_center Numeric or \code{NULL}. Time step of the first
#'   disturbance event centre. Defaults to evenly spaced centres from 0.22 to 0.78 of n_time.
#' @param disturbance_width Numeric in \eqn{(0, 1)}. Width of each Gaussian
#'   forcing pulse as a fraction of \code{n_time}.
#' @param seasonal_amp Numeric \eqn{\geq 0}. Amplitude of additive seasonal
#'   forcing overlaid on the hydrological disturbance signal.
#' @param seasonal_phase Numeric. Phase offset (radians) of the seasonal forcing.
#' @param history_strength Numeric in \eqn{[0, 1]}. Scales the synthetic
#'   memory state at initialization and during disturbance. That state affects
#'   Fe crystallisation, accessibility and generated microbial descriptors; it
#'   is not a measured fraction of community carry-over.
#' @param rescue Character; one of \code{"none"} (default), \code{"capacity"},
#'   \code{"connectivity"}, or \code{"kinetics"}. Simulates a targeted
#'   scenario modification (capacity at initialization; alpha/k throughout):
#'   \describe{
#'     \item{\code{"capacity"}}{Fe(III) inventory replenishment.}
#'     \item{\code{"connectivity"}}{Increases alpha_accept and alpha_donate in selected Fe/Mn/N/S/C rate expressions and calculated accessible capacity. Effects on any recovery outcome must be evaluated, not assumed beneficial.}
#'     \item{\code{"kinetics"}}{Increases both calculated accessible-capacity rates. k_accept additionally gates crystalline-Fe reduction; k_donate has no direct process-rate gate. This intervention is not a general exchange-kinetics model.}
#'   }
#' @param event_tau_h Positive numeric. Disturbance timescale \eqn{\tau} (h)
#'   used internally by the accessible-capacity calculation; passed to
#'   \code{\link{rri_accessible_capacity}}.
#' @param sequencing_depth Positive numeric. Mean library size for taxonomic
#'   count data (\code{micro_data} block), modelled as a negative-binomial
#'   process.
#' @param metat_depth Positive numeric. Mean library size for metatranscript
#'   counts (\code{micro_metat_counts}), modelled separately from taxonomic
#'   counts with a higher biological variance.
#' @param decoupling Numeric in \eqn{[0, 1]}. Cross-domain stochastic
#'   decoupling parameter for selected noise terms; other stochastic terms remain.
#' @param zero_inflation Numeric in \eqn{[0, 1]}. Structural-zero probability
#'   for taxonomic count features (simulates taxa absent from some samples).
#' @param MNAR_strength Numeric in \eqn{[0, 1]}. Maximum missing-not-at-random
#'   (MNAR) probability for Eh values under strongly reducing conditions.
#' @param Eh_dropout_threshold Numeric. Eh (mV) below which the MNAR dropout
#'   probability begins to rise; is an artificial missingness design, not a platinum-electrode detection limit.
#' @param micro_mean,micro_slope,micro_lambda_min,micro_lambda_max
#'   Backward-compatible parameters controlling mean and slope of the
#'   log-linear model for taxonomic count intensity. See legacy documentation.
#' @param stochastic_reassembly Logical. If \code{TRUE} (default), adds
#'   additional stochastic variation; this is not an explicit succession model.
#' @param include_graph Logical. If \code{TRUE} and \pkg{igraph} is installed,
#'   returns an \code{igraph} random graph object independent of the generated community in
#'   \code{$graph}.
#' @param depth_labels Character vector of length \code{n_depth} or
#'   \code{NULL}. Custom labels for depth strata. Defaults to
#'   \code{"D1"}, \code{"D2"}, \ldots
#'
#' @details alpha and k are prescribed internal state variables, not inferred parameters:
#' they are computed from forcing, pore structure and memory at the start of
#' each time step. Alpha multiplies selected Fe/Mn/N/S/C expressions as an
#' accessibility factor, whereas k_accept directly gates only the model's
#' crystalline-Fe reduction expression. A low alpha suppresses selected rates;
#' restoring alpha changes selected rates, but does not guarantee functional recovery.
#' k additionally gates the crystalline-Fe exchange rate. The latent_truth
#' vector is a constructed index sharing ingredients with soil_data;
#' omit those ingredients from observable-only benchmarks. history_pair denotes
#' a shared random effect, not experimentally matched twins. Low n_time can leave
#' no adequate baseline or recovery; inspect forcing and analyse each event separately.
#' Mineral crystallisation is continuous in this implementation, not an
#' event-only reoxidation ratchet. k does not explicitly decrease with crystallinity.
#' ROL, ROS and several microbial features are synthetic descriptors, not calibrated
#' physical fluxes. Do not infer quantitative field rates from their labels alone.
#' Q_accept now counts one electron per crystalline Fe(III); the previous 0.22
#' factor mixed accessibility into inventory. DOC reducing equivalents remain
#' an illustrative coefficient dependent on assumed carbon oxidation state.
#' @return List of identifiers, data blocks, latent states, flux descriptors,
#' balance checks, metadata and legacy views. The metadata records hidden columns,
#' seed and RNG configuration. With seed supplied, the caller RNG is restored.
#' @importFrom stats rnorm runif rlnorm rnbinom
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(n_plot = 2, n_depth = 2, n_plant = 2,
#'                                    n_time = 20, seed = 42)
#'   nrow(sim$id)  # 2 x 2 x 2 x 20 = 160 rows
#'   names(sim)    # top-level list elements
#' }
#' @export
simulate_redox_holobiont <- function(
  n_plot = 4,
  n_depth = 2,
  n_plant = 6,
  n_time = 30,
  p_micro = 60,
  seed = 123,
  scenario = c("flood_drain", "drought_rewet"),
  n_cycles = 2L,
  disturbance_strength = 0.65,
  disturbance_center = NULL,
  disturbance_width = 0.08,
  seasonal_amp = 0.08,
  seasonal_phase = 0,
  history_strength = 0.55,
  rescue = c("none", "capacity", "connectivity", "kinetics"),
  event_tau_h = 24,
  sequencing_depth = 2e5,
  metat_depth = 5e5,
  decoupling = 0.25,
  zero_inflation = 0.20,
  MNAR_strength = 0.30,
  Eh_dropout_threshold = 100,
  micro_mean = 8,
  micro_slope = 3,
  micro_lambda_min = 1e-8,
  micro_lambda_max = 1e6,
  stochastic_reassembly = TRUE,
  include_graph = FALSE,
  depth_labels = NULL
) {
  clamp <- function(x, lo = 0, hi = 1) pmin(hi, pmax(lo, x))
  logistic <- function(x) 1 / (1 + exp(-x))
  assert_integer <- function(x, name, minimum = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x > .Machine$integer.max || x != floor(x)) {
      stop("`", name, "` must be one integer >= ", minimum, ".", call. = FALSE)
    }
  }
  assert_scalar <- function(x, name, lower = -Inf, upper = Inf,
                            lower_open = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
      if (lower_open) x > lower else x >= lower
    ok <- ok && x <= upper
    if (!ok) stop("Invalid `", name, "`.", call. = FALSE)
  }

  if (!is.null(seed)) {
    assert_integer(seed, "seed", 0L)
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(as.integer(seed))
  }
  assert_integer(n_plot, "n_plot")
  assert_integer(n_depth, "n_depth")
  assert_integer(n_plant, "n_plant")
  assert_integer(n_time, "n_time", 4L)
  assert_integer(p_micro, "p_micro")
  assert_integer(n_cycles, "n_cycles")
  assert_scalar(disturbance_strength, "disturbance_strength", 0, 1)
  assert_scalar(disturbance_width, "disturbance_width", 0, 1, TRUE)
  assert_scalar(seasonal_amp, "seasonal_amp", 0, Inf)
  assert_scalar(seasonal_phase, "seasonal_phase")
  assert_scalar(history_strength, "history_strength", 0, 1)
  assert_scalar(decoupling, "decoupling", 0, 1)
  assert_scalar(zero_inflation, "zero_inflation", 0, 1)
  assert_scalar(MNAR_strength, "MNAR_strength", 0, 1)
  assert_scalar(Eh_dropout_threshold, "Eh_dropout_threshold")
  assert_scalar(event_tau_h, "event_tau_h", 0, Inf, TRUE)
  assert_scalar(sequencing_depth, "sequencing_depth", 0, Inf, TRUE)
  assert_scalar(metat_depth, "metat_depth", 0, Inf, TRUE)
  assert_scalar(micro_mean, "micro_mean", 0, Inf, TRUE)
  assert_scalar(micro_slope, "micro_slope", 0, Inf)
  assert_scalar(micro_lambda_min, "micro_lambda_min", 0, Inf, TRUE)
  assert_scalar(micro_lambda_max, "micro_lambda_max", 0, Inf, TRUE)
  if (micro_lambda_max <= micro_lambda_min) {
    stop("`micro_lambda_max` must exceed `micro_lambda_min`.", call. = FALSE)
  }
  if (!is.logical(stochastic_reassembly) || length(stochastic_reassembly) != 1L ||
    is.na(stochastic_reassembly)) {
    stop("`stochastic_reassembly` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_graph) || length(include_graph) != 1L ||
    is.na(include_graph)) {
    stop("`include_graph` must be TRUE or FALSE.", call. = FALSE)
  }

  scenario <- match.arg(scenario)
  rescue <- match.arg(rescue)
  if (!is.null(depth_labels) &&
    (!is.character(depth_labels) || length(depth_labels) != n_depth ||
      any(!nzchar(depth_labels)) || anyDuplicated(depth_labels))) {
    stop("`depth_labels` must contain n_depth unique non-empty labels.",
      call. = FALSE
    )
  }

  depth_levels <- if (is.null(depth_labels)) {
    paste0("D", seq_len(n_depth))
  } else {
    depth_labels
  }
  id <- expand.grid(
    plot = paste0("P", seq_len(n_plot)),
    depth = depth_levels,
    plant_id = paste0("Plant", seq_len(n_plant)),
    time = seq_len(n_time),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  n <- nrow(id)
  id$row_id <- seq_len(n)
  id$unit_id <- interaction(id$plot, id$depth, id$plant_id,
    drop = TRUE, lex.order = TRUE
  )
  plant_index <- match(id$plant_id, paste0("Plant", seq_len(n_plant)))
  id$history_pair <- paste(id$plot, id$depth,
    paste0("Pair", ceiling(plant_index / 2)),
    sep = "_"
  )
  id$history <- ifelse(plant_index %% 2 == 0, "preconditioned", "naive")
  id$scenario <- scenario
  id$rescue <- rescue

  # Event forcing. Multiple cycles are represented by the maximum of separated
  # Gaussian pulses; rewetting/drainage recovery is a lagged, wider pulse.
  time_values <- seq_len(n_time)
  width_time <- max(0.75, disturbance_width * n_time)
  if (is.null(disturbance_center)) {
    centers <- seq(0.22 * n_time, 0.78 * n_time, length.out = n_cycles)
  } else {
    assert_scalar(disturbance_center, "disturbance_center")
    spacing <- if (n_cycles > 1L) max(2 * width_time, n_time / (n_cycles + 1)) else 0
    centers <- disturbance_center + spacing * (seq_len(n_cycles) - 1L)
    centers <- clamp(centers, 1, n_time)
  }
  pulse_matrix <- vapply(centers, function(cc) {
    exp(-0.5 * ((time_values - cc) / width_time)^2)
  }, numeric(n_time))
  if (is.null(dim(pulse_matrix))) pulse_matrix <- matrix(pulse_matrix, ncol = 1L)
  event_time <- apply(pulse_matrix, 1, max)
  cycle_time <- max.col(pulse_matrix, ties.method = "first")
  recovery_matrix <- vapply(centers + 1.6 * width_time, function(cc) {
    exp(-0.5 * ((time_values - cc) / (1.35 * width_time))^2)
  }, numeric(n_time))
  if (is.null(dim(recovery_matrix))) recovery_matrix <- matrix(recovery_matrix, ncol = 1L)
  recovery_time <- apply(recovery_matrix, 1, max)
  phase_time <- ifelse(event_time >= 0.35, "disturbance",
    ifelse(recovery_time >= 0.25, "recovery", "baseline")
  )

  depth_index <- match(id$depth, depth_levels)
  depth_scaled <- (depth_index - 1) / max(1, n_depth - 1)
  plot_index <- match(id$plot, paste0("P", seq_len(n_plot)))
  event <- event_time[id$time]
  recovery <- recovery_time[id$time]
  seasonal <- seasonal_amp * sin(2 * pi * id$time / n_time + seasonal_phase)
  plot_wet <- stats::rnorm(n_plot, 0, 0.035)[plot_index]
  base_wfps <- 0.57 + 0.15 * depth_scaled + plot_wet + seasonal

  if (scenario == "flood_drain") {
    WFPS <- base_wfps + 0.33 * disturbance_strength * event -
      0.08 * disturbance_strength * recovery
  } else {
    WFPS <- base_wfps - 0.34 * disturbance_strength * event +
      0.28 * disturbance_strength * recovery
  }
  WFPS <- clamp(WFPS, 0.18, 0.99)
  bulk_density <- clamp(1.15 + 0.18 * depth_scaled +
    stats::rnorm(n, 0, 0.025), 0.8, 1.65)
  porosity <- clamp(1 - bulk_density / 2.65, 0.30, 0.70)
  air_filled_porosity <- porosity * (1 - WFPS)
  pore_connectivity <- clamp(
    exp(-((WFPS - 0.68) / 0.28)^2) *
      (0.55 + 0.45 * air_filled_porosity / pmax(porosity, 1e-6)), 0, 1
  )
  temperature_C <- 20 + 4 * sin(2 * pi * id$time / n_time - 0.5) -
    1.5 * depth_scaled + stats::rnorm(n, 0, 0.35)
  matric_potential_kPa <- -0.8 * exp(6.5 * (0.82 - WFPS))
  water_table_cm <- 35 - 70 * WFPS + 8 * depth_scaled

  forcing <- data.frame(
    event_intensity = event,
    recovery_intensity = recovery,
    WFPS = WFPS,
    air_filled_porosity = air_filled_porosity,
    matric_potential_kPa = matric_potential_kPa,
    water_table_cm = water_table_cm,
    temperature_C = temperature_C,
    pore_connectivity = pore_connectivity,
    cycle = cycle_time[id$time],
    phase = phase_time[id$time],
    stringsAsFactors = FALSE
  )
  id$cycle <- forcing$cycle
  id$phase <- forcing$phase
  id$event_intensity <- forcing$event_intensity
  id$WFPS <- forcing$WFPS
  id$water_table_cm <- forcing$water_table_cm

  # Allocate longitudinal states.
  state_names <- c(
    "FeIII_poor_crystalline", "FeIII_crystalline", "FeII",
    "FeS", "MnIV", "MnIII", "MnII", "NO3", "NH4", "SO4",
    "sulfide", "CH4", "DOC", "humic_EAC", "humic_EDC", "memory",
    "root_biomass", "aerenchyma"
  )
  state <- matrix(NA_real_,
    nrow = n, ncol = length(state_names),
    dimnames = list(NULL, state_names)
  )
  flux_names <- c(
    "FeIII_reduction", "FeII_oxidation", "Fe_crystallisation",
    "FeS_formation", "FeS_oxidation", "MnIV_to_MnIII",
    "MnIII_to_MnII", "MnII_oxidation", "denitrification", "DNRA",
    "nitrification", "sulfate_reduction", "sulfide_oxidation",
    "methanogenesis", "methane_oxidation", "CO2_flux", "N2O_flux",
    "CH4_flux"
  )
  flux_mat <- matrix(0,
    nrow = n, ncol = length(flux_names),
    dimnames = list(NULL, flux_names)
  )
  plant_names <- c(
    "SPAD", "FvFm", "PhiPSII", "NPQ", "ROL",
    "root_biomass", "root_length_density", "root_porosity",
    "aerenchyma", "Fe_plaque", "Mn_plaque", "ROS_load"
  )
  plant_mat <- matrix(NA_real_,
    nrow = n, ncol = length(plant_names),
    dimnames = list(NULL, plant_names)
  )
  process_names <- c(
    "porewater_O2", "reduction_signal", "oxidation_signal",
    "alpha_accept", "alpha_donate", "k_accept", "k_donate",
    "Q_accept", "Q_donate", "Cacc_EAC", "Cacc_EDC",
    "Cacc_total", "Cacc_fraction", "net_oxidative_balance",
    "Eh", "pH"
  )
  process <- matrix(NA_real_,
    nrow = n, ncol = length(process_names),
    dimnames = list(NULL, process_names)
  )
  fe_error_all <- rep(NA_real_, n)
  mn_error_all <- rep(NA_real_, n)

  unit_levels <- levels(id$unit_id)
  pair_levels <- unique(id$history_pair)
  pair_random <- stats::rnorm(length(pair_levels), 0, 1)
  names(pair_random) <- pair_levels
  plant_random <- stats::rnorm(length(unit_levels), 0, 0.12)
  names(plant_random) <- unit_levels

  # State-transition loop for each plant-depth trajectory.
  for (uu in unit_levels) {
    ix <- which(id$unit_id == uu)
    ix <- ix[order(id$time[ix])]
    first <- ix[1]
    dep <- depth_scaled[first]
    pair_eff <- unname(pair_random[id$history_pair[first]])
    unit_eff <- unname(plant_random[uu])
    preconditioned <- id$history[first] == "preconditioned"

    fe_total0 <- max(40, 145 + 48 * dep + 12 * pair_eff)
    mn_total0 <- max(4, 18 + 7 * dep + 2 * pair_eff)
    s <- c(
      FeIII_poor_crystalline = 0.48 * fe_total0,
      FeIII_crystalline = 0.38 * fe_total0,
      FeII = 0.12 * fe_total0,
      FeS = 0.02 * fe_total0,
      MnIV = 0.58 * mn_total0,
      MnIII = 0.20 * mn_total0,
      MnII = 0.22 * mn_total0,
      NO3 = 5.5 + stats::runif(1, 0, 2),
      NH4 = 3.0 + stats::runif(1, 0, 1.5),
      SO4 = 9.0 + stats::runif(1, 0, 3),
      sulfide = 0.3 + stats::runif(1, 0, 0.3),
      CH4 = 0.08 + stats::runif(1, 0, 0.08),
      DOC = max(3, 15 + 4 * pair_eff),
      humic_EAC = max(3, 22 + 4 * pair_eff),
      humic_EDC = max(3, 14 + 3 * pair_eff),
      memory = if (preconditioned) 0.30 * history_strength else 0.03,
      root_biomass = max(0.12, 0.42 + 0.06 * pair_eff + 0.04 * unit_eff),
      aerenchyma = if (preconditioned) 0.24 else 0.12
    )
    if (rescue == "capacity") {
      s["FeIII_poor_crystalline"] <- s["FeIII_poor_crystalline"] * 1.20
      s["humic_EAC"] <- s["humic_EAC"] * 1.20
    }
    fe_total_initial <- sum(s[c(
      "FeIII_poor_crystalline", "FeIII_crystalline",
      "FeII", "FeS"
    )])
    mn_total_initial <- sum(s[c("MnIV", "MnIII", "MnII")])

    for (jj in seq_along(ix)) {
      rr <- ix[jj]
      wet_stress <- clamp((forcing$WFPS[rr] - 0.78) / 0.18)
      dry_stress <- clamp((0.43 - forcing$WFPS[rr]) / 0.25)
      stress <- clamp(wet_stress + dry_stress)

      aerenchyma_target <- clamp(0.10 + 0.62 * wet_stress +
        0.20 * s["memory"])
      s["aerenchyma"] <- clamp(0.78 * s["aerenchyma"] +
        0.22 * aerenchyma_target)
      growth <- 0.025 * (1 - stress) * (1 - jj / (1.6 * length(ix)))
      loss <- 0.020 * stress * (1 + 0.6 * s["memory"])
      s["root_biomass"] <- max(0.04, s["root_biomass"] + growth - loss)
      root_porosity <- clamp(0.11 + 0.48 * s["aerenchyma"], 0.08, 0.58)
      root_length_density <- 0.65 + 4.2 * s["root_biomass"]
      ROL <- max(0, 0.12 + 2.5 * s["root_biomass"] * root_porosity *
        (0.35 + forcing$air_filled_porosity[rr] /
          pmax(forcing$air_filled_porosity[rr] + 0.06, 1e-6)))

      porewater_O2 <- clamp(0.29 * (1 - forcing$WFPS[rr])^1.7 +
        0.045 * ROL, 0.001, 0.32)
      ox_signal <- clamp(porewater_O2 / 0.28 + 0.10 * ROL, 0, 1.4)
      red_signal <- clamp(forcing$WFPS[rr]^2.2 *
        (1 - porewater_O2 / 0.32) *
        (s["DOC"] / (s["DOC"] + 12)), 0, 1)
      red_signal <- clamp(red_signal + stats::rnorm(1, 0, 0.03 * decoupling))
      ox_signal <- clamp(ox_signal + stats::rnorm(1, 0, 0.03 * decoupling), 0, 1.5)

      # --- Causal state variables: alpha and k computed BEFORE process transitions ---
      # alpha (connectivity) = fraction of each pool connected to operative
      # electron-transfer pathways under current pore topology. Low alpha means
      # donors cannot reach acceptors regardless of thermodynamic favourability.
      # k (kinetics) = effective exchange rate. Low k means connected but
      # slowly exchanging pools (e.g., crystalline Fe) contribute little per step.
      # Both are computed here from current forcing and carry-over memory, then
      # used as multipliers in selected process-rate expressions below.
      alpha_accept <- clamp(0.12 + 0.57 * forcing$pore_connectivity[rr] +
        0.10 * ROL - 0.20 * s["memory"])
      alpha_donate <- clamp(0.16 + 0.48 * forcing$pore_connectivity[rr] +
        0.18 * red_signal - 0.10 * s["memory"])
      if (rescue == "connectivity") {
        alpha_accept <- clamp(alpha_accept + 0.20)
        alpha_donate <- clamp(alpha_donate + 0.12)
      }
      k_accept <- 0.012 + 0.16 * ox_signal + 0.045 * red_signal
      k_donate <- 0.010 + 0.11 * red_signal + 0.055 * ox_signal
      if (rescue == "kinetics") {
        k_accept <- 1.6 * k_accept
        k_donate <- 1.4 * k_donate
      }
      # Normalised kinetics factor: 1.0 at reference exchange rate (0.20 h^-1),
      # used to additionally gate crystalline-Fe exchange (kinetic barrier).
      k_norm <- clamp(k_accept / 0.20)

      # Fe transitions. Simultaneous outflows are scaled to the source pool,
      # which preserves total Fe across Fe(III), Fe(II), and FeS.
      # alpha_accept gates acceptor-side transitions (Fe(III) reduction, FeS oxidation);
      # k_norm additionally gates crystalline-Fe exchange (kinetic barrier).
      # alpha_accept also gates Fe(II) oxidation: Fe(II) must reach O2-bearing domains.
      fe_red_pc <- s["FeIII_poor_crystalline"] *
        (0.006 + 0.095 * red_signal) * alpha_accept
      fe_red_cr <- s["FeIII_crystalline"] *
        (0.001 + 0.012 * red_signal) * alpha_accept * k_norm
      fe_cryst <- s["FeIII_poor_crystalline"] *
        (0.001 + 0.007 * s["memory"] * red_signal)
      out_pc <- fe_red_pc + fe_cryst
      if (out_pc > s["FeIII_poor_crystalline"]) {
        scale_pc <- s["FeIII_poor_crystalline"] / out_pc
        fe_red_pc <- fe_red_pc * scale_pc
        fe_cryst <- fe_cryst * scale_pc
      }
      fe_red_cr <- min(fe_red_cr, s["FeIII_crystalline"])
      fe_ox <- s["FeII"] * (0.006 + 0.16 * ox_signal) * alpha_accept
      fes_form <- s["FeII"] * 0.025 * red_signal *
        s["sulfide"] / (s["sulfide"] + 1)
      out_fe2 <- fe_ox + fes_form
      if (out_fe2 > s["FeII"]) {
        scale_fe2 <- s["FeII"] / out_fe2
        fe_ox <- fe_ox * scale_fe2
        fes_form <- fes_form * scale_fe2
      }
      fes_ox <- min(s["FeS"], s["FeS"] * 0.08 * ox_signal * alpha_accept)
      s["FeIII_poor_crystalline"] <- s["FeIII_poor_crystalline"] -
        fe_red_pc - fe_cryst + fe_ox + fes_ox
      s["FeIII_crystalline"] <- s["FeIII_crystalline"] - fe_red_cr + fe_cryst
      s["FeII"] <- s["FeII"] + fe_red_pc + fe_red_cr - fe_ox - fes_form
      s["FeS"] <- s["FeS"] + fes_form - fes_ox

      # Mn(IV) <-> Mn(III) <-> Mn(II), also mass conserving.
      # alpha_accept gates Mn(IV/III) reduction and Mn(II/III) oxidation:
      # manganese-reducing and -oxidising bacteria require donor-acceptor contact.
      mn_red_iv <- min(s["MnIV"],
        s["MnIV"] * (0.007 + 0.075 * red_signal) * alpha_accept)
      mn_red_iii <- s["MnIII"] * (0.005 + 0.065 * red_signal) * alpha_accept
      mn_ox_iii  <- s["MnIII"] * (0.003 + 0.075 * ox_signal) * alpha_accept
      if (mn_red_iii + mn_ox_iii > s["MnIII"]) {
        scale_mn3 <- s["MnIII"] / (mn_red_iii + mn_ox_iii)
        mn_red_iii <- mn_red_iii * scale_mn3
        mn_ox_iii <- mn_ox_iii * scale_mn3
      }
      mn_ox_ii <- min(s["MnII"],
        s["MnII"] * (0.004 + 0.11 * ox_signal) * alpha_accept)
      s["MnIV"] <- s["MnIV"] - mn_red_iv + mn_ox_iii
      s["MnIII"] <- s["MnIII"] + mn_red_iv - mn_red_iii + mn_ox_ii - mn_ox_iii
      s["MnII"] <- s["MnII"] + mn_red_iii - mn_ox_ii

      # N, S, methane, and dissolved organic pools. Rates are recorded as
      # process fluxes; pools are bounded at zero.
      # alpha_accept gates nitrification (O2 connectivity), denitrification/DNRA
      # (NO3 connectivity), and sulfate reduction (SO4 connectivity).
      # alpha_donate gates methanogenesis (donor connectivity to methanogens).
      denit <- min(s["NO3"],
        0.24 * red_signal * s["NO3"] / (s["NO3"] + 2) * alpha_accept)
      dnra <- min(
        s["NO3"] - denit,
        0.11 * red_signal * s["NO3"] / (s["NO3"] + 2) * alpha_accept
      )
      nitrif <- min(s["NH4"],
        0.16 * ox_signal * s["NH4"] / (s["NH4"] + 2) * alpha_accept)
      s["NO3"] <- max(0, s["NO3"] - denit - dnra + nitrif)
      s["NH4"] <- max(0, s["NH4"] + dnra - nitrif + 0.025 * s["DOC"])
      sulfate_red <- min(s["SO4"],
        0.10 * red_signal * s["SO4"] / (s["SO4"] + 3) * alpha_accept)
      sulfide_ox <- min(s["sulfide"],
        0.16 * ox_signal * s["sulfide"] / (s["sulfide"] + 0.8) * alpha_accept)
      s["SO4"] <- max(0, s["SO4"] - sulfate_red + sulfide_ox)
      s["sulfide"] <- max(0, s["sulfide"] + sulfate_red - sulfide_ox -
        0.5 * fes_form + 0.5 * fes_ox)
      methanogenesis <- 0.10 * red_signal^2 * s["DOC"] / (s["DOC"] + 10) *
        1 / (1 + s["NO3"] + 0.25 * s["SO4"]) * alpha_donate
      methane_oxidation <- min(s["CH4"],
        0.25 * ox_signal * s["CH4"] / (s["CH4"] + 0.2) * alpha_accept)
      s["CH4"] <- max(0, s["CH4"] + methanogenesis - methane_oxidation)
      carbon_use <- 0.12 * red_signal + 0.05 * ox_signal + 0.2 * methanogenesis
      s["DOC"] <- max(1, s["DOC"] - carbon_use + 0.10 * s["root_biomass"])
      s["humic_EAC"] <- max(0.1, s["humic_EAC"] + 0.04 * ox_signal -
        0.05 * red_signal)
      s["humic_EDC"] <- max(0.1, s["humic_EDC"] + 0.05 * red_signal -
        0.04 * ox_signal)

      disturbance_load <- clamp(abs(red_signal - 0.35) * event[rr] +
        0.45 * stress)
      memory_gain <- history_strength * (0.035 * disturbance_load +
        0.025 * fe_cryst /
          pmax(fe_total_initial, 1e-6))
      memory_decay <- 0.018 * (1 - event[rr]) * (1 - s["FeIII_crystalline"] /
        fe_total_initial)
      s["memory"] <- clamp(s["memory"] + memory_gain - max(0, memory_decay))

      Q_accept <- s["FeIII_poor_crystalline"] + s["FeIII_crystalline"] +
        2 * s["MnIV"] + s["MnIII"] + 5 * s["NO3"] +
        s["humic_EAC"] + 8 * s["SO4"]
      Q_donate <- s["FeII"] + s["FeS"] * 9 + 2 * s["MnII"] + s["MnIII"] +
        8 * s["sulfide"] + 8 * s["CH4"] + s["humic_EDC"] +
        0.5 * s["DOC"] + 8 * s["NH4"]
      Cacc_EAC <- Q_accept * alpha_accept * (1 - exp(-k_accept * event_tau_h))
      Cacc_EDC <- Q_donate * alpha_donate * (1 - exp(-k_donate * event_tau_h))
      Cacc_total <- Cacc_EAC + Cacc_EDC
      total_Q <- Q_accept + Q_donate
      Cacc_fraction <- Cacc_total / pmax(total_Q, 1e-9)
      net_balance <- Cacc_EAC - Cacc_EDC
      Eh <- 470 * (porewater_O2 / 0.32) - 185 * log1p(Q_donate / pmax(Q_accept, 1)) -
        90 * red_signal + 35 + stats::rnorm(1, 0, 12)
      pH <- clamp(6.65 - 0.28 * red_signal + 0.12 * ox_signal -
        0.08 * sulfide_ox + stats::rnorm(1, 0, 0.06), 4.5, 8.5)

      ROS_load <- clamp(0.18 + 0.70 * stress + 0.30 * s["memory"] -
        0.12 * s["aerenchyma"] + stats::rnorm(1, 0, 0.04), 0, 1.5)
      SPAD <- 44 - 11 * ROS_load + 1.8 * s["root_biomass"] +
        stats::rnorm(1, 0, 0.9)
      FvFm <- clamp(0.835 - 0.18 * ROS_load + stats::rnorm(1, 0, 0.009), 0.35, 0.86)
      PhiPSII <- clamp(0.52 - 0.26 * ROS_load + stats::rnorm(1, 0, 0.015), 0.05, 0.65)
      NPQ <- max(0.05, 0.55 + 1.45 * ROS_load + stats::rnorm(1, 0, 0.08))
      Fe_plaque <- max(0, 0.24 * ROL * s["FeII"] * (0.5 + ox_signal))
      Mn_plaque <- max(0, 0.16 * ROL * s["MnII"] * (0.5 + ox_signal))

      state[rr, ] <- s[state_names]
      plant_mat[rr, ] <- c(
        SPAD, FvFm, PhiPSII, NPQ, ROL,
        s["root_biomass"], root_length_density, root_porosity,
        s["aerenchyma"], Fe_plaque, Mn_plaque, ROS_load
      )
      process[rr, ] <- c(
        porewater_O2, red_signal, ox_signal, alpha_accept,
        alpha_donate, k_accept, k_donate, Q_accept, Q_donate,
        Cacc_EAC, Cacc_EDC, Cacc_total, Cacc_fraction,
        net_balance, Eh, pH
      )
      flux_mat[rr, ] <- c(
        fe_red_pc + fe_red_cr, fe_ox, fe_cryst, fes_form, fes_ox,
        mn_red_iv, mn_red_iii, mn_ox_ii + mn_ox_iii, denit, dnra,
        nitrif, sulfate_red, sulfide_ox, methanogenesis, methane_oxidation,
        0.18 * red_signal * s["DOC"] + methane_oxidation,
        0.018 * denit * (1 + s["memory"]),
        max(0, methanogenesis - methane_oxidation)
      )
    }

    fe_error <- rowSums(state[ix, c(
      "FeIII_poor_crystalline", "FeIII_crystalline",
      "FeII", "FeS"
    ), drop = FALSE]) - fe_total_initial
    mn_error <- rowSums(state[ix, c("MnIV", "MnIII", "MnII"), drop = FALSE]) -
      mn_total_initial
    fe_error_all[ix] <- fe_error
    mn_error_all[ix] <- mn_error
  }

  # Synthetic DNA and MetaT observations. DNA counts are generated from current
  # state suitability and memory; this is not an explicit population-growth or
  # slow-turnover model. Transcript counts respond to current suitability,
  # conditional on the generated DNA count and library size; neither is flux.
  gene_names <- c(
    "mtrA", "omcS", "cyc2", "mnxG", "narG", "napA", "nirK",
    "nirS", "norB", "nosZ", "nrfA", "amoA_AOA", "amoA_AOB",
    "nxrB", "dsrA", "dsrB", "mcrA", "pmoA"
  )
  red <- process[, "reduction_signal"]
  ox <- clamp(process[, "oxidation_signal"])
  mem <- state[, "memory"]
  activity <- cbind(
    mtrA = clamp(0.08 + 0.85 * red * state[, "FeIII_poor_crystalline"] /
      (state[, "FeIII_poor_crystalline"] + 40)),
    omcS = clamp(0.06 + 0.78 * red * forcing$pore_connectivity),
    cyc2 = clamp(0.05 + 0.85 * ox * state[, "FeII"] / (state[, "FeII"] + 15)),
    mnxG = clamp(0.04 + 0.78 * ox * state[, "MnII"] / (state[, "MnII"] + 4)),
    narG = clamp(0.05 + 0.82 * red * state[, "NO3"] / (state[, "NO3"] + 2)),
    napA = clamp(0.07 + 0.65 * red * state[, "NO3"] / (state[, "NO3"] + 1)),
    nirK = clamp(0.04 + 0.76 * red),
    nirS = clamp(0.04 + 0.80 * red),
    norB = clamp(0.04 + 0.72 * red),
    nosZ = clamp(0.06 + 0.62 * red * (1 - 0.35 * mem)),
    nrfA = clamp(0.03 + 0.74 * red * state[, "DOC"] / (state[, "DOC"] + 10)),
    amoA_AOA = clamp(0.05 + 0.72 * ox * state[, "NH4"] / (state[, "NH4"] + 2)),
    amoA_AOB = clamp(0.04 + 0.80 * ox * state[, "NH4"] / (state[, "NH4"] + 3)),
    nxrB = clamp(0.04 + 0.68 * ox),
    dsrA = clamp(0.02 + 0.86 * red * state[, "SO4"] / (state[, "SO4"] + 4)),
    dsrB = clamp(0.02 + 0.84 * red * state[, "SO4"] / (state[, "SO4"] + 4)),
    mcrA = clamp(0.01 + 0.92 * red^2 /
      (1 + state[, "NO3"] + 0.2 * state[, "SO4"])),
    pmoA = clamp(0.02 + 0.90 * ox * state[, "CH4"] / (state[, "CH4"] + 0.15))
  )
  # as.vector() + matrix() guarantees a proper 2-D matrix even when pmin/pmax
  # drops the dim attribute in some R builds.
  activity <- matrix(
    clamp(
      as.vector(activity) +
        stats::rnorm(n * length(gene_names), 0, 0.035 * decoupling)
    ),
    nrow = n, ncol = length(gene_names),
    dimnames = list(NULL, gene_names)
  )

  gene_base <- stats::runif(length(gene_names), 4.3, 5.8)
  names(gene_base) <- gene_names
  dna_mat <- matrix(0,
    nrow = n, ncol = length(gene_names),
    dimnames = list(NULL, gene_names)
  )
  for (gg in seq_along(gene_names)) {
    log10_mu <- gene_base[gg] + 0.70 * (activity[, gg] - 0.4) +
      0.28 * mem + if (isTRUE(stochastic_reassembly)) {
        stats::rnorm(n, 0, 0.10 + 0.08 * decoupling)
      } else {
        0
      }
    mu <- clamp(10^log10_mu, 10, 5e8)
    dna_mat[, gg] <- stats::rnbinom(n, mu = mu, size = 18)
  }
  micro_gene_abundance <- as.data.frame(dna_mat, check.names = FALSE)

  metat_library_target <- pmax(1000, round(stats::rlnorm(
    n, log(metat_depth) - 0.5 * 0.28^2, 0.28
  )))
  metat_mat <- matrix(0,
    nrow = n, ncol = length(gene_names),
    dimnames = list(NULL, gene_names)
  )
  for (rr in seq_len(n)) {
    weights <- sqrt(dna_mat[rr, ] + 1) * exp(3.2 * activity[rr, ])
    mu <- metat_library_target[rr] * weights / sum(weights)
    metat_mat[rr, ] <- stats::rnbinom(length(gene_names), mu = mu, size = 10)
  }
  micro_metat_counts <- as.data.frame(metat_mat, check.names = FALSE)
  raw_total_ratio <- rowSums(metat_mat) / pmax(rowSums(dna_mat), 1)
  micro_metat_metadata <- data.frame(
    target_library_size = metat_library_target,
    observed_library_size = rowSums(metat_mat),
    raw_total_RNA_DNA_ratio = raw_total_ratio,
    RNA_DNA_activity_ratio = raw_total_ratio
  )

  micro_traits <- data.frame(
    EET_reduction = rowMeans(activity[, c("mtrA", "omcS"), drop = FALSE]),
    Fe_oxidation = activity[, "cyc2"],
    Mn_oxidation = activity[, "mnxG"],
    denitrification = rowMeans(activity[, c(
      "narG", "napA", "nirK", "nirS",
      "norB", "nosZ"
    ), drop = FALSE]),
    DNRA = activity[, "nrfA"],
    nitrification = rowMeans(activity[, c("amoA_AOA", "amoA_AOB", "nxrB"),
      drop = FALSE
    ]),
    sulfate_reduction = rowMeans(activity[, c("dsrA", "dsrB"), drop = FALSE]),
    methanogenesis = activity[, "mcrA"],
    methane_oxidation = activity[, "pmoA"]
  )

  # Optional ASV-like community for backward compatibility. It is generated
  # from heterogeneous redox niches, not used to create the hidden state.
  tax_library_target <- pmax(100, round(stats::rlnorm(
    n, log(sequencing_depth * micro_mean / 8) - 0.5 * 0.35^2, 0.35
  )))
  niche_optimum <- stats::runif(p_micro, 0, 1)
  niche_width <- stats::runif(p_micro, 0.12, 0.35)
  redox_position <- clamp(0.65 * red + 0.25 * mem + 0.10 * (1 - ox))
  micro_mat <- matrix(0, nrow = n, ncol = p_micro)
  for (rr in seq_len(n)) {
    niche <- exp(-0.5 * ((redox_position[rr] - niche_optimum) / niche_width)^2)
    niche <- niche * exp(stats::rnorm(
      p_micro, 0,
      if (stochastic_reassembly) 0.35 else 0.12
    ))
    mu <- tax_library_target[rr] * niche / sum(niche)
    mu <- clamp(mu, micro_lambda_min, micro_lambda_max)
    micro_mat[rr, ] <- stats::rnbinom(p_micro,
      mu = mu, size =
        max(0.5, 5 / (1 + micro_slope))
    )
  }
  if (zero_inflation > 0) {
    zeros <- matrix(stats::runif(n * p_micro) < zero_inflation,
      nrow = n, ncol = p_micro
    )
    micro_mat[zeros] <- 0
  }
  micro_data <- as.data.frame(micro_mat)
  names(micro_data) <- paste0("ASV", seq_len(p_micro))

  # Soil data and legacy views.
  soil_data <- data.frame(
    Eh = process[, "Eh"],
    pH = process[, "pH"],
    WFPS = forcing$WFPS,
    porewater_O2_mmol_L = process[, "porewater_O2"],
    O2_supply_mmol_kg = process[, "porewater_O2"] *
      (porosity * forcing$WFPS) / bulk_density,
    FeIII_poor_crystalline_mmol_kg = state[, "FeIII_poor_crystalline"],
    FeIII_crystalline_mmol_kg = state[, "FeIII_crystalline"],
    FeII_mmol_kg = state[, "FeII"],
    FeS_mmol_kg = state[, "FeS"],
    MnIV_mmol_kg = state[, "MnIV"],
    MnIII_mmol_kg = state[, "MnIII"],
    MnII_mmol_kg = state[, "MnII"],
    NO3_mmol_kg = state[, "NO3"],
    NH4_mmol_kg = state[, "NH4"],
    SO4_mmol_kg = state[, "SO4"],
    sulfide_mmol_kg = state[, "sulfide"],
    CH4_mmol_kg = state[, "CH4"],
    DOC_mmolC_kg = state[, "DOC"],
    humic_EAC_mmol_e_kg = state[, "humic_EAC"],
    humic_EDC_mmol_e_kg = state[, "humic_EDC"],
    EAC = process[, "Q_accept"],
    EDC = process[, "Q_donate"],
    Cacc_EAC = process[, "Cacc_EAC"],
    Cacc_EDC = process[, "Cacc_EDC"],
    Cacc_total = process[, "Cacc_total"],
    Cacc_fraction = process[, "Cacc_fraction"],
    net_oxidative_balance = process[, "net_oxidative_balance"],
    alpha_accept = process[, "alpha_accept"],
    alpha_donate = process[, "alpha_donate"],
    k_accept_h = process[, "k_accept"],
    k_donate_h = process[, "k_donate"],
    pore_connectivity = forcing$pore_connectivity,
    bulk_density_g_cm3 = bulk_density,
    porosity = porosity,
    Fe2.Fe3 = state[, "FeII"] /
      pmax(state[, "FeIII_poor_crystalline"] + state[, "FeIII_crystalline"], 1e-9),
    Mn2.Mn4 = state[, "MnII"] / pmax(state[, "MnIV"], 1e-9),
    NH4.NO3 = state[, "NH4"] / pmax(state[, "NO3"], 1e-9)
  )
  dropout_prob <- clamp(MNAR_strength * logistic((Eh_dropout_threshold -
    soil_data$Eh) / 35))
  soil_data$Eh[stats::runif(n) < dropout_prob] <- NA_real_
  Eh_stability <- soil_data
  plant_data <- as.data.frame(plant_mat)
  fluxes <- as.data.frame(flux_mat)

  latent_state <- data.frame(
    Q_accept = process[, "Q_accept"],
    Q_donate = process[, "Q_donate"],
    alpha_accept = process[, "alpha_accept"],
    alpha_donate = process[, "alpha_donate"],
    k_accept_h = process[, "k_accept"],
    k_donate_h = process[, "k_donate"],
    memory = state[, "memory"],
    redox_position = redox_position,
    Cacc_EAC = process[, "Cacc_EAC"],
    Cacc_EDC = process[, "Cacc_EDC"],
    net_oxidative_balance = process[, "net_oxidative_balance"]
  )
  # latent_truth v3.1: prescribed simulator target — NO plant observables.
  # Formula: z(t) = 0.26*Bdir + 0.24*(Cacc/ΣQ) + 0.18*sqrt(α_acc*α_don) +
  #                 0.16*sqrt(k*_acc * k*_don) + 0.16*Mbalance
  # Bdir: bounded [0,1] transform of event-direction-specific accessible-capacity
  #   balance (Cacc_EAC - Cacc_EDC). Reference limits frozen at ±10 mmol e⁻ kg⁻¹.
  #   These are fixed illustrative limits; values outside them are truncated,
  #   so users must inspect saturation in every simulation design.
  # sqrt(α_acc * α_don): geometric mean of accept/donate connectivity (both ∈ [0,1]).
  # k* = k/0.20 h⁻¹ truncated to [0,1]; sqrt(k*_acc * k*_don): geometric mean.
  # Mbalance = 1 - memory: lower accumulated memory → better resilience.
  # Version 3.1 removes functional_readiness leakage and capacity_position from
  # the target. Weights sum to 1.0.
  Bdir_ref_lo <- -10   # frozen reference limit (mmol e⁻ kg⁻¹)
  Bdir_ref_hi <-  10   # frozen reference limit (mmol e⁻ kg⁻¹)
  Bdir <- clamp(
    (process[, "net_oxidative_balance"] - Bdir_ref_lo) /
      (Bdir_ref_hi - Bdir_ref_lo)
  )
  alpha_gm <- sqrt(pmax(process[, "alpha_accept"] * process[, "alpha_donate"], 0))
  k_gm     <- sqrt(
    clamp(process[, "k_accept"] / 0.20) *
    clamp(process[, "k_donate"] / 0.20)
  )
  Mbalance <- clamp(1 - state[, "memory"])
  latent_truth <- clamp(
    0.26 * Bdir +
    0.24 * process[, "Cacc_fraction"] +
    0.18 * alpha_gm +
    0.16 * k_gm +
    0.16 * Mbalance
  )

  conservation_checks <- data.frame(
    check = c(
      "maximum_absolute_Fe_mass_balance_error",
      "maximum_absolute_Mn_mass_balance_error",
      "minimum_simulated_pool", "nonnegative_pool_check"
    ),
    value = c(
      max(abs(fe_error_all), na.rm = TRUE),
      max(abs(mn_error_all), na.rm = TRUE),
      min(state[, setdiff(state_names, c(
        "memory", "root_biomass",
        "aerenchyma"
      ))], na.rm = TRUE),
      as.numeric(all(state[, setdiff(state_names, c(
        "memory", "root_biomass",
        "aerenchyma"
      ))] >= -1e-10))
    )
  )

  graph <- NULL
  if (isTRUE(include_graph) && requireNamespace("igraph", quietly = TRUE)) {
    graph <- if (isTRUE(stochastic_reassembly)) {
      igraph::sample_pa(max(20, p_micro), m = 2, directed = FALSE)
    } else {
      igraph::sample_smallworld(1, max(20, p_micro), 4, 0.08)
    }
    graph <- igraph::simplify(graph)
  }

  list(
    id = id,
    forcing = forcing,
    latent_state = latent_state,
    soil_data = soil_data,
    plant_data = plant_data,
    micro_gene_abundance = micro_gene_abundance,
    micro_metat_counts = micro_metat_counts,
    micro_metat_metadata = micro_metat_metadata,
    micro_traits = micro_traits,
    fluxes = fluxes,
    conservation_checks = conservation_checks,
    ROS_flux = plant_data,
    Eh_stability = Eh_stability,
    micro_data = micro_data,
    latent_truth = latent_truth,
    graph = graph,
    metadata = list(
      simulator_version = "3.1",
      seed = seed,
      RNGkind = RNGkind(),
      time_unit = "one day per transition",
      model_status = "Illustrative transitions; only Fe and Mn balances checked; alpha enters selected rates and k_accept gates crystalline-Fe reduction",
      intervention_scope = "Capacity changes initial Fe(III). Connectivity changes accessibility terms in selected Fe/Mn/N/S/C rates. Kinetics changes calculated accessible capacity; k_accept also gates crystalline-Fe reduction. Evaluate outcomes rather than assuming rescue.",
      validation_scope = "latent_truth v3.1 is a prescribed synthetic target: z(t) = 0.26*Bdir + 0.24*(Cacc/sumQ) + 0.18*sqrt(alpha_accept*alpha_donate) + 0.16*sqrt(k_accept_norm*k_donate_norm) + 0.16*Mbalance. No plant observable enters the target. Agreement is an internal diagnostic, not empirical or predictive validation.",
      hidden_input_columns = c("Cacc_EAC", "Cacc_EDC", "Cacc_total", "Cacc_fraction",
        "net_oxidative_balance", "alpha_accept", "alpha_donate", "k_accept_h", "k_donate_h"),
      scenario = scenario,
      rescue = rescue,
      n_cycles = n_cycles,
      event_centers = centers,
      event_tau_h = event_tau_h,
      Fe_Mn_units = "mmol kg-1 dry soil",
      EAC_EDC_units = "mmol electron equivalents kg-1 dry soil",
      gene_units = "synthetic DNA count scale labelled copies g-1 dry soil; not assay-calibrated",
      MetaT_units = "negative-binomial transcript counts",
      RNA_DNA_ratio_scope = "Legacy RNA_DNA_activity_ratio is a raw total-count ratio confounded by the separate observation scales; it is not an activity estimate",
      caution = paste(
        "Synthetic data test computational recovery and falsifiable model",
        "predictions; they do not constitute empirical ecological validation."
      )
    )
  )
}
