# Simulate illustrative soil-plant-microbe redox trajectories

Synthetic daily trajectories with explicit Fe/Mn redistribution, plant
indicators, gene abundance and transcript/count observation models. Only
Fe and Mn inventories have closed-balance checks. C, N, S and oxygen
budgets are not fully balanced. Model parameters are illustrative, not
fitted. Gene abundance indicates potential; transcript counts are
observations, not flux.

## Usage

``` r
simulate_redox_holobiont(
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
  sequencing_depth = 2e+05,
  metat_depth = 5e+05,
  decoupling = 0.25,
  zero_inflation = 0.2,
  MNAR_strength = 0.3,
  Eh_dropout_threshold = 100,
  micro_mean = 8,
  micro_slope = 3,
  micro_lambda_min = 1e-08,
  micro_lambda_max = 1e+06,
  stochastic_reassembly = TRUE,
  include_graph = FALSE,
  depth_labels = NULL
)
```

## Arguments

- n_plot:

  Positive integer. Number of plots (spatial replicates).

- n_depth:

  Positive integer. Number of depth strata per plot.

- n_plant:

  Positive integer. Number of plants per plot-depth unit. Must be \\\geq
  1\\. Plants contribute to ROL and ROS signals.

- n_time:

  Positive integer \\\geq 4\\. Number of daily time steps.

- p_micro:

  Positive integer. Number of ASV-like microbial taxonomic features
  generated alongside functional gene data.

- seed:

  Integer or `NULL`. Random seed passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) before simulation;
  `NULL` means no seeding (non-reproducible). All manuscript figures use
  explicit seeds.

- scenario:

  Character; one of `"flood_drain"` (default) or `"drought_rewet"`.
  Determines the shape of the hydrological forcing function and the sign
  of the dominant redox transition.

- n_cycles:

  Positive integer. Number of forcing pulses.

- disturbance_strength:

  Numeric in \\\[0, 1\]\\. Event severity. Controls peak WFPS, anaerobic
  volume fraction, and the amplitude of soil redox transitions.

- disturbance_center:

  Numeric or `NULL`. Time step of the first disturbance event centre.
  Defaults to evenly spaced centres from 0.22 to 0.78 of n_time.

- disturbance_width:

  Numeric in \\(0, 1)\\. Width of each Gaussian forcing pulse as a
  fraction of `n_time`.

- seasonal_amp:

  Numeric \\\geq 0\\. Amplitude of additive seasonal forcing overlaid on
  the hydrological disturbance signal.

- seasonal_phase:

  Numeric. Phase offset (radians) of the seasonal forcing.

- history_strength:

  Numeric in \\\[0, 1\]\\. Scales the synthetic memory state at
  initialization and during disturbance. That state affects Fe
  crystallisation, accessibility and generated microbial descriptors; it
  is not a measured fraction of community carry-over.

- rescue:

  Character; one of `"none"` (default), `"capacity"`, `"connectivity"`,
  or `"kinetics"`. Simulates a targeted scenario modification (capacity
  at initialization; alpha/k throughout):

  `"capacity"`

  :   Fe(III) inventory replenishment.

  `"connectivity"`

  :   Increases alpha_accept and alpha_donate in selected Fe/Mn/N/S/C
      rate expressions and calculated accessible capacity. Effects on
      any recovery outcome must be evaluated, not assumed beneficial.

  `"kinetics"`

  :   Increases both calculated accessible-capacity rates. k_accept
      additionally gates crystalline-Fe reduction; k_donate has no
      direct process-rate gate. This intervention is not a general
      exchange-kinetics model.

- event_tau_h:

  Positive numeric. Disturbance timescale \\\tau\\ (h) used internally
  by the accessible-capacity calculation; passed to
  [`rri_accessible_capacity`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md).

- sequencing_depth:

  Positive numeric. Mean library size for taxonomic count data
  (`micro_data` block), modelled as a negative-binomial process.

- metat_depth:

  Positive numeric. Mean library size for metatranscript counts
  (`micro_metat_counts`), modelled separately from taxonomic counts with
  a higher biological variance.

- decoupling:

  Numeric in \\\[0, 1\]\\. Cross-domain stochastic decoupling parameter
  for selected noise terms; other stochastic terms remain.

- zero_inflation:

  Numeric in \\\[0, 1\]\\. Structural-zero probability for taxonomic
  count features (simulates taxa absent from some samples).

- MNAR_strength:

  Numeric in \\\[0, 1\]\\. Maximum missing-not-at-random (MNAR)
  probability for Eh values under strongly reducing conditions.

- Eh_dropout_threshold:

  Numeric. Eh (mV) below which the MNAR dropout probability begins to
  rise; is an artificial missingness design, not a platinum-electrode
  detection limit.

- micro_mean, micro_slope, micro_lambda_min, micro_lambda_max:

  Backward-compatible parameters controlling mean and slope of the
  log-linear model for taxonomic count intensity. See legacy
  documentation.

- stochastic_reassembly:

  Logical. If `TRUE` (default), adds additional stochastic variation;
  this is not an explicit succession model.

- include_graph:

  Logical. If `TRUE` and igraph is installed, returns an `igraph` random
  graph object independent of the generated community in `$graph`.

- depth_labels:

  Character vector of length `n_depth` or `NULL`. Custom labels for
  depth strata. Defaults to `"D1"`, `"D2"`, ...

## Value

List of identifiers, data blocks, latent states, flux descriptors,
balance checks, metadata and legacy views. The metadata records hidden
columns, seed and RNG configuration. With seed supplied, the caller RNG
is restored.

## Details

alpha and k are prescribed internal state variables, not inferred
parameters: they are computed from forcing, pore structure and memory at
the start of each time step. Alpha multiplies selected Fe/Mn/N/S/C
expressions as an accessibility factor, whereas k_accept directly gates
only the model's crystalline-Fe reduction expression. A low alpha
suppresses selected rates; restoring alpha changes selected rates, but
does not guarantee functional recovery. k additionally gates the
crystalline-Fe exchange rate. The latent_truth vector is a constructed
index sharing ingredients with soil_data; omit those ingredients from
observable-only benchmarks. history_pair denotes a shared random effect,
not experimentally matched twins. Low n_time can leave no adequate
baseline or recovery; inspect forcing and analyse each event separately.
Mineral crystallisation is continuous in this implementation, not an
event-only reoxidation ratchet. k does not explicitly decrease with
crystallinity. ROL, ROS and several microbial features are synthetic
descriptors, not calibrated physical fluxes. Do not infer quantitative
field rates from their labels alone. Q_accept now counts one electron
per crystalline Fe(III); the previous 0.22 factor mixed accessibility
into inventory. DOC reducing equivalents remain an illustrative
coefficient dependent on assumed carbon oxidation state.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim <- simulate_redox_holobiont(n_plot = 2, n_depth = 2, n_plant = 2,
                                   n_time = 20, seed = 42)
  nrow(sim$id)  # 2 x 2 x 2 x 20 = 160 rows
  names(sim)    # top-level list elements
} # }
```
