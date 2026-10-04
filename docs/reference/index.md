# Package index

## Simulation and bundled example inputs

Generate illustrative soil–plant–microbe trajectories and the declared
reservoir template used by the demonstrations.

- [`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
  : Simulate illustrative soil-plant-microbe redox trajectories
- [`rri_simulation_demo()`](https://mghotbi.github.io/HRRI/reference/rri_simulation_demo.md)
  : Reproducible observable-only HRRI demonstration
- [`rri_default_reservoirs()`](https://mghotbi.github.io/HRRI/reference/rri_default_reservoirs.md)
  : Illustrative reservoir parameter template
- [`HRRI_workflow_example`](https://mghotbi.github.io/HRRI/reference/HRRI_workflow_example.md)
  [`redoxrri_example`](https://mghotbi.github.io/HRRI/reference/HRRI_workflow_example.md)
  : Worked example: a complete HRRI workflow

## Mechanistic calculation from declared inputs

Quantities computed from inventories, accessibility and exchange rates
that the user supplies. These are never estimated from the scores.

- [`rri_accessible_capacity()`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md)
  : Event-window accessible reservoir capacities
- [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  : Stoichiometric O\\\_2\\ Demand from Reduced-Pool Inventories

## Domain scoring from aligned observations

Build oriented soil, plant and microbial scores, with coverage and
missingness retained rather than recoded as zero.

- [`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  : Score observed soil, plant and microbial panels
- [`rri_pipeline_st()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md)
  : Exploratory domain-score integration (legacy interface)
- [`attach_hrri_ids()`](https://mghotbi.github.io/HRRI/reference/attach_hrri_ids.md)
  : Attach design identifiers to a score table with explicit alignment
  checks
- [`rri_reference_scores()`](https://mghotbi.github.io/HRRI/reference/rri_reference_scores.md)
  : Score departures from an explicitly defined reference
- [`rri_capacity_index()`](https://mghotbi.github.io/HRRI/reference/rri_capacity_index.md)
  : Oxidative-oriented soil feature composite
- [`rri_root_physio()`](https://mghotbi.github.io/HRRI/reference/rri_root_physio.md)
  : Exploratory root-trait composite
- [`rri_micro_functional_score()`](https://mghotbi.github.io/HRRI/reference/rri_micro_functional_score.md)
  : Construct an explicitly weighted microbial guild contrast
- [`rri_kinetics_score()`](https://mghotbi.github.io/HRRI/reference/rri_kinetics_score.md)
  : Descriptive recovery speed score
- [`rri_connectivity_score()`](https://mghotbi.github.io/HRRI/reference/rri_connectivity_score.md)
  : Cross-domain association or graph-topology summary
- [`rri_compensation_index()`](https://mghotbi.github.io/HRRI/reference/rri_compensation_index.md)
  : Cross-domain asynchrony diagnostic

## Recovery and memory descriptors

Descriptive diagnostics for a single disturbance, reported with fit
status and estimability so an unsupported quantity stays missing.

- [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md)
  : Descriptive recovery metrics for a single disturbance
- [`rri_memory_index()`](https://mghotbi.github.io/HRRI/reference/rri_memory_index.md)
  : Persistent-displacement and loop-area diagnostic

## Agreement, sensitivity and benchmarking

Compare a score with a reference target while respecting clustering, and
probe how the result depends on weights and specification.

- [`rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/rri_accuracy.md)
  [`print(`*`<rri_accuracy>`*`)`](https://mghotbi.github.io/HRRI/reference/rri_accuracy.md)
  : Agreement between a score and a reference target, respecting
  clustering
- [`rri_latent_correlation()`](https://mghotbi.github.io/HRRI/reference/rri_latent_correlation.md)
  : Correlation with a Simulator-Defined Target
- [`rri_property_scores()`](https://mghotbi.github.io/HRRI/reference/rri_property_scores.md)
  : Summarise supported diagnostics without fabricating missing
  properties
- [`rri_domain_influence()`](https://mghotbi.github.io/HRRI/reference/rri_domain_influence.md)
  : Realised influence of each domain on the composite score
- [`rri_sensitivity()`](https://mghotbi.github.io/HRRI/reference/rri_sensitivity.md)
  : Sensitivity to domain aggregation weights
- [`benchmark_hrri()`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
  [`print(`*`<hrri_benchmark>`*`)`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
  : Benchmark diagnostic agreement with a simulator-defined target
- [`hrri_infer_architecture()`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`print(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`summary(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`as.data.frame(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`validate_architecture()`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`print(`*`<hrri_arch_validation>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  : Fit a conditional capacity-recovery curve (legacy function name)

## Figures

Publication figures used in the manuscript, plus the general-purpose
diagnostic plots and the shared theme.

- [`plot_rri_framework()`](https://mghotbi.github.io/HRRI/reference/plot_rri_framework.md)
  : Complementary mechanistic and observational routes
- [`plot_rri_identifiability()`](https://mghotbi.github.io/HRRI/reference/plot_rri_identifiability.md)
  : Identifiability, sensitivity and oxygen-accounting illustrations
- [`plot_rri_timeseries()`](https://mghotbi.github.io/HRRI/reference/plot_rri_timeseries.md)
  : Aligned time series with separate physical units
- [`plot_rri_recovery_diagnostics()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_diagnostics.md)
  : Aligned score trajectories and recovery-descriptor availability
- [`plot_rri_properties()`](https://mghotbi.github.io/HRRI/reference/plot_rri_properties.md)
  : Profile of Available HRRI Diagnostics
- [`plot_rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/plot_rri_accuracy.md)
  : Four-panel diagnostic figure for an accuracy assessment
- [`plot_rri_recovery_map()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_map.md)
  : Plot RRI Recovery Map
- [`plot_rri_recovery_landscape()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_landscape.md)
  : Plot a recovery landscape from RRI perturbation-recovery metrics
- [`plot_rri_state_space()`](https://mghotbi.github.io/HRRI/reference/plot_rri_state_space.md)
  : Plot domain-score space with correctly matched trajectory
  diagnostics
- [`plot_rri_validation()`](https://mghotbi.github.io/HRRI/reference/plot_rri_validation.md)
  : Scatter of HRRI score against a simulator-defined target
- [`plot_RRI_ternary()`](https://mghotbi.github.io/HRRI/reference/plot_RRI_ternary.md)
  : Ternary Plot of Relative Domain Scores
- [`plot_hrri_benchmark()`](https://mghotbi.github.io/HRRI/reference/plot_hrri_benchmark.md)
  : Plot descriptive benchmark agreement
- [`plot_rri_simulation_demo()`](https://mghotbi.github.io/HRRI/reference/plot_rri_simulation_demo.md)
  : Draw the reproducible simulation demonstration
- [`theme_ems()`](https://mghotbi.github.io/HRRI/reference/theme_ems.md)
  : EMS plotting theme
