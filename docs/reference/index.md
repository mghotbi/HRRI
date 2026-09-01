# Package index

## Simulation

Generate illustrative multi-domain trajectories with closed Fe and Mn
inventories under flood–drain or drought–rewet forcing.

- [`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
  : Simulate illustrative soil-plant-microbe redox trajectories
- [`rri_simulation_demo()`](https://mghotbi.github.io/HRRI/reference/rri_simulation_demo.md)
  : Reproducible observable-only HRRI demonstration
- [`plot_rri_simulation_demo()`](https://mghotbi.github.io/HRRI/reference/plot_rri_simulation_demo.md)
  : Draw the reproducible simulation demonstration

## Capacity and stoichiometry

Event-window accessible capacity and electron-balance calculations from
explicitly supplied inventories and kinetic parameters.

- [`rri_accessible_capacity()`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md)
  : Event-window accessible reservoir capacities
- [`rri_default_reservoirs()`](https://mghotbi.github.io/HRRI/reference/rri_default_reservoirs.md)
  : Illustrative reservoir parameter template
- [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  : Stoichiometric O\\\_2\\ Demand from Reduced-Pool Inventories
- [`rri_capacity_index()`](https://mghotbi.github.io/HRRI/reference/rri_capacity_index.md)
  : Oxidative-oriented soil feature composite

## Domain scoring

Integrate plant, soil and microbial observations into composite scores,
with identifier alignment and coverage diagnostics.

- [`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  : Score observed soil, plant and microbial panels
- [`rri_pipeline_st()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md)
  : Exploratory domain-score integration (legacy interface)
- [`rri_reference_scores()`](https://mghotbi.github.io/HRRI/reference/rri_reference_scores.md)
  : Score departures from an explicitly defined reference
- [`attach_hrri_ids()`](https://mghotbi.github.io/HRRI/reference/attach_hrri_ids.md)
  : Attach design identifiers to a score table with explicit alignment
  checks
- [`rri_micro_functional_score()`](https://mghotbi.github.io/HRRI/reference/rri_micro_functional_score.md)
  : Construct an explicitly weighted microbial guild contrast
- [`rri_root_physio()`](https://mghotbi.github.io/HRRI/reference/rri_root_physio.md)
  : Exploratory root-trait composite

## Recovery signatures

Describe decline and return of a higher-is-better score across a single
disturbance event.

- [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md)
  : Descriptive recovery metrics for a single disturbance
- [`rri_memory_index()`](https://mghotbi.github.io/HRRI/reference/rri_memory_index.md)
  : Persistent-displacement and loop-area diagnostic
- [`rri_kinetics_score()`](https://mghotbi.github.io/HRRI/reference/rri_kinetics_score.md)
  : Descriptive recovery speed score
- [`rri_connectivity_score()`](https://mghotbi.github.io/HRRI/reference/rri_connectivity_score.md)
  : Cross-domain association or graph-topology summary
- [`rri_compensation_index()`](https://mghotbi.github.io/HRRI/reference/rri_compensation_index.md)
  : Cross-domain asynchrony diagnostic
- [`rri_property_scores()`](https://mghotbi.github.io/HRRI/reference/rri_property_scores.md)
  : Summarise supported diagnostics without fabricating missing
  properties
- [`rri_sensitivity()`](https://mghotbi.github.io/HRRI/reference/rri_sensitivity.md)
  : Sensitivity to domain aggregation weights

## Benchmarking and inference

Internal consistency checks against prescribed synthetic targets. These
are not empirical validation.

- [`benchmark_hrri()`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
  [`print(`*`<hrri_benchmark>`*`)`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
  : Benchmark diagnostic agreement with a simulator-defined target
- [`plot_hrri_benchmark()`](https://mghotbi.github.io/HRRI/reference/plot_hrri_benchmark.md)
  : Plot descriptive benchmark agreement
- [`hrri_infer_architecture()`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`print(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`summary(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`as.data.frame(`*`<hrri_arch>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`validate_architecture()`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  [`print(`*`<hrri_arch_validation>`*`)`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  : Fit a conditional capacity-recovery curve (legacy function name)
- [`rri_latent_correlation()`](https://mghotbi.github.io/HRRI/reference/rri_latent_correlation.md)
  : Correlation with a Simulator-Defined Target

## Visualisation

Publication-oriented diagnostic plots.

- [`plot_rri_timeseries()`](https://mghotbi.github.io/HRRI/reference/plot_rri_timeseries.md)
  : Aligned time series with separate physical units
- [`plot_rri_state_space()`](https://mghotbi.github.io/HRRI/reference/plot_rri_state_space.md)
  : Plot domain-score space with correctly matched trajectory
  diagnostics
- [`plot_RRI_ternary()`](https://mghotbi.github.io/HRRI/reference/plot_RRI_ternary.md)
  : Ternary Plot of Relative Domain Scores
- [`plot_rri_recovery_map()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_map.md)
  : Plot RRI Recovery Map
- [`plot_rri_recovery_landscape()`](https://mghotbi.github.io/HRRI/reference/plot_rri_recovery_landscape.md)
  : Plot a recovery landscape from RRI perturbation-recovery metrics
- [`plot_rri_properties()`](https://mghotbi.github.io/HRRI/reference/plot_rri_properties.md)
  : Radar Chart of Available HRRI Diagnostics
- [`plot_rri_validation()`](https://mghotbi.github.io/HRRI/reference/plot_rri_validation.md)
  : Scatter of HRRI score against a simulator-defined target
- [`theme_ems()`](https://mghotbi.github.io/HRRI/reference/theme_ems.md)
  : EMS plotting theme

## Data

Example dataset shipped with the package.

- [`redoxrri_example`](https://mghotbi.github.io/HRRI/reference/redoxrri_example.md)
  : Example workflow for RedoxRRI
