# Changelog

## HRRI 0.99.1

### Documentation and infrastructure

- Added
  [`attach_hrri_ids()`](https://mghotbi.github.io/HRRI/reference/attach_hrri_ids.md):
  joins design identifiers onto a score table using a unique `row_id` or
  a complete observation key, falling back to row order only with an
  explicit warning. Records the method used in the `id_alignment`
  attribute.
- Completed roxygen documentation: every exported function now carries
  `@importFrom` declarations for all external calls and a worked
  `@examples` block. Regenerated `NAMESPACE` and all help pages.
- Added a styled HTML vignette theme and a pkgdown site configuration
  with the function reference grouped by topic.
- Added GitHub Actions workflows for multi-platform `R CMD check` and
  pkgdown deployment.
- Fixed a parse error in `test-property-scores.R` and corrected two test
  expectations that compared per-row `domain_coverage` and `n_domains`
  columns against scalars.

### Scientific corrections

- Reframed the simulator benchmark as diagnostic agreement with a
  prescribed synthetic target, not empirical validation, prediction or
  interval coverage.
- Reframed
  [`hrri_infer_architecture()`](https://mghotbi.github.io/HRRI/reference/hrri_infer_architecture.md)
  as a conditional recovery-curve fit. Its fitted amplitude and rate are
  no longer described as mechanistic accessibility and reservoir
  exchange kinetics.
- Corrected oxygen-demand documentation: the calculation is a
  stoichiometric potential demand and the comparison is to a specified
  oxygen stock, not air-saturated porewater or a thermodynamic recovery
  ceiling.
- Clarified simulator scope: only Fe and Mn have closed-balance checks;
  parameters and synthetic flux-like descriptors are not field
  calibrated.
- Included Mn(III) in the simulator electron-donating inventory.

### Interface and reliability

- [`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  now forwards domain weights and accepts partial soil, plant or
  microbial panels.
- Missing domains remain missing. Per-row coverage, observed-domain
  count and effective weights are returned.
- Known hidden simulator states are excluded from observation-only
  scoring.
- Benchmark corruption is reproducible by seed, failed seeds are
  reported, and score rows are explicitly aligned to identifiers.
- Fixed the single-row dimension drop in
  [`rri_memory_index()`](https://mghotbi.github.io/HRRI/reference/rri_memory_index.md).
- Added strict numeric, duplicate-name, key-alignment and graph checks.
- Made the illustrative reservoir template default to the core bulk
  EAC/EDC schema; phase-resolved reservoirs now require explicit column
  names.
- Added explicit +1/-1 orientation for root-trait composites and reject
  reuse of one measurement as multiple weighted components.
- Changed the time-series plot to separate panels for unlike units.
- Corrected validation-plot RMSE to compare the score directly with its
  target.

### Packaging

- Removed generated reports, figures, macOS metadata and author-only
  scripts from the installable source tree; the original upload retains
  those files.
- Removed unreachable optional dependencies and corrected DESCRIPTION
  metadata.
- Added submission audit and validation instructions.
