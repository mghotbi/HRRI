# Changelog

## HRRI 1.0.2

### Fixes

- `tests/testthat/test-api-surface.R` required a `\value` section on
  every Rd topic. That is wrong for a narrative page:
  `HRRI_workflow_example` has no `\usage` and returns nothing, so the
  test failed on a correctly written topic. `\value` is now required
  only of function topics, which is the rule `R CMD check` itself
  applies.

### Note for anyone upgrading over an existing source tree

`unzip -o` overwrites files but never deletes them. Release 1.0.0
renamed `R/data-redoxrri_example.R` to `R/HRRI_workflow_example.R`, so a
tree updated by extracting over the old one keeps both, and both `.Rd`
files claim `\alias{redoxrri_example}`. Duplicate aliases make
`R CMD INSTALL` fail, and
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
then renders vignettes against whatever older copy of the package is
already installed, producing a confusing “could not find function” for a
function that is plainly present in the sources.

`sync_package.sh`, distributed alongside the release, compares a tree
against `package_manifest.txt` and reports or removes the leftovers.

## HRRI 1.0.1

### Fixes

- [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  no longer fails obscurely when every reduced-species column the caller
  names is absent from `soil_df`. The numeric coercion ran before the
  validity check, and with no matching columns it attempted a
  zero-column replacement, so `[<-.data.frame` raised a cryptic error
  instead of the informative “none of the specified reduced-species
  columns were found” message a few lines below. The coercion is now
  guarded.

### Documentation

- [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  cited “Ghotbi et al. (2026)” in its details and example while the
  reference entry carries no year, the mechanistic review being
  unpublished. The in-text citations now read “submitted”, matching the
  reference.

## HRRI 1.0.0

First release candidate. No change to any calculation; this release
completes the packaging, documentation and infrastructure.

### New infrastructure

- `_pkgdown.yml` rebuilt from the actual Rd topics and grouped by task
  rather than alphabetically: scoring pipelines, the four hidden states,
  recovery signatures, model checking, figures, simulation. All 36
  topics are listed exactly once, which is what the previous index
  failed at when
  [`rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/rri_accuracy.md),
  [`plot_rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/plot_rri_accuracy.md)
  and
  [`rri_domain_influence()`](https://mghotbi.github.io/HRRI/reference/rri_domain_influence.md)
  were added without entries and the site build stopped.

- GitHub Actions workflows for `R CMD check` (five platform and version
  combinations, run with `--as-cran --run-donttest`) and for the pkgdown
  site. Vignettes are built in CI, so a repeat of the ggtern
  incompatibility that once broke vignette building would surface before
  release.

- `cran-comments.md`, `.gitignore`, and a `.Rbuildignore` that keeps all
  of the above out of the source tarball.

- `tests/testthat/test-api-surface.R` checks the public interface
  itself: every export resolves, S3 methods are registered rather than
  merely defined, argument names of the entry points are pinned, every
  topic has a value and a runnable example, no example is hidden behind
  `\dontrun`, and seeding leaves the caller’s `.Random.seed` exactly as
  found. These guard the failure modes that `R CMD check` cannot see.

### Documentation

- [`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
  keeps `seed = 123` rather than `NULL`, so that an unseeded call to an
  illustrative generator is reproducible. The `@param` entry now says
  so, and explains why it is safe: the RNG kind and `.Random.seed` are
  saved on entry and restored on exit, so the caller’s stream is
  untouched. The diagnostic functions default to `seed = NULL`.

- The topic `redoxrri_example` is renamed `HRRI_workflow_example`. Its
  former name implied a dataset, but the package ships no `data/`
  directory and `data(redoxrri_example)` would have failed; the page now
  states plainly that it is a worked example rather than an object. The
  old name is kept as an alias, so existing links still resolve.

- Remaining references to the package’s former name, “RedoxRRI”,
  replaced with “HRRI” throughout the help pages.

- `inst/CITATION` now matches the current software manuscript: four
  authors and the title “HRRI: direction-aware diagnostics for
  soil-plant-microbiome redox recovery across hydroclimatic
  disturbances”, consistent with README and the vignettes.

## HRRI 0.99.12

### Documentation

- Manuscript statuses corrected throughout. The mechanistic review is
  submitted to Soil Biology & Biochemistry, and the theory paper to
  Communications Earth & Environment; the software paper remains in
  preparation. Applied in README, the workflow vignette, `inst/CITATION`
  and the
  [`rri_o2_demand()`](https://mghotbi.github.io/HRRI/reference/rri_o2_demand.md)
  help page, which cites the review for its stoichiometric table.

- `inst/CITATION` had the software paper as “Manuscript submitted” while
  README and the vignette said “in preparation”. Aligned to in
  preparation, matching the wording confirmed for README. A single line
  in `inst/CITATION` reverses this if the paper has in fact been
  submitted.

## HRRI 0.99.11

### Fixes

- Added `.Rbuildignore`, clearing both remaining `R CMD check` NOTEs:
  `.github` under “checking for hidden files and directories”, and
  `LICENSE.md` and `_pkgdown.yml` under “checking top-level files”.
  These are development and website files that belong in the repository
  but not in a source tarball.

  Each directory is listed in two forms, `^docs$` and `^docs/`. A
  pattern anchored with a trailing `$` matches the directory entry but
  not the files inside it, which is why an earlier `^\.github$` left the
  hidden-files NOTE in place. The patterns were checked against a list
  of representative paths before release: every development file is
  excluded and all fourteen sampled package files survive.

## HRRI 0.99.10

### Documentation

- The theory manuscript is recorded with its confirmed author list and
  status: Ghotbi, Kolody, Ghotbi and Holtgrewe-Stukenbrock, “A Theory of
  Hydroclimatic Redox Resilience”, submitted to Communications Earth &
  Environment. Updated in `inst/CITATION`, README and the workflow
  vignette.

### Known inconsistency

- `inst/CITATION` still lists the package manuscript under its earlier
  title with five authors including Marco Guerreiro, while README and
  the workflow vignette follow the current draft: four authors, and the
  title “HRRI: direction-aware diagnostics for soil-plant-microbiome
  redox recovery across hydroclimatic disturbances”. The author list is
  not something to change by inference, so both forms are left in place
  until confirmed. The review manuscript’s status is likewise recorded
  as in preparation because it has not been confirmed.

## HRRI 0.99.9

### Documentation fixes in `rri_o2_demand()`

- The interpretation section was headed “the 25-fold contrast” while its
  own text, and the example, said 26. The worked example is 33.75 mmol
  O2 per kg against a 1.3 mmol per kg stock, a ratio of 25.96, so 26 is
  right and the heading was wrong.

- `@return` did not describe what the function returns. Five elements
  were undocumented (`n_species_observed`, `species_coverage`,
  `interpretation`, `ch4_unit_used`, `acetate_basis_used`), and three
  columns of `components` were listed under names the code does not use:
  `mean_inventory` is really `mean_inventory_mmol`,
  `fraction_of_total_demand` is really `fraction_total_demand`, and
  `n_observed` was missing. `stoich_table` was described as “the full
  stoichiometric table” when it holds only the species actually found.
  All eleven returned elements are now documented under their real
  names. `R CMD check` cannot detect this class of error.

- The reference to the mechanistic review was updated to the current
  draft: four authors, and “disturbances” rather than “disturbance”.

- Verified DOIs added to the two rate-law references, Stumm and
  Lee (1961) <doi:10.1021/ie50614a030> and Millero, Sotolongo and
  Izaguirre (1987) <doi:10.1016/0016-7037(87)90093-7>.

### Verification

The stoichiometry was checked against the electron balance rather than
taken on trust. All eight coefficients are correct at 4 electrons per
O2: Fe(II) 1 e- (0.25), Mn(II) 2 e- (0.50), HS- 8 e- (2.00), FeS 9 e-
(2.25), FeS2 15 e- (3.75), NH4+ 8 e- (2.00), CH4 8 e- (2.00), acetate 8
e- (2.00) and 1.00 on a carbon basis. The pyrite equation balances on
Fe, S, O and H, and releases the 4 mol H+ per mol FeS2 the text claims.

## HRRI 0.99.8

### Documentation

- Every documented function now carries an explicit `@title`.
  Twenty-nine roxygen blocks had been relying on the implicit first-line
  title; six S3 methods documented with `@rdname` or `@describeIn`
  correctly do not have one, since they inherit the parent topic’s title
  and a second `@title` would conflict.

- `DESCRIPTION` cites six published methods references in the form CRAN
  asks for, `authors (year) <doi:...>` with no space after `doi:`. Every
  DOI was resolved against Crossref before being added; none is a
  placeholder.

  - Sander, Hofstetter and Gorski (2015) <doi:10.1021/acs.est.5b00006>,
    mediated electrochemical determination of EAC and EDC.
  - Kluepfel, Piepenbrock, Kappler and Sander (2014)
    <doi:10.1038/ngeo2084>, regeneration of electron-accepting capacity
    across anoxic periods.
  - Thompson, Chadwick, Rancourt and Chorover (2006)
    <doi:10.1016/j.gca.2005.12.005>, iron-oxide crystallinity increases
    under redox oscillation, the empirical basis for the mineralogical
    ratchet.
  - Keiluweit, Wanzek, Kleber, Nico and Fendorf (2017)
    <doi:10.1038/s41467-017-01406-6>, anaerobic microsites in aerobic
    soil.
  - Lin (1989) <doi:10.2307/2532051> and Kobayashi and Salam (2000)
    <doi:10.2134/agronj2000.922345x>, the agreement statistics.

- The three companion manuscripts are described in README, the workflow
  vignette and `inst/CITATION`, each marked in preparation, with a note
  that entries are provisional until publication. They are deliberately
  absent from `DESCRIPTION`: CRAN asks that field to carry only
  references a reader can retrieve, so unpublished work belongs in these
  files instead.

- README gained a table mapping each published reference to the specific
  part of HRRI that rests on it, so a reader can see why each is cited
  rather than taking the list on trust.

## HRRI 0.99.6

### Fixes

- [`rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/rri_accuracy.md)
  no longer modifies the user’s global environment.
  [`set.seed()`](https://rdrr.io/r/base/Random.html) writes
  `.Random.seed` into `.GlobalEnv`, so calling it without restoring the
  previous state silently resets a random stream the caller may have
  been relying on. The RNG kind and seed are now saved on entry and put
  back with [`on.exit()`](https://rdrr.io/r/base/on.exit.html), matching
  the pattern
  [`simulate_redox_holobiont()`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md)
  and
  [`benchmark_hrri()`](https://mghotbi.github.io/HRRI/reference/benchmark_hrri.md)
  have always used. This was the only place in the package that seeded
  without restoring.

### Documentation

- Audited against the CRAN Cookbook “Code Issues” chapter. Clean on
  every recipe: no `T`/`F` for `TRUE`/`FALSE`; every
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is behind a
  user-supplied argument defaulting to `NULL`; all
  [`cat()`](https://rdrr.io/r/base/cat.html)/[`print()`](https://rdrr.io/r/base/print.html)
  calls sit inside `print`/`summary` methods or behind a
  `verbose`/`print` argument;
  [`par()`](https://rdrr.io/r/graphics/par.html) is saved with
  `no.readonly = TRUE` and restored via
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html); nothing writes to
  the home filespace or the working directory; no
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html) detritus, no
  [`installed.packages()`](https://rdrr.io/r/utils/installed.packages.html),
  no `options(warn = -1)`, no parallelism, no software installation from
  package code.

- Both vignettes now restore the ggplot2 theme they set.
  [`theme_set()`](https://ggplot2.tidyverse.org/reference/get_theme.html)
  changes state that persists for the rest of the session; vignettes
  build in their own process so nothing outside was affected, but it is
  the same courtesy CRAN asks for with
  [`par()`](https://rdrr.io/r/graphics/par.html) and
  [`options()`](https://rdrr.io/r/base/options.html).

- All `\dontrun{}` example wrappers replaced with `\donttest{}`.
  `\dontrun` is for code that *cannot* run – code needing credentials,
  external software or hardware. None of these examples was in that
  category, and the wrapper meant `R CMD check` never executed them, so
  they could drift out of step with the code without anyone noticing.
  They are now run under `R CMD check --run-donttest` and by CRAN.

- The
  [`plot_RRI_ternary()`](https://mghotbi.github.io/HRRI/reference/plot_RRI_ternary.md)
  example is additionally guarded on the ggplot2 version. Unguarded it
  would error under ggplot2 \>= 4.0.0, where ggtern cannot be used; it
  now skips instead.

- `DESCRIPTION` cites the methodological reference and the two sources
  for the agreement statistics, as CRAN policy asks. The framework
  citation is a repository URL until the manuscript has a DOI, at which
  point it should be replaced with `<doi:...>`.

## HRRI 0.99.5

### Fixes

- Vignette building failed with

      Error in `plot_theme()`:
      ! The `tern.axis.ticks.length.major` theme element must be a <rel> object.

  The cause is **ggtern**, not the figure that happened to trigger it.
  ggtern patches ggplot2’s internal element tree when its namespace is
  *loaded*, not when a ternary plot is drawn. Under ggplot2 \>= 4.0.0,
  which rewrote the plot object on S7, that patch leaves elements in a
  form ggplot2 rejects, and every subsequent `ggplot` in the session
  fails in `plot_theme()`.

  The previous guard called
  [`requireNamespace("ggtern")`](https://rdrr.io/r/base/ns-load.html)
  and wrapped the plot in [`try()`](https://rdrr.io/r/base/try.html).
  That is too late:
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) loads the
  namespace, which is the damage, and
  [`try()`](https://rdrr.io/r/base/try.html) only catches the ternary
  plot’s own failure. Since `R CMD build` builds vignettes in one R
  session, the gallery vignette’s probe brought down the workflow
  vignette that followed it.

  The ggplot2 version is now checked *before* ggtern is touched, in a
  new internal helper `hrri_ggtern_ok()`, which never calls
  [`requireNamespace("ggtern")`](https://rdrr.io/r/base/ns-load.html) on
  an incompatible ggplot2. Applied in
  [`plot_RRI_ternary()`](https://mghotbi.github.io/HRRI/reference/plot_RRI_ternary.md),
  the gallery vignette, the `redoxrri_example` example and
  `tests/testthat/test-plot_rri_ternary.R` – the test had the same
  defect, because `skip_if_not_installed()` also loads the namespace and
  would have poisoned every later test file.

  Where ggtern cannot be used the ternary panel is skipped with an
  explanation, and the same composition is reported numerically, as
  before.

## HRRI 0.99.4

### New

- [`rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/rri_accuracy.md)
  assesses how closely a score tracks a reference target using
  statistics appropriate to repeated observations of the same units. A
  pooled row-wise correlation over a longitudinal panel is not a valid
  precision claim: 480 rows that are 12 trajectories sampled 40 times
  carry far less information than 480 independent observations. The
  function reports the intra-class correlation, the design effect
  `1 + (m - 1) * ICC`, and the effective sample size, and obtains
  intervals by resampling whole trajectories. It also reports the naive
  row-resampled interval alongside, so the cost of ignoring the
  clustering is shown rather than asserted.

  Agreement is reported as Lin’s concordance correlation coefficient
  next to Pearson’s `r`. The two differ whenever the score is
  miscalibrated: a score equal to twice the target plus a constant has
  `r = 1` and concordance near 0.3. Calibration intercept and slope are
  given so the direction of the miscalibration is visible. Mean squared
  error is partitioned into squared bias, variance mismatch and lack of
  correlation (Kobayashi and Salam 2000); the three components sum to
  MSE exactly, and the residual is returned as an attribute so a future
  regression in the arithmetic would be caught.

  The permutation null exchanges whole trajectories rather than rows, so
  each unit keeps its temporal shape. Permuting rows would destroy the
  shared event-driven signal and make almost any score appear
  significant.

- [`plot_rri_accuracy()`](https://mghotbi.github.io/HRRI/reference/plot_rri_accuracy.md)
  draws the assessment as a four-panel figure: calibration against the
  1:1 line, a Bland-Altman agreement plot with limits computed from
  cluster means, the row-resampled and trajectory-resampled bootstrap
  distributions side by side, and the exact error partition. Panels are
  selectable, and the figure degrades to a list of plots rather than
  failing when **patchwork** is absent.

### Fixes

- `man/rri_domain_influence.Rd` was missing from the 0.99.3 tarball: the
  function was exported and had a complete roxygen block, but no Rd file
  was generated, which is an “undocumented code objects” warning under
  `R CMD check`. Added.
- Removed `vignettes/HRRI.R`, an orphaned knitr output with no
  corresponding `.Rmd`, which `R CMD check` flags as a file in
  `vignettes/` that builds no vignette.

### Documentation

- The workflow vignette’s validation section is rewritten. It previously
  reported a single pooled correlation against `latent_truth`; it now
  reports cluster-aware agreement, shows the figure, and states plainly
  that a target and a score sharing a generator establish internal
  consistency of the estimator and not empirical validation.

## HRRI 0.99.3

### New

- [`rri_domain_influence()`](https://mghotbi.github.io/HRRI/reference/rri_domain_influence.md)
  reports the realised contribution of each domain to the composite, as
  an exact variance decomposition (phi_d = w_d Cov(S_d, R) / Var(R),
  summing to one). Declared weights and realised influence diverge
  whenever domain scores differ in dispersion, covary with one another,
  or are missing for some rows; this function makes that visible instead
  of leaving it to be inferred from a figure.

### Fixes

- [`plot_rri_properties()`](https://mghotbi.github.io/HRRI/reference/plot_rri_properties.md)
  no longer draws a literal “NA” label at the tip of an unavailable
  axis. The label was placed at the same radius as the axis name and
  overlapped it, rendering “Capacity” as “CNAacity”. Unavailable axes
  now show the short spoke only.
- [`plot_rri_properties()`](https://mghotbi.github.io/HRRI/reference/plot_rri_properties.md)
  dropped two spurious warnings: an invalid `label.size` argument passed
  to
  [`annotate()`](https://ggplot2.tidyverse.org/reference/annotate.html),
  and a removed-rows warning from plotting an `NA` vertex.
- [`rri_property_scores()`](https://mghotbi.github.io/HRRI/reference/rri_property_scores.md)
  now reports which input is missing when a property is unavailable.
  Capacity requires `soil_df`; omitting it is the usual reason Capacity
  returns `NA`. The worked example shows both the complete and the
  partial call.
- `R/plot_rri_recovery_landscape.R` used a literal em dash inside a
  string, which is non-ASCII in executable code. Replaced with `\u2014`;
  the rendered glyph is unchanged.

### Documentation

- [`rri_recovery_metrics()`](https://mghotbi.github.io/HRRI/reference/rri_recovery_metrics.md)
  gains a section explaining why `k_recovery` and `t_half` are
  frequently `NA`. The rate fit requires at least `min_pts` recovery
  observations that fall after the trough, below baseline, and before
  the first baseline crossing; short windows routinely leave fewer. The
  worked example now uses a 40-day record with a 10-day event so the fit
  succeeds, and contrasts it with a short window where the rate is
  correctly withheld.

## HRRI 0.99.2

### Scientific change

- Memory (`M`) is now a holobiont state rather than a mineralogical one.
  It accumulates from hydrological event load, the Fe-crystallinity
  ratchet, persistent plant acclimation (aerenchyma displacement), and
  microbial community displacement. Component weights sum to the
  previous total (0.060), so overall accumulation rate is unchanged
  while its causes are distributed across the three domains.
- Added `micro_legacy` as a simulator state: it accrues under sustained
  reduction and relaxes more slowly under oxic recovery, giving the
  microbial component an asymmetric (hysteretic) legacy.
- `latent_state` now returns `micro_legacy` and `plant_legacy` so the
  memory decomposition is auditable.
- Memory decay remains keyed to Fe crystallinity, which is the least
  reversible component and sets the floor on memory loss.
- **This changes `latent_truth`.** All figures, benchmarks and any
  stored simulation outputs must be regenerated.

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
