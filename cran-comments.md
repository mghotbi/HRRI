## Test environments

* local: macOS 26.6 (aarch64-apple-darwin24.4.0), R 4.5.2
* win-builder: R-devel, R-release            <!-- DELETE IF NOT RUN -->
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest,
  windows-latest                             <!-- DELETE IF NOT RUN -->

## R CMD check results

0 errors | 0 warnings | 1 note

> New submission

This is a first submission.

Locally an additional note reports that math rendering was not checked because
the 'V8' package is unavailable on the test machine. This concerns the check
environment rather than the package.

## Notes on content

* One example is conditional. `plot_RRI_ternary()` needs the suggested package
  **ggtern**, which patches ggplot2's internal element tree when its namespace
  is loaded. Under ggplot2 >= 4.0.0 that patch makes every subsequent plot in
  the session fail, so the example, the vignette and the function all test the
  ggplot2 version before ggtern is loaded and skip the ternary panel rather
  than error. No example uses `\dontrun{}`; all run under `--run-donttest`.

* Random number generation. Functions that seed expose the seed as an argument
  and restore both the RNG kind and `.Random.seed` via `on.exit()`, so a
  caller's stream is returned exactly as found. This is asserted in
  `tests/testthat/test-api-surface.R`.

  `simulate_redox_holobiont()` is an illustrative data generator and keeps a
  fixed default of `seed = 123` rather than `NULL`, so that an unseeded call
  reproduces the output shown in the vignettes. Since the caller's stream is
  restored on exit, the user's session state is unchanged. Functions that
  consume randomness for inference, such as `rri_accuracy()`, default to
  `seed = NULL` and do not seed unless asked.

* The package writes nothing to the user's home filespace or working
  directory, creates no files in the temporary directory, and does not
  otherwise modify the global environment.

* References. The `Description` field cites only published, resolvable
  sources. Three companion manuscripts describing the framework are not yet
  published; they are recorded in `inst/CITATION`, README and the vignettes
  instead, and will be added to `DESCRIPTION` with DOIs on publication.
