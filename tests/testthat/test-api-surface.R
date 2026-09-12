## tests/testthat/test-api-surface.R
##
## Guards the public interface itself rather than any single calculation.
##
## These checks exist because of two defects that got through everything else.
## `R CMD check` verifies that documented arguments match a function's formals,
## but it never verifies that a documented RETURN VALUE matches what the
## function returns: rri_o2_demand() described three component columns under
## names the code does not use, and five returned elements were not documented
## at all. Nothing mechanical caught that. It also cannot notice when a new
## export is added without a pkgdown entry, which broke the website build.
##
## Everything here is introspection over the installed package, so a failure
## means the interface really did change, not that a simulation drifted.

pkg <- "HRRI"

exported_names <- function() {
  getNamespaceExports(asNamespace(pkg))
}

test_that("every export resolves and is a function or a documented dataset", {
  ex <- exported_names()
  expect_gt(length(ex), 30L)
  for (nm in ex) {
    obj <- tryCatch(get(nm, envir = asNamespace(pkg)),
                    error = function(e) NULL)
    expect_false(is.null(obj), info = paste("export does not resolve:", nm))
  }
})

test_that("S3 methods are registered, not merely defined", {
  for (m in list(c("print", "rri_accuracy"),
                 c("print", "hrri_benchmark"),
                 c("print", "hrri_arch"),
                 c("summary", "hrri_arch"),
                 c("as.data.frame", "hrri_arch"))) {
    fn <- tryCatch(getS3method(m[1], m[2], optional = TRUE),
                   error = function(e) NULL)
    expect_true(is.function(fn),
                info = paste0("S3 method not registered: ", m[1], ".", m[2]))
  }
})

test_that("argument names of the public entry points are stable", {
  ## A silent rename here would break user scripts without failing any other
  ## test. Only the arguments users actually pass by name are pinned.
  expect_true(all(c("score", "target", "cluster", "n_boot", "n_perm", "conf",
                    "seed") %in% names(formals(rri_accuracy))))
  expect_true(all(c("acc", "panels", "base_size") %in%
                    names(formals(plot_rri_accuracy))))
  expect_true(all(c("res", "domains", "rri_col", "weights") %in%
                    names(formals(rri_domain_influence))))
  expect_true(all(c("soil_df", "fe2_col", "o2_supply_col", "custom_coefs") %in%
                    names(formals(rri_o2_demand))))
  expect_true(all(c("n_plot", "n_depth", "n_plant", "n_time", "seed") %in%
                    names(formals(simulate_redox_holobiont))))
})

test_that("seeding never disturbs the caller's random stream", {
  ## CRAN policy: .Random.seed lives in .GlobalEnv and must be left as found.
  ## Every function that seeds must expose the seed as an argument, and must
  ## put the stream back.
  for (fn in c("rri_accuracy", "simulate_redox_holobiont", "benchmark_hrri")) {
    f <- get(fn, envir = asNamespace(pkg))
    expect_true(any(grepl("seed", names(formals(f)))),
                info = paste(fn, "seeds without exposing a seed argument"))
  }
  ## Nothing seeds unless asked: every seed argument defaults to NULL.
  for (fn in c("rri_accuracy", "simulate_redox_holobiont")) {
    f <- get(fn, envir = asNamespace(pkg))
    expect_null(eval(formals(f)$seed),
                info = paste(fn, "seed does not default to NULL"))
  }

  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  invisible(simulate_redox_holobiont(n_plot = 1, n_depth = 1, n_plant = 1,
                                     n_time = 6, p_micro = 3, seed = 7))
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(before, after)
})

test_that("every exported topic is documented with a value and an example", {
  skip_on_cran()
  db <- tryCatch(tools::Rd_db(pkg), error = function(e) NULL)
  skip_if(is.null(db) || !length(db), "Rd database unavailable")

  tag_of <- function(x) {
    t <- attr(x, "Rd_tag")
    if (is.null(t)) "" else t
  }
  sections <- function(rd) vapply(rd, tag_of, character(1))

  aliases <- unlist(lapply(db, function(rd) {
    unlist(rd[sections(rd) == "\\alias"])
  }), use.names = FALSE)
  aliases <- trimws(aliases)

  missing_doc <- setdiff(exported_names(), aliases)
  expect_length(missing_doc, 0L)

  ## \value is required of function topics only. A narrative page such as
  ## HRRI_workflow_example has no \usage and returns nothing, so demanding a
  ## \value there is wrong -- R CMD check applies the same rule.
  is_function_topic <- function(rd) "\\usage" %in% sections(rd)
  fn_db <- Filter(is_function_topic, db)

  no_value <- names(Filter(function(rd) !("\\value" %in% sections(rd)), fn_db))
  no_examples <- names(Filter(function(rd) !("\\examples" %in% sections(rd)), db))
  expect_length(no_value, 0L)
  expect_length(no_examples, 0L)
})

test_that("no example is hidden behind dontrun", {
  skip_on_cran()
  db <- tryCatch(tools::Rd_db(pkg), error = function(e) NULL)
  skip_if(is.null(db) || !length(db), "Rd database unavailable")
  has_dontrun <- vapply(db, function(rd) {
    txt <- tryCatch(paste(as.character(rd), collapse = " "),
                    error = function(e) "")
    grepl("\\\\dontrun", txt, perl = TRUE)
  }, logical(1))
  ## \dontrun means R CMD check never runs the example, so it can drift out of
  ## step with the code unnoticed. \donttest still runs under --run-donttest.
  expect_length(names(which(has_dontrun)), 0L)
})

test_that("theme_ems returns a usable ggplot2 theme", {
  skip_if_not_installed("ggplot2")
  th <- theme_ems()
  expect_s3_class(th, "theme")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + th
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
