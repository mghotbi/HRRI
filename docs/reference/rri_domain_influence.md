# Realised influence of each domain on the composite score

A declared weight is not the same thing as realised influence. Two
domains given equal weight contribute unequally to the composite
whenever their scores differ in dispersion, in how strongly they covary
with the other domains, or in how often they are missing. This function
reports what each domain actually contributed, so that a claim such as
"the index is driven by soil" can be checked rather than inferred from
the pattern of a figure.

## Usage

``` r
rri_domain_influence(
  res,
  domains = c("Physio", "Soil", "Micro"),
  rri_col = "RRI",
  weights = NULL
)
```

## Arguments

- res:

  An `RRI` object from
  [`rri_pipeline()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  or
  [`rri_pipeline_st()`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md).

- domains:

  Character vector of domain score columns. Defaults to
  `c("Physio", "Soil", "Micro")`.

- rri_col:

  Name of the composite column. Default `"RRI"`.

- weights:

  Optional named numeric vector of the nominal weights used to build the
  composite. If `NULL` (default) the function tries
  `res$effective_weights`, then `res$meta$weights`, and otherwise
  reports realised influence without a nominal comparison.

## Value

A list with three elements.

- `influence`:

  One row per domain: `mean`, `sd`, `n_missing`, `cor_with_rri`,
  `nominal_weight`, `realised_share` and `ratio`.

- `covariance`:

  Pairwise correlations between domain scores. Strong cross-domain
  correlation means influence cannot be attributed cleanly to one
  domain.

- `notes`:

  Character vector of diagnostics worth acting on.

## Details

**How realised share is computed.** For weights \\w_d\\ and domain
scores \\S_d\\, the composite is \\R = \sum_d w_d S_d\\. Because
\\\mathrm{Var}(R) = \sum_d w_d \mathrm{Cov}(S_d, R)\\, the quantity

\$\$\phi_d = w_d \\ \mathrm{Cov}(S_d, R) / \mathrm{Var}(R)\$\$

is an exact decomposition: the \\\phi_d\\ sum to one. Each domain's
share therefore includes its own variance *and* its share of the
covariance it has with the other domains. This is the appropriate
attribution when domains are correlated, which they generally are under
a shared forcing.

**How to read `ratio`.** `ratio` is realised share divided by nominal
weight. A value near 1 means the domain influenced the composite about
as much as intended. Values above roughly 1.3 or below roughly 0.7
indicate that the declared weights are not delivering the intended
balance, usually for one of three reasons: the domain score is more (or
less) dispersed than the others after scaling; it covaries strongly with
the others, so it absorbs shared variance; or it is missing for many
rows, so per-row weight renormalisation quietly redistributes its
weight.

**What this does not establish.** A high realised share is a statement
about the score, not about the ecosystem. It does not show that the
domain is mechanistically more important, and it is not evidence that
the composite is wrong. It shows only where the variance in this
particular composite came from, for this cohort, under these weights and
this scaling.

## See also

[`rri_sensitivity()`](https://mghotbi.github.io/HRRI/reference/rri_sensitivity.md)
for the effect of alternative weight grids;
[`rri_compensation_index()`](https://mghotbi.github.io/HRRI/reference/rri_compensation_index.md)
for cross-domain asynchrony.

## Examples

``` r
sim <- simulate_redox_holobiont(
  n_plot = 2, n_depth = 2, n_plant = 3, n_time = 40,
  p_micro = 25, seed = 2026
)

res <- suppressWarnings(rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
))

infl <- rri_domain_influence(res)
infl$influence
#>   domain      mean        sd n_missing cor_with_rri nominal_weight
#> 1 Physio 0.4617876 0.2553977         0    0.6537610             NA
#> 2   Soil 0.4924688 0.3274460         0    0.2994660             NA
#> 3  Micro 0.4552182 0.2205703         0    0.5443328             NA
#>   realised_share ratio
#> 1             NA    NA
#> 2             NA    NA
#> 3             NA    NA
infl$notes
#> [1] "Nominal weights could not be recovered from `res`; supply `weights` to compare realised influence against intent."
```
