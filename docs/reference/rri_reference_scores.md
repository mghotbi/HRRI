# Score departures from an explicitly defined reference

An optional, transparent alternative to latent-axis scoring. Each
feature score is max(0, 1 - abs(value - target) / tolerance). The result
measures proximity to the declared reference, not validated ecosystem
functioning or a universal resilience scale.

## Usage

``` r
rri_reference_scores(
  data,
  reference,
  id = NULL,
  domain_weights = c(Physio = 0.4, Soil = 0.35, Micro = 0.25),
  min_coverage = 0.5,
  na_policy = c("available", "complete")
)
```

## Arguments

- data:

  Numeric feature data frame; rows must be aligned with id.

- reference:

  Data frame with feature, domain, target, tolerance, weight. Domains
  are Physio, Soil or Micro. Tolerance is a positive distance from
  target at which the feature score reaches zero. Reference rows for
  unmeasured features may be retained to report coverage against a
  common panel.

- id:

  Optional aligned identifiers.

- domain_weights:

  Named, non-negative domain weights.

- min_coverage:

  Minimum weighted within-domain feature coverage.

- na_policy:

  Use available domains or require every positive-weight domain.

## Value

An RRI object with fixed-reference domain scores, feature scores,
coverage, effective row-specific domain weights, and reference metadata.

## Examples

``` r
dat <- data.frame(Eh = c(100, 200, 150), pH = c(5.5, 6.0, 5.8))
ref <- data.frame(feature = c("Eh", "pH"), domain = c("Soil", "Soil"),
                  target = c(150, 5.8), tolerance = c(100, 0.5), weight = 1)
rri_reference_scores(dat, ref)$row_scores
#>   Physio Soil Micro  RRI domain_coverage n_domains
#> 1     NA 0.45    NA 0.45            0.35         1
#> 2     NA 0.55    NA 0.55            0.35         1
#> 3     NA 1.00    NA 1.00            0.35         1
```
