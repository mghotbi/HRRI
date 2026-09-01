# Construct an explicitly weighted microbial guild contrast

Summarises supplied guild measurements or proxies. The weights define a
contrast, not a universal ordering of microbial resilience.
Denitrification, sulfate reduction and methanogenesis may be beneficial
or detrimental depending on the specified ecosystem function and
disturbance.

## Usage

``` r
rri_micro_functional_score(micro_traits, weights = NULL, scale = TRUE)
```

## Arguments

- micro_traits:

  Numeric data frame with comparable, justified guild scales. Gene
  abundance or expression does not by itself measure process rate.

- weights:

  Named finite signed weights. Default legacy weights are illustrative
  only and trigger a warning; supply scientifically justified weights.

- scale:

  If TRUE, scale finite contrasts within the supplied cohort to
  `[0, 1]`. Constants map to 0.5 and wholly unobserved rows remain NA.

## Value

Numeric vector with coverage and raw_contrast attributes. The raw
contrast is a signed weighted sum divided by the available absolute
weight. Missing guilds are not zeros; changing availability changes the
estimand.

## See also

rri_reference_scores

## Examples

``` r
x <- data.frame(EET_reduction = c(0.2, 0.4, NA),
                methanogenesis = c(0.3, 0.1, NA))
rri_micro_functional_score(x, weights = c(EET_reduction = 1, methanogenesis = -1))
#> [1]  0  1 NA
#> attr(,"coverage")
#> [1] 1 1 0
#> attr(,"raw_contrast")
#> [1] -0.05  0.15    NA
```
