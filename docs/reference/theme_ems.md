# EMS plotting theme

A simple ggplot theme used for HRRI visualizations.

## Usage

``` r
theme_ems(base_size = 12)
```

## Arguments

- base_size:

  Base font size

## Value

A ggplot2 theme object

## Examples

``` r
# \donttest{
  library(ggplot2)
#> Warning: package ‘ggplot2’ was built under R version 4.5.2
  ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
    geom_point() + theme_ems()

# }
```
