# EMS plotting theme

A simple ggplot theme used for RedoxRRI visualizations.

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
if (FALSE) { # \dontrun{
  library(ggplot2)
  ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) +
    geom_point() + theme_ems()
} # }
```
