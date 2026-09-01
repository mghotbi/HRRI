# Exploratory root-trait composite

A weighted standardized trait summary. Trait direction is
context-dependent: greater ROL, porosity, aerenchyma or SRL is not
universally better plant performance or greater oxygen delivery to every
root region.

## Usage

``` r
rri_root_physio(
  plant_df,
  biomass_col = NULL,
  length_col = NULL,
  rol_col = NULL,
  aerenchyma_col = NULL,
  porosity_col = NULL,
  srl_col = NULL,
  w_biomass = 0.2,
  w_length = 0.2,
  w_rol = 0.3,
  w_aerenchyma = 0.2,
  w_porosity = 0.05,
  w_srl = 0.05,
  directions = NULL,
  scaling = c("pnorm", "minmax")
)
```

## Arguments

- plant_df:

  Numeric root-trait measurements.

- biomass_col, length_col:

  Columns for root biomass and length/density.

- rol_col:

  Column for measured radial oxygen loss, with consistent units.

- aerenchyma_col, porosity_col:

  Optional aeration trait columns; avoid double weighting correlated
  measures of the same anatomical attribute.

- srl_col:

  Column for specific root length.

- w_biomass, w_length, w_rol, w_aerenchyma, w_porosity, w_srl:

  Non-negative weights.

- directions:

  Optional named numeric vector assigning +1 (larger maps to a larger
  score) or -1 (larger maps to a smaller score) to each selected
  measurement column. NULL uses +1 for backward compatibility and warns.

- scaling:

  pnorm (normal-CDF scaling, not a calibrated probability) or minmax.

## Value

root_physio_score and trait contributions; unobserved values stay NA.

## Examples

``` r
df <- data.frame(ROL = c(0.5, 1.0, 1.5), biomass = c(10, 20, 30))
rri_root_physio(df, rol_col = "ROL", biomass_col = "biomass",
  directions = c(ROL = 1, biomass = 1))$root_physio_score
#> [1] 0.1586553 0.5000000 0.8413447
```
