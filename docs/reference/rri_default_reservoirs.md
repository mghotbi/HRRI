# Illustrative reservoir parameter template

Returns example values, not calibrated mineral-specific constants. The
default uses the core bulk EAC and EDC columns. Supply phase-resolved
column names only when those inventories are nonoverlapping and
expressed in electron-equivalent units.

## Usage

``` r
rri_default_reservoirs(
  eac_ferrihydrite_col = NULL,
  eac_goethite_col = NULL,
  eac_structural_col = NULL,
  edc_humic_fast_col = NULL,
  edc_humic_slow_col = NULL,
  edc_eac_col = "EAC",
  edc_edc_col = "EDC"
)
```

## Arguments

- eac_ferrihydrite_col, eac_goethite_col, eac_structural_col:

  EAC phase columns.

- edc_humic_fast_col, edc_humic_slow_col:

  EDC fraction columns.

- edc_eac_col, edc_edc_col:

  Bulk EAC and EDC fallback columns. Argument names are retained for
  backward compatibility.

## Value

Named list for rri_accessible_capacity; k uses inverse hours.

## Examples

``` r
rri_default_reservoirs()
#> $bulk_EAC
#> $bulk_EAC$Q_col
#> [1] "EAC"
#> 
#> $bulk_EAC$alpha
#> [1] 0.5
#> 
#> $bulk_EAC$k
#> [1] 0.2
#> 
#> $bulk_EAC$type
#> [1] "EAC"
#> 
#> 
#> $bulk_EDC
#> $bulk_EDC$Q_col
#> [1] "EDC"
#> 
#> $bulk_EDC$alpha
#> [1] 0.45
#> 
#> $bulk_EDC$k
#> [1] 0.15
#> 
#> $bulk_EDC$type
#> [1] "EDC"
#> 
#> 
```
