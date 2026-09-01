# Event-window accessible reservoir capacities

Computes Q \* alpha \* (1-exp(-k\*tau)) for declared reservoirs. Q must
be nonoverlapping electron-equivalent inventories with explicit reaction
endpoints. EAC and EDC are returned separately. Their sum is an
inventory descriptor, not electron flux; their difference is not an
oxygen budget.

## Usage

``` r
rri_accessible_capacity(
  soil_df,
  reservoirs,
  tau = 24,
  normalise = TRUE,
  return_components = FALSE
)
```

## Arguments

- soil_df:

  Numeric reservoir measurements.

- reservoirs:

  Nonempty named list of Q_col, alpha, k and type (EAC/EDC). alpha in
  `[0, 1]` and k \>= 0 may be scalars, row vectors or column names.
  Parameters must be specified for the relevant process and conditions;
  they are not identifiable separately from one accessible-capacity
  observation.

- tau:

  Non-negative duration, scalar or row vector; units reciprocal to k.

- normalise:

  Divide by the sum of observed inventories. This produces an accessible
  fraction, not absolute capacity or guaranteed comparability.

- return_components:

  Include reservoir contribution summaries.

## Value

Capacities, observed subtotal, fraction and reservoir coverage. Missing
reservoir types remain NA. Partial rows are labelled observed subtotals;
absence is not zero. ck_limited is retained as NA because 0.30 is not a
validated threshold. Negative inventories are treated as missing.

## Details

Default reservoir parameters are illustrative scenario values only.
Fe(II) oxidation rates must not be assigned as Fe(III) reduction
constants. For capacity estimation use experimentally constrained
process-specific rates.

## Examples

``` r
df <- data.frame(EAC = c(10, 20, 30), EDC = c(5, 8, 12))
res <- rri_accessible_capacity(
  df,
  reservoirs = list(
    bulk_EAC = list(Q_col = "EAC", alpha = 0.5, k = 0.2, type = "EAC"),
    bulk_EDC = list(Q_col = "EDC", alpha = 0.45, k = 0.15, type = "EDC")
  ),
  tau = 24
)
res$cacc
#> [1] 0.4764915 0.4792620 0.4792620
```
