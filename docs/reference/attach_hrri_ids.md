# Attach design identifiers to a score table with explicit alignment checks

Joins experimental design identifiers onto a pipeline score table.
Alignment is established by a shared unique `row_id`, or by a complete
shared observation key, or - only as a last resort and with an explicit
warning - by preserved input row order. Matching row counts alone do not
establish alignment, so the method actually used is recorded in the
`id_alignment` attribute of the returned data frame.

## Usage

``` r
attach_hrri_ids(scores, id, key = NULL)
```

## Arguments

- scores:

  Data frame of row-level scores, typically `res$row_scores` from
  [`rri_pipeline`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md)
  or
  [`rri_pipeline_st`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md).

- id:

  Data frame of design identifiers, typically `sim$id` from
  [`simulate_redox_holobiont`](https://mghotbi.github.io/HRRI/reference/simulate_redox_holobiont.md).

- key:

  Optional character vector naming the observation key columns to join
  on. `NULL` (default) selects a key automatically: `row_id` when
  present and unique in both inputs, otherwise the intersection of
  `c("plot", "depth", "plant_id", "time")` present in both inputs.

## Value

`scores` with the non-conflicting columns of `id` attached. The
`id_alignment` attribute is a list giving the `method` (`"row_id"`,
`"observation_key"` or `"row_order"`), the `key` columns used, and the
number of rows matched. Shared identifier columns already present in
`scores` are checked for conflicts rather than silently overwritten.

## See also

[`rri_pipeline`](https://mghotbi.github.io/HRRI/reference/rri_pipeline.md),
[`rri_pipeline_st`](https://mghotbi.github.io/HRRI/reference/rri_pipeline_st.md)

## Examples

``` r
scores <- data.frame(row_id = 1:4, RRI = c(0.4, 0.6, 0.5, 0.7))
ids <- data.frame(row_id = 1:4,
                  plot = c("P1", "P1", "P2", "P2"),
                  time = c(1, 2, 1, 2))
out <- attach_hrri_ids(scores, ids)
attr(out, "id_alignment")$method
#> [1] "row_id"
head(out)
#>   plot time row_id RRI
#> 1   P1    1      1 0.4
#> 2   P1    2      2 0.6
#> 3   P2    1      3 0.5
#> 4   P2    2      4 0.7
```
