# Ternary Plot of Relative Domain Scores

Creates a ternary diagram of the relative magnitudes of plant, soil and
microbial domain scores after closure to a unit sum. These coordinates
are display quantities, not fractions of causal buffering capacity.
Points are filled according to the corresponding composite `RRI` value.

## Usage

``` r
plot_RRI_ternary(
  ternary_df,
  point_size = 5,
  point_alpha = 0.9,
  palette = "plasma",
  show_subtitle = TRUE,
  show_centroid = TRUE,
  centroid_shape = 23,
  centroid_size = 1.4,
  tolerance = 1e-06,
  renormalize = FALSE,
  centroid_method = c("auto", "simplex_mean", "aitchison_mean")
)
```

## Arguments

- ternary_df:

  A data frame containing compositional columns `Physio`, `Soil`,
  `Micro`, and `RRI`.

- point_size:

  Numeric; size of ternary points.

- point_alpha:

  Numeric between 0 and 1 controlling point transparency.

- palette:

  Character; viridis palette option.

- show_subtitle:

  Logical; display system-level RRI mean in subtitle.

- show_centroid:

  Logical; add compositional centroid marker.

- centroid_shape:

  Numeric; ggplot2 shape for centroid marker.

- centroid_size:

  Numeric multiplier for centroid size.

- tolerance:

  Numeric; tolerance used for compositional closure checks.

- renormalize:

  Logical; if TRUE, renormalises rows that do not sum to 1.

- centroid_method:

  Character; one of `"auto"`, `"simplex_mean"`, or `"aitchison_mean"`.

## Value

A `ggtern` object.

## Details

Closure removes absolute score magnitude: rows with proportional domain
scores occupy the same position even when their composite scores differ.
Do not infer mechanistic allocation, causal contribution or substitution
from this plot. If clr-transformed coordinates are attached as an
attribute (`"clr"`), the centroid can be computed using the Aitchison
mean. Otherwise, a simplex arithmetic mean is used.

## Examples

``` r
if (FALSE) { # \dontrun{
# Optional packages ggtern and viridis are required.
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 2,
  n_plant = 2,
  n_time = 8,
  p_micro = 6,
  seed = 1234
)

# ---- Compute RedoxRRI ----
res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)

# ---- Extract compositional scores ----
ternary_df <- res$row_scores_comp

# ---- Plot ternary allocation ----
p <- plot_RRI_ternary(
  ternary_df,
  point_size = 3,
  show_centroid = TRUE
)
} # }
```
