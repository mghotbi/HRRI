# Stoichiometric O\\\_2\\ Demand from Reduced-Pool Inventories

Computes the complete-oxidation O\\\_2\\ demand for a given rhizosphere
reduced-pool inventory and compares it to the specified O\\\_2\\ stock
on the same dry-soil mass basis:

\$\$O_2^{\mathrm{demand}} = \sum\_{j} n_j \cdot s_j\$\$

where \\n_j\\ is the molar inventory (mmol kg\\^{-1}\\) of reduced
species \\j\\ and \\s_j\\ is the stoichiometric O\\\_2\\ coefficient for
complete oxidation to the specified endpoint.

Stoichiometric coefficients follow the electron balance table of Ghotbi,
Ghotbi, Stukenbrock, Mühling and Spielvogel (2026) (Box 1 of the
mechanistic review):

|  |  |  |
|----|----|----|
| **Reduced pool** | **Endpoint** | **O\\\_2\\ (mol mol\\^{-1}\\)** |
| Fe\\^{2+}\\ | Fe(III) oxyhydroxide | 0.25 |
| Mn\\^{2+}\\ | MnO\\\_2\\ | 0.50 |
| HS\\^{-}\\ | SO\\\_4^{2-}\\ | 2.00 |
| FeS (mackinawite) | Fe(III) + SO\\\_4^{2-}\\ | 2.25 |
| FeS\\\_2\\ (pyrite) | Fe(III) + 2SO\\\_4^{2-}\\ | 3.75 |
| NH\\\_4^+\\ | NO\\\_3^-\\ | 2.00 |
| CH\\\_4\\ | CO\\\_2\\ | 2.00 |
| Acetate equivalents | CO\\\_2\\ | 2.00 |

The *O\\\_2\\ deficit ratio* (`O2_deficit_ratio`) is
\\O_2^{\mathrm{demand}} / O_2^{\mathrm{supply}}\\: values \\\>1\\
indicate that the specified O\\\_2\\ stock is smaller than the demand.
This stock ratio does not determine recovery or account for continuing
O2 delivery.

## Usage

``` r
rri_o2_demand(
  soil_df,
  fe2_col = NULL,
  mn2_col = NULL,
  hs_col = NULL,
  fes_col = NULL,
  fes2_col = NULL,
  nh4_col = NULL,
  ch4_col = NULL,
  acetate_col = NULL,
  ch4_unit = c("mmol_kg", "umol_kg"),
  acetate_basis = c("acetate", "carbon"),
  o2_supply_col = NULL,
  custom_coefs = NULL,
  bulk_density = NULL,
  theta_v = NULL,
  particle_density = 2.65,
  return_components = TRUE
)
```

## Arguments

- soil_df:

  Data frame with soil chemistry (rows = samples).

- fe2_col:

  Character or `NULL`. Column for Fe\\^{2+}\\ (mmol kg\\^{-1}\\).
  Default stoichiometric coefficient: 0.25.

- mn2_col:

  Character or `NULL`. Column for Mn\\^{2+}\\ (mmol kg\\^{-1}\\).
  Coefficient: 0.50.

- hs_col:

  Character or `NULL`. Column for HS\\^{-}\\ (mmol kg\\^{-1}\\).
  Coefficient: 2.00.

- fes_col:

  Character or `NULL`. Column for FeS/mackinawite (mmol kg\\^{-1}\\).
  Coefficient: 2.25.

- fes2_col:

  Character or `NULL`. Column for FeS\\\_2\\/pyrite (mmol kg\\^{-1}\\).
  Coefficient: 3.75.

- nh4_col:

  Character or `NULL`. Column for NH\\\_4^+\\ (mmol kg\\^{-1}\\).
  Coefficient: 2.00.

- ch4_col:

  Character or `NULL`. Column for CH\\\_4\\.

- acetate_col:

  Character or `NULL`. Column for dissolved organic matter expressed as
  acetate or acetate-carbon equivalents.

- ch4_unit:

  Character. Unit of `ch4_col`: `"mmol_kg"` (default) or `"umol_kg"`.

- acetate_basis:

  Character. `"acetate"` (default; 2 mol O2 per mol acetate) or
  `"carbon"` (1 mol O2 per mol acetate-C).

- o2_supply_col:

  Character or `NULL`. Column for an explicitly defined O\\\_2\\
  inventory (mmol O\\\_2\\ kg\\^{-1}\\). If `NULL`, deficit ratio is not
  computed.

- custom_coefs:

  Optional named numeric vector to override or extend the built-in
  stoichiometric coefficients. Names must match the argument names above
  (e.g., `c(fe2_col = 0.25)`). Useful for system-specific endpoint
  assumptions.

- bulk_density:

  Numeric of length one or nrow(soil_df). Soil bulk density (g
  cm\\^{-3}\\) for volumetric conversion of demand to mmol O\\\_2\\
  L\\^{-1}\\ porewater. Set to `NULL` to skip volumetric conversion.

- theta_v:

  Numeric of length one or nrow(soil_df), or `NULL`. Volumetric water
  content (L water L\\^{-1}\\ bulk soil). When omitted and
  `bulk_density` is supplied, saturated porosity is estimated as
  `1 - bulk_density / particle_density`.

- particle_density:

  Numeric. Particle density in g cm\\^{-3}\\ used only for the
  saturated-porosity estimate. Default 2.65.

- return_components:

  Logical. If `TRUE` (default), return a per-species contribution
  matrix.

## Value

A list:

- `o2_demand`:

  Numeric vector (mmol O\\\_2\\ kg\\^{-1} dry soil\\) of
  complete-oxidation O\\\_2\\ demand per sample.

- `o2_deficit_ratio`:

  Per-sample \\O_2^{\mathrm{demand}} / O_2^{\mathrm{supply}}\\. Values
  \\\>1\\ indicate demand exceeds the specified stock. `NA` when
  `o2_supply_col` is absent.

- `o2_demand_vol`:

  Volumetric O\\\_2\\ demand (mmol L\\^{-1}\\ porewater); `NA` if
  `bulk_density` is `NULL`.

- `n_species_used`:

  Integer. Number of reduced-pool columns found.

- `components`:

  Data frame (one row per species) with: `species`, `stoich_coef`,
  `mean_inventory`, `mean_o2_contribution`, `fraction_of_total_demand`.
  Returned only when `return_components = TRUE`.

- `stoich_table`:

  Data frame of the full stoichiometric table used, including any custom
  overrides.

## Details

**Interpretation — the 25-fold contrast.**

The mechanistic review (Ghotbi *et al.*, 2026) provides a worked
example: a Fe-rich rhizosphere containing 50 mmol Fe(II) kg\\^{-1}\\, 5
mmol FeS kg\\^{-1}\\, 2 mmol Mn(II) kg\\^{-1}\\, 2 mmol NH\\\_4^+\\
kg\\^{-1}\\, 2 mmol acetate kg\\^{-1}\\, and 0.5 mmol CH\\\_4\\
kg\\^{-1}\\ has a complete-oxidation O\\\_2\\ ceiling of \\\approx 34\\
mmol O\\\_2\\ kg\\^{-1}\\, versus only \\\approx 1.3\\ mmol O\\\_2\\
kg\\^{-1}\\ in an illustrative initial pore-gas stock: about a 26-fold
contrast. The assumed stock is not air-saturated porewater. Continuing
atmospheric and root O2 delivery can replenish it. Accessibility,
reaction kinetics and transport determine realized demand during an
event.

**Pyrite stoichiometry.**

Complete pyrite oxidation to sulfate and Fe(III) oxyhydroxide releases 4
mol H\\^+\\ mol\\^{-1}\\ FeS\\\_2\\ and consumes 3.75 mol O\\\_2\\:
FeS\\\_2\\ + 15/4 O\\\_2\\ + 7/2 H\\\_2\\O → Fe(OH)\\\_3\\ +
2H\\\_2\\SO\\\_4\\. Partial oxidation to sulfur intermediates (S\\^0\\,
thiosulfate) requires fewer moles; adjust via `custom_coefs`.

**pH coupling.**

Under the applicable rate law, homogeneous abiotic Fe(II) oxidation
increases \\\approx\\100-fold per unit pH rise (Stumm & Lee, 1961;
Millero *et al.*, 1987). The stoichiometric demand computed here is for
complete oxidation and is independent of pH, but actual O\\\_2\\
consumption rates will be pH-modulated. This function reports a
stoichiometric potential demand, not a thermodynamic limit or rate.

## References

Ghotbi, M., Ghotbi, M., Stukenbrock, E. H., Mühling, K. H., &
Spielvogel, S. (2026). Rhizosphere redox recovery after hydrological
disturbance: mechanisms across the soil–plant–microbiome continuum.
*Manuscript submitted*.

Stumm, W., & Lee, G. F. (1961). Oxygenation of ferrous iron. *Industrial
& Engineering Chemistry*, 53, 143–146.

Millero, F. J., Sotolongo, S., & Izaguirre, M. (1987). The oxidation
kinetics of Fe(II) in seawater. *Geochimica et Cosmochimica Acta*, 51,
793–801.

## See also

[`rri_accessible_capacity`](https://mghotbi.github.io/HRRI/reference/rri_accessible_capacity.md),
[`rri_capacity_index`](https://mghotbi.github.io/HRRI/reference/rri_capacity_index.md),
[`rri_root_physio`](https://mghotbi.github.io/HRRI/reference/rri_root_physio.md)

## Examples

``` r
## Reproduce the worked example from Ghotbi et al. (2026) Box 1
worked_example <- data.frame(
  Fe2 = 50.0, # mmol kg-1
  FeS = 5.0,
  Mn2 = 2.0,
  NH4 = 2.0,
  acetate = 2.0,
  CH4 = 0.5, # mmol kg-1; set ch4_unit = "umol_kg" for umol input
  O2_pw = 1.3 # assumed initial O2 stock, mmol/kg; NOT air-saturated porewater
)

demand <- rri_o2_demand(
  soil_df = worked_example,
  fe2_col = "Fe2",
  mn2_col = "Mn2",
  fes_col = "FeS",
  nh4_col = "NH4",
  acetate_col = "acetate",
  ch4_col = "CH4",
  o2_supply_col = "O2_pw",
  return_components = TRUE
)

demand$o2_demand # should be ~34 mmol O2 kg-1
#> [1] 33.75
demand$o2_deficit_ratio # should be ~26
#> [1] 25.96154
demand$components
#>   species n_observed stoich_coef mean_inventory_mmol mean_o2_contribution
#> 1     Fe2          1        0.25                50.0                12.50
#> 2     Mn2          1        0.50                 2.0                 1.00
#> 3     FeS          1        2.25                 5.0                11.25
#> 4     NH4          1        2.00                 2.0                 4.00
#> 5     CH4          1        2.00                 0.5                 1.00
#> 6 acetate          1        2.00                 2.0                 4.00
#>   fraction_total_demand
#> 1                0.3704
#> 2                0.0296
#> 3                0.3333
#> 4                0.1185
#> 5                0.0296
#> 6                0.1185
```
