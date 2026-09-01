#' @title Stoichiometric O\eqn{_2} Demand from Reduced-Pool Inventories
#'
#' @description
#' Computes the complete-oxidation O\eqn{_2} demand for a given rhizosphere
#' reduced-pool inventory and compares it to the specified O\eqn{_2}
#' stock on the same dry-soil mass basis:
#'
#' \deqn{O_2^{\mathrm{demand}} = \sum_{j} n_j \cdot s_j}
#'
#' where \eqn{n_j} is the molar inventory (mmol kg\eqn{^{-1}}) of reduced
#' species \eqn{j} and \eqn{s_j} is the stoichiometric O\eqn{_2} coefficient
#' for complete oxidation to the specified endpoint.
#'
#' Stoichiometric coefficients follow the electron balance table of Ghotbi,
#' Ghotbi, Stukenbrock, Mühling and Spielvogel (2026) (Box 1 of the
#' mechanistic review):
#' \tabular{llr}{
#'   \strong{Reduced pool} \tab \strong{Endpoint} \tab \strong{O\eqn{_2} (mol mol\eqn{^{-1}})}\cr
#'   Fe\eqn{^{2+}} \tab Fe(III) oxyhydroxide \tab 0.25 \cr
#'   Mn\eqn{^{2+}} \tab MnO\eqn{_2} \tab 0.50 \cr
#'   HS\eqn{^{-}} \tab SO\eqn{_4^{2-}} \tab 2.00 \cr
#'   FeS (mackinawite) \tab Fe(III) + SO\eqn{_4^{2-}} \tab 2.25 \cr
#'   FeS\eqn{_2} (pyrite) \tab Fe(III) + 2SO\eqn{_4^{2-}} \tab 3.75 \cr
#'   NH\eqn{_4^+} \tab NO\eqn{_3^-} \tab 2.00 \cr
#'   CH\eqn{_4} \tab CO\eqn{_2} \tab 2.00 \cr
#'   Acetate equivalents \tab CO\eqn{_2} \tab 2.00 \cr
#' }
#'
#' The \emph{O\eqn{_2} deficit ratio} (\code{O2_deficit_ratio}) is
#' \eqn{O_2^{\mathrm{demand}} / O_2^{\mathrm{supply}}}: values \eqn{>1}
#' indicate that the specified O\eqn{_2} stock is smaller than the demand.
#' This stock ratio does not determine recovery or account for continuing O2 delivery.
#'
#' @param soil_df Data frame with soil chemistry (rows = samples).
#' @param fe2_col     Character or \code{NULL}. Column for Fe\eqn{^{2+}}
#'   (mmol kg\eqn{^{-1}}). Default stoichiometric coefficient: 0.25.
#' @param mn2_col     Character or \code{NULL}. Column for Mn\eqn{^{2+}}
#'   (mmol kg\eqn{^{-1}}). Coefficient: 0.50.
#' @param hs_col      Character or \code{NULL}. Column for HS\eqn{^{-}}
#'   (mmol kg\eqn{^{-1}}). Coefficient: 2.00.
#' @param fes_col     Character or \code{NULL}. Column for FeS/mackinawite
#'   (mmol kg\eqn{^{-1}}). Coefficient: 2.25.
#' @param fes2_col    Character or \code{NULL}. Column for FeS\eqn{_2}/pyrite
#'   (mmol kg\eqn{^{-1}}). Coefficient: 3.75.
#' @param nh4_col     Character or \code{NULL}. Column for NH\eqn{_4^+}
#'   (mmol kg\eqn{^{-1}}). Coefficient: 2.00.
#' @param ch4_col     Character or \code{NULL}. Column for CH\eqn{_4}.
#' @param ch4_unit Character. Unit of \code{ch4_col}: \code{"mmol_kg"}
#'   (default) or \code{"umol_kg"}.
#' @param acetate_col Character or \code{NULL}. Column for dissolved organic
#'   matter expressed as acetate or acetate-carbon equivalents.
#' @param acetate_basis Character. \code{"acetate"} (default; 2 mol O2 per
#'   mol acetate) or \code{"carbon"} (1 mol O2 per mol acetate-C).
#' @param o2_supply_col Character or \code{NULL}. Column for an explicitly defined O\eqn{_2}
#'   inventory (mmol O\eqn{_2} kg\eqn{^{-1}}). If \code{NULL}, deficit
#'   ratio is not computed.
#' @param custom_coefs Optional named numeric vector to override or extend the
#'   built-in stoichiometric coefficients. Names must match the argument names
#'   above (e.g., \code{c(fe2_col = 0.25)}). Useful for system-specific
#'   endpoint assumptions.
#' @param bulk_density Numeric of length one or nrow(soil_df). Soil bulk
#'   density (g cm\eqn{^{-3}}) for
#'   volumetric conversion of demand to mmol O\eqn{_2} L\eqn{^{-1}} porewater.
#'   Set to \code{NULL} to skip volumetric conversion.
#' @param theta_v Numeric of length one or nrow(soil_df), or \code{NULL}.
#'   Volumetric water content (L water L\eqn{^{-1}} bulk soil). When omitted
#'   and \code{bulk_density} is supplied, saturated porosity is estimated as
#'   \code{1 - bulk_density / particle_density}.
#' @param particle_density Numeric. Particle density in g cm\eqn{^{-3}} used
#'   only for the saturated-porosity estimate. Default 2.65.
#' @param return_components Logical. If \code{TRUE} (default), return a
#'   per-species contribution matrix.
#'
#' @details
#' \strong{Interpretation --- the 25-fold contrast.}
#'
#' The mechanistic review (Ghotbi \emph{et al.}, 2026) provides a worked
#' example: a Fe-rich rhizosphere containing 50 mmol Fe(II) kg\eqn{^{-1}},
#' 5 mmol FeS kg\eqn{^{-1}}, 2 mmol Mn(II) kg\eqn{^{-1}}, 2 mmol NH\eqn{_4^+}
#' kg\eqn{^{-1}}, 2 mmol acetate kg\eqn{^{-1}}, and 0.5 mmol CH\eqn{_4}
#' kg\eqn{^{-1}} has a complete-oxidation O\eqn{_2} ceiling of
#' \eqn{\approx 34} mmol O\eqn{_2} kg\eqn{^{-1}}, versus only \eqn{\approx 1.3}
#' mmol O\eqn{_2} kg\eqn{^{-1}} in an illustrative initial pore-gas stock:
#' about a 26-fold contrast. The assumed stock is not air-saturated porewater.
#' Continuing atmospheric and root O2 delivery can replenish it. Accessibility,
#' reaction kinetics and transport determine realized demand during an event.
#'
#' \strong{Pyrite stoichiometry.}
#'
#' Complete pyrite oxidation to sulfate and Fe(III) oxyhydroxide releases
#' 4 mol H\eqn{^+} mol\eqn{^{-1}} FeS\eqn{_2} and consumes 3.75 mol O\eqn{_2}:
#' FeS\eqn{_2} + 15/4 O\eqn{_2} + 7/2 H\eqn{_2}O →
#' Fe(OH)\eqn{_3} + 2H\eqn{_2}SO\eqn{_4}.
#' Partial oxidation to sulfur intermediates (S\eqn{^0}, thiosulfate) requires
#' fewer moles; adjust via \code{custom_coefs}.
#'
#' \strong{pH coupling.}
#'
#' Under the applicable rate law, homogeneous abiotic Fe(II) oxidation increases \eqn{\approx}100-fold per unit pH rise
#' (Stumm & Lee, 1961; Millero \emph{et al.}, 1987). The stoichiometric
#' demand computed here is for complete oxidation and is independent of pH,
#' but actual O\eqn{_2} consumption rates will be pH-modulated. This function
#' reports a stoichiometric potential demand, not a thermodynamic limit or rate.
#'
#' @return A list:
#' \describe{
#'   \item{\code{o2_demand}}{Numeric vector (mmol O\eqn{_2} kg\eqn{^{-1}
#'     dry soil}) of complete-oxidation O\eqn{_2} demand per sample.}
#'   \item{\code{o2_deficit_ratio}}{Per-sample \eqn{O_2^{\mathrm{demand}} /
#'     O_2^{\mathrm{supply}}}. Values \eqn{>1} indicate demand exceeds the specified stock.
#'     \code{NA} when \code{o2_supply_col} is absent.}
#'   \item{\code{o2_demand_vol}}{Volumetric O\eqn{_2} demand (mmol L\eqn{^{-1}}
#'     porewater); \code{NA} if \code{bulk_density} is \code{NULL}.}
#'   \item{\code{n_species_used}}{Integer. Number of reduced-pool columns found.}
#'   \item{\code{components}}{Data frame (one row per species) with:
#'     \code{species}, \code{stoich_coef}, \code{mean_inventory},
#'     \code{mean_o2_contribution}, \code{fraction_of_total_demand}.
#'     Returned only when \code{return_components = TRUE}.}
#'   \item{\code{stoich_table}}{Data frame of the full stoichiometric table
#'     used, including any custom overrides.}
#' }
#'
#' @references
#' Ghotbi, M., Ghotbi, M., Stukenbrock, E. H., Mühling, K. H., &
#' Spielvogel, S. (2026). Rhizosphere redox recovery after hydrological
#' disturbance: mechanisms across the soil--plant--microbiome continuum.
#' \emph{Manuscript submitted}.
#'
#' Stumm, W., & Lee, G. F. (1961). Oxygenation of ferrous iron.
#' \emph{Industrial & Engineering Chemistry}, 53, 143--146.
#'
#' Millero, F. J., Sotolongo, S., & Izaguirre, M. (1987). The oxidation kinetics
#' of Fe(II) in seawater. \emph{Geochimica et Cosmochimica Acta}, 51, 793--801.
#'
#' @seealso \code{\link{rri_accessible_capacity}},
#'   \code{\link{rri_capacity_index}}, \code{\link{rri_root_physio}}
#'
#' @examples
#' ## Reproduce the worked example from Ghotbi et al. (2026) Box 1
#' worked_example <- data.frame(
#'   Fe2 = 50.0, # mmol kg-1
#'   FeS = 5.0,
#'   Mn2 = 2.0,
#'   NH4 = 2.0,
#'   acetate = 2.0,
#'   CH4 = 0.5, # mmol kg-1; set ch4_unit = "umol_kg" for umol input
#'   O2_pw = 1.3 # assumed initial O2 stock, mmol/kg; NOT air-saturated porewater
#' )
#'
#' demand <- rri_o2_demand(
#'   soil_df = worked_example,
#'   fe2_col = "Fe2",
#'   mn2_col = "Mn2",
#'   fes_col = "FeS",
#'   nh4_col = "NH4",
#'   acetate_col = "acetate",
#'   ch4_col = "CH4",
#'   o2_supply_col = "O2_pw",
#'   return_components = TRUE
#' )
#'
#' demand$o2_demand # should be ~34 mmol O2 kg-1
#' demand$o2_deficit_ratio # should be ~26
#' demand$components
#'
#' @export
rri_o2_demand <- function(
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
) {
  if (!is.data.frame(soil_df)) soil_df <- as.data.frame(soil_df)
  ch4_unit <- match.arg(ch4_unit)
  acetate_basis <- match.arg(acetate_basis)

  ## -- stoichiometric coefficients (Box 1, Ghotbi et al. 2026) -----------------
  # O2 mol consumed per mol reduced species for complete oxidation
  stoich <- c(
    fe2_col     = 0.25, # Fe2+ → Fe(III)(OH)3
    mn2_col     = 0.50, # Mn2+ → MnO2
    hs_col      = 2.00, # HS- → SO4 2-
    fes_col     = 2.25, # FeS → Fe(III) + SO4 2-
    fes2_col    = 3.75, # FeS2 → Fe(III) + 2 SO4 2-
    nh4_col     = 2.00, # NH4+ → NO3-  (nitrification: 2 O2 mol-1)
    ch4_col     = 2.00, # CH4 → CO2 + H2O
    acetate_col = if (acetate_basis == "acetate") 2.00 else 1.00
  )

  if (!is.null(custom_coefs)) {
    if (!is.numeric(custom_coefs) || is.null(names(custom_coefs)) ||
        anyDuplicated(names(custom_coefs)) || anyNA(names(custom_coefs)))
      stop("custom_coefs must be uniquely named numeric coefficients.");
    bad <- setdiff(names(custom_coefs), names(stoich))
    if (length(bad) > 0) {
      warning("custom_coefs names not recognised and will be ignored: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
    stoich[intersect(names(custom_coefs), names(stoich))] <- custom_coefs[
      intersect(names(custom_coefs), names(stoich))
    ]
  }

  if (any(!is.finite(stoich)) || any(stoich < 0))
    stop("Stoichiometric coefficients must be finite and non-negative.")
  ## -- match column arguments to stoich names ------------------------------------
  col_args <- list(
    fe2_col     = fe2_col,
    mn2_col     = mn2_col,
    hs_col      = hs_col,
    fes_col     = fes_col,
    fes2_col    = fes2_col,
    nh4_col     = nh4_col,
    ch4_col     = ch4_col,
    acetate_col = acetate_col
  )

  specified <- Filter(Negate(is.null), col_args)
  if (any(!vapply(specified, function(x) is.character(x) && length(x) == 1L &&
                  !is.na(x) && nzchar(x), logical(1)))) stop("Species columns must be single names.")
  if (anyDuplicated(unlist(specified))) stop("A reservoir column cannot be counted twice.")
  numeric_cols <- intersect(c(unlist(specified), o2_supply_col), names(soil_df))
  soil_df[, numeric_cols] <- .rri_numeric_df(soil_df[, numeric_cols, drop = FALSE])
  present <- Filter(
    function(x) !is.null(x) && x %in% names(soil_df),
    col_args
  )

  if (length(present) == 0) {
    stop(
      "None of the specified reduced-species columns were found in `soil_df`. ",
      "Provide at least one of: fe2_col, mn2_col, hs_col, fes_col, fes2_col, ",
      "nh4_col, ch4_col, acetate_col.",
      call. = FALSE
    )
  }

  missing_cols <- Filter(function(x) !is.null(x) && !x %in% names(soil_df), col_args)
  if (length(missing_cols) > 0) {
    warning("The following columns were not found and will be skipped: ",
      paste(unlist(missing_cols), collapse = ", "),
      call. = FALSE
    )
  }

  n <- nrow(soil_df)
  arg_names <- names(present)
  n_sp <- length(present)
  demand_mat <- matrix(NA_real_, nrow = n, ncol = n_sp)
  inv_mat <- matrix(NA_real_, nrow = n, ncol = n_sp)
  sp_labels <- character(n_sp)
  coefs <- numeric(n_sp)

  for (i in seq_len(n_sp)) {
    aname <- arg_names[i]
    col_name <- present[[i]]
    inv_i <- as.numeric(soil_df[[col_name]])
    inv_i[!is.finite(inv_i) | inv_i < 0] <- NA_real_

    if (aname == "ch4_col" && ch4_unit == "umol_kg") inv_i <- inv_i / 1000

    coef_i <- stoich[aname]
    demand_mat[, i] <- inv_i * coef_i
    inv_mat[, i] <- inv_i
    sp_labels[i] <- col_name
    coefs[i] <- coef_i
  }

  o2_demand <- rowSums(demand_mat, na.rm = TRUE)
  o2_demand[rowSums(is.finite(demand_mat)) == 0L] <- NA_real_

  ## -- deficit ratio against porewater O2 supply --------------------------------
  o2_deficit_ratio <- rep(NA_real_, n)
  if (!is.null(o2_supply_col) && o2_supply_col %in% names(soil_df)) {
    o2_supply <- as.numeric(soil_df[[o2_supply_col]])
    o2_supply[!is.finite(o2_supply) | o2_supply < 0] <- NA_real_
    o2_deficit_ratio <- o2_demand / o2_supply
    o2_deficit_ratio[is.nan(o2_deficit_ratio)] <- NA_real_
  } else if (!is.null(o2_supply_col)) {
    warning("`o2_supply_col` '", o2_supply_col,
      "' not found in soil_df; deficit ratio not computed.",
      call. = FALSE
    )
  }

  ## -- volumetric conversion (mmol O2 L-1 porewater) ----------------------------
  o2_demand_vol <- rep(NA_real_, n)
  if (!is.null(bulk_density)) {
    if (!is.numeric(bulk_density) || !length(bulk_density) %in% c(1L, n)) {
      stop("`bulk_density` must be NULL or numeric with length 1 or nrow(soil_df).",
        call. = FALSE
      )
    }
    bd <- if (length(bulk_density) == 1L) {
      rep(bulk_density, n)
    } else {
      as.numeric(bulk_density)
    }
    if (any(!is.finite(bd)) || any(bd <= 0)) {
      stop("All `bulk_density` values must be finite and > 0.", call. = FALSE)
    }
    if (!is.numeric(particle_density) || length(particle_density) != 1L ||
      !is.finite(particle_density) || particle_density <= max(bd)) {
      stop("`particle_density` must be a finite scalar greater than all bulk densities.",
        call. = FALSE
      )
    }
    theta <- if (is.null(theta_v)) {
      1 - bd / particle_density
    } else if (is.numeric(theta_v) && length(theta_v) %in% c(1L, n)) {
      if (length(theta_v) == 1L) rep(theta_v, n) else as.numeric(theta_v)
    } else {
      stop("`theta_v` must be NULL or numeric with length 1 or nrow(soil_df).",
        call. = FALSE
      )
    }
    if (any(!is.finite(theta)) || any(theta <= 0 | theta > 1)) {
      stop("All `theta_v` values must be finite and in (0, 1].", call. = FALSE)
    }

    # bulk_density in g cm-3 is numerically kg dry soil L-1 bulk soil.
    # Dividing by L water L-1 bulk soil gives mmol O2 L-1 porewater.
    if (any(theta > 1 - bd / particle_density + 1e-10))
      stop("theta_v cannot exceed the porosity implied by the supplied densities.")
    o2_demand_vol <- o2_demand * bd / theta
  }

  ## -- components table ---------------------------------------------------------
  mean_demand <- vapply(seq_len(ncol(demand_mat)),
                        function(j) .rri_mean(demand_mat[, j]), numeric(1))
  total_mean <- sum(mean_demand, na.rm = TRUE)
  frac_demand <- if (total_mean > 0) mean_demand / total_mean else rep(NA_real_, n_sp)

  stoich_table <- data.frame(
    species_arg = arg_names,
    column_used = sp_labels,
    stoich_coef_O2 = coefs,
    endpoint = c(
      fe2_col     = "Fe(III) oxyhydroxide",
      mn2_col     = "MnO2",
      hs_col      = "SO4 2-",
      fes_col     = "Fe(III) + SO4 2-",
      fes2_col    = "Fe(III) + 2 SO4 2-",
      nh4_col     = "NO3-",
      ch4_col     = "CO2",
      acetate_col = "CO2"
    )[arg_names],
    stringsAsFactors = FALSE
  )

  out <- list(
    o2_demand = o2_demand,
    o2_deficit_ratio = o2_deficit_ratio,
    o2_demand_vol = o2_demand_vol,
    n_species_used = n_sp,
    n_species_observed = rowSums(is.finite(demand_mat)),
    species_coverage = rowSums(is.finite(demand_mat)) / sum(!vapply(col_args, is.null, logical(1))),
    interpretation = "Stoichiometric demand subtotal; supply ratio uses a stock, not an oxygen delivery rate",
    ch4_unit_used = ch4_unit,
    acetate_basis_used = acetate_basis,
    stoich_table = stoich_table
  )

  if (return_components) {
    out$components <- data.frame(
      species               = sp_labels,
      n_observed            = colSums(is.finite(inv_mat)),
      stoich_coef           = coefs,
      mean_inventory_mmol   = round(vapply(seq_len(ncol(inv_mat)),
        function(j) .rri_mean(inv_mat[, j]), numeric(1)), 3),
      mean_o2_contribution  = round(mean_demand, 3),
      fraction_total_demand = round(frac_demand, 4),
      stringsAsFactors      = FALSE
    )
  }

  out
}
