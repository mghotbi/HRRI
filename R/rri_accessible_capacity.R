#' Event-window accessible reservoir capacities
#'
#' @description Computes Q * alpha * (1-exp(-k*tau)) for declared reservoirs.
#' Q must be nonoverlapping electron-equivalent inventories with explicit reaction
#' endpoints. EAC and EDC are returned separately. Their sum is an inventory
#' descriptor, not electron flux; their difference is not an oxygen budget.
#' @param soil_df Numeric reservoir measurements.
#' @param reservoirs Nonempty named list of Q_col, alpha, k and type (EAC/EDC).
#' alpha in `[0, 1]` and k >= 0 may be scalars, row vectors or column names. Parameters
#' must be specified for the relevant process and conditions; they are not
#' identifiable separately from one accessible-capacity observation.
#' @param tau Non-negative duration, scalar or row vector; units reciprocal to k.
#' @param normalise Divide by the sum of observed inventories. This produces
#' an accessible fraction, not absolute capacity or guaranteed comparability.
#' @param return_components Include reservoir contribution summaries.
#' @return Capacities, observed subtotal, fraction and reservoir coverage.
#' Missing reservoir types remain NA. Partial rows are labelled observed subtotals;
#' absence is not zero. ck_limited is retained as NA because 0.30 is not a
#' validated threshold. Negative inventories are treated as missing.
#' @details Default reservoir parameters are illustrative scenario values only.
#' Fe(II) oxidation rates must not be assigned as Fe(III) reduction constants.
#' For capacity estimation use experimentally constrained process-specific rates.
#' @examples
#' df <- data.frame(EAC = c(10, 20, 30), EDC = c(5, 8, 12))
#' res <- rri_accessible_capacity(
#'   df,
#'   reservoirs = list(
#'     bulk_EAC = list(Q_col = "EAC", alpha = 0.5, k = 0.2, type = "EAC"),
#'     bulk_EDC = list(Q_col = "EDC", alpha = 0.45, k = 0.15, type = "EDC")
#'   ),
#'   tau = 24
#' )
#' res$cacc
#' @export
rri_accessible_capacity <- function(
  soil_df,
  reservoirs,
  tau = 24,
  normalise = TRUE,
  return_components = FALSE
) {
  ## -- input checks ------------------------------------------------------------
  if (!is.data.frame(soil_df)) soil_df <- as.data.frame(soil_df)
  if (!is.list(reservoirs) || length(reservoirs) == 0) {
    stop("`reservoirs` must be a non-empty named list.", call. = FALSE)
  }

  if (is.null(names(reservoirs)) || any(!nzchar(names(reservoirs))) || anyDuplicated(names(reservoirs)))
    stop("reservoirs must have unique nonempty names.", call. = FALSE)
  n <- nrow(soil_df)
  if (n < 1L || anyDuplicated(names(soil_df)))
    stop("soil_df must have rows and unique column names.")
  if (!is.numeric(tau) || !length(tau) %in% c(1L, n) ||
    any(!is.finite(tau)) || any(tau < 0)) {
    stop("`tau` must be non-negative and have length 1 or nrow(soil_df).",
      call. = FALSE
    )
  }
  tau_vec <- if (length(tau) == 1L) rep(as.numeric(tau), n) else as.numeric(tau)

  # alpha and k may be fixed reservoir parameters, row-wise numeric vectors,
  # or names of columns in soil_df. This is required when connectivity and
  # exchange kinetics are hidden states that change through time.
  resolve_parameter <- function(x, field, reservoir_name) {
    if (is.character(x) && length(x) == 1L && x %in% names(soil_df)) {
      if (!is.numeric(soil_df[[x]]))
        stop("Reservoir '", reservoir_name, "': parameter column '", x,
             "' must be numeric.")
      out <- soil_df[[x]]
    } else if (is.numeric(x) && length(x) %in% c(1L, n)) {
      out <- if (length(x) == 1L) rep(as.numeric(x), n) else as.numeric(x)
    } else {
      stop("Reservoir '", reservoir_name, "': `", field,
        "` must be a scalar, a length-n vector, or a soil_df column name.",
        call. = FALSE
      )
    }
    if (any(!is.finite(out))) {
      stop("Reservoir '", reservoir_name, "': `", field,
        "` contains non-finite values.",
        call. = FALSE
      )
    }
    out
  }

  ## -- resolve which reservoirs are present in the data ------------------------
  for (nm in names(reservoirs)) {
    r <- reservoirs[[nm]]
    if (!is.list(r) || !all(c("Q_col","alpha","k","type") %in% names(r)) ||
        !is.character(r$Q_col) || length(r$Q_col) != 1L || is.na(r$Q_col) ||
        !is.character(r$type) || length(r$type) != 1L ||
        !r$type %in% c("EAC","EDC"))
      stop("Invalid reservoir specification: ", nm)
  }
  present <- Filter(function(r) r$Q_col %in% names(soil_df), reservoirs)

  if (length(present) == 0) {
    stop(
      "No requested reservoir inventory columns were found in soil_df.",
      call. = FALSE
    )
  }

  q_cols <- vapply(present, function(r) r$Q_col, character(1))
  if (anyDuplicated(q_cols)) stop("The same inventory column cannot be counted twice.")
  missing_res <- setdiff(names(reservoirs), names(present))
  if (length(missing_res) > 0) {
    warning(
      "The following reservoirs were skipped because Q_col was absent: ",
      paste(missing_res, collapse = ", "),
      call. = FALSE
    )
  }

  ## -- compute per-reservoir contributions -------------------------------------
  # Cacc = sum_i Q_i * alpha_i * (1 - exp(-k_i * tau)); parameters and tau
  # may be fixed or row specific, as declared by the caller.

  reservoir_names <- names(present)
  n_res <- length(present)
  contribution_mat <- matrix(NA_real_, nrow = n, ncol = n_res)
  total_Q_mat <- matrix(NA_real_, nrow = n, ncol = n_res)

  alpha_vec <- numeric(n_res)
  k_vec <- numeric(n_res)
  type_vec <- character(n_res)
  sat_frac <- numeric(n_res) # (1 - exp(-k * tau))
  mean_Q <- numeric(n_res)
  mean_contr <- numeric(n_res)

  for (i in seq_len(n_res)) {
    r <- present[[i]]
    if (!is.numeric(soil_df[[r$Q_col]]))
      stop("Reservoir '", reservoir_names[i], "': Q_col must be numeric.")
    Q_i <- soil_df[[r$Q_col]]
    Q_i[!is.finite(Q_i) | Q_i < 0] <- NA_real_

    alpha_i <- resolve_parameter(r$alpha, "alpha", reservoir_names[i])
    k_i <- resolve_parameter(r$k, "k", reservoir_names[i])
    if (any(alpha_i < 0 | alpha_i > 1)) {
      stop("Reservoir '", reservoir_names[i],
        "': all alpha values must be in [0, 1].",
        call. = FALSE
      )
    }
    if (any(k_i < 0)) {
      stop("Reservoir '", reservoir_names[i],
        "': all k values must be >= 0.",
        call. = FALSE
      )
    }

    sf_i <- -expm1(-k_i * tau_vec)

    contrib_i <- Q_i * alpha_i * sf_i
    contribution_mat[, i] <- contrib_i
    total_Q_mat[, i] <- Q_i

    alpha_vec[i] <- .rri_mean(alpha_i)
    k_vec[i] <- .rri_mean(k_i)
    type_vec[i] <- r$type
    sat_frac[i] <- .rri_mean(sf_i)
    mean_Q[i] <- .rri_mean(Q_i)
    mean_contr[i] <- .rri_mean(contrib_i)
  }

  ## -- aggregate to Cacc --------------------------------------------------------
  cacc_raw <- rowSums(contribution_mat, na.rm = TRUE)
  total_Q <- rowSums(total_Q_mat, na.rm = TRUE)
  all_missing <- rowSums(is.finite(total_Q_mat)) == 0L
  cacc_raw[all_missing] <- NA_real_
  total_Q[all_missing] <- NA_real_

  # EAC and EDC components separately (before normalisation)
  is_eac <- type_vec == "EAC"
  sum_type <- function(j) {
    if (!any(j)) return(rep(NA_real_, n))
    m <- contribution_mat[, j, drop = FALSE]
    z <- rowSums(m, na.rm = TRUE)
    z[rowSums(is.finite(m)) == 0] <- NA_real_
    z
  }
  cacc_eac <- sum_type(is_eac)
  cacc_edc <- sum_type(!is_eac)
  net_balance <- cacc_eac - cacc_edc

  # Fraction of total inventory that is accessible
  cacc_fraction <- ifelse(total_Q > 0, cacc_raw / total_Q, NA_real_)

  # Theoretical maximum (all alpha = 1, all sat = 1)
  cacc_max <- rowSums(total_Q_mat, na.rm = TRUE)
  cacc_max[cacc_max == 0] <- NA_real_

  if (normalise) {
    cacc_out <- cacc_raw / pmax(cacc_max, .Machine$double.eps)
    cacc_out <- pmax(0, pmin(1, cacc_out))
  } else {
    cacc_out <- cacc_raw
  }

  ## -- connectivity/kinetics limitation flag -----------------------------------
  ck_limited <- rep(NA, n) # no empirically justified diagnostic cutoff

  ## -- build output -------------------------------------------------------------
  out <- list(
    cacc                  = cacc_out,
    cacc_raw              = cacc_raw,
    cacc_eac              = cacc_eac,
    cacc_edc              = cacc_edc,
    net_oxidative_balance = net_balance,
    cacc_fraction         = cacc_fraction,
    total_inventory       = total_Q,
    n_reservoirs_observed = rowSums(is.finite(total_Q_mat)),
    n_reservoirs_requested = length(reservoirs),
    reservoir_coverage = rowSums(is.finite(total_Q_mat)) / length(reservoirs),
    interpretation = "Observed reservoir subtotal; EAC and EDC are distinct capacities, not a transfer flux",
    ck_limited            = ck_limited,
    tau_used              = if (length(tau) == 1L) tau else tau_vec
  )

  if (return_components) {
    out$components <- data.frame(
      reservoir           = reservoir_names,
      type                = type_vec,
      alpha               = round(alpha_vec, 4),
      k_per_tau_unit      = round(k_vec, 4),
      saturation_fraction = round(sat_frac, 4),
      mean_Q              = round(mean_Q, 3),
      mean_contribution   = round(mean_contr, 3),
      stringsAsFactors    = FALSE
    )
  }

  out
}


#' Illustrative reservoir parameter template
#' @description Returns example values, not calibrated mineral-specific constants.
#' The default uses the core bulk EAC and EDC columns. Supply phase-resolved
#' column names only when those inventories are nonoverlapping and expressed in
#' electron-equivalent units.
#' @param eac_ferrihydrite_col,eac_goethite_col,eac_structural_col EAC phase columns.
#' @param edc_humic_fast_col,edc_humic_slow_col EDC fraction columns.
#' @param edc_eac_col,edc_edc_col Bulk EAC and EDC fallback columns. Argument
#' names are retained for backward compatibility.
#' @return Named list for rri_accessible_capacity; k uses inverse hours.
#' @examples
#' rri_default_reservoirs()
#' @export
rri_default_reservoirs <- function(
  eac_ferrihydrite_col = NULL,
  eac_goethite_col = NULL,
  eac_structural_col = NULL,
  edc_humic_fast_col = NULL,
  edc_humic_slow_col = NULL,
  edc_eac_col = "EAC", # bulk fallback
  edc_edc_col = "EDC" # bulk fallback
) {
  res <- list()

  ## EAC reservoirs --- ordered by decreasing k (most reactive first)
  if (!is.null(eac_ferrihydrite_col)) {
    res$ferrihydrite <- list(
      Q_col = eac_ferrihydrite_col,
      alpha = 0.85, # high surface area, well-connected porewater
      k     = 1.00, # illustrative h-1; not inferred from Fe(II) oxidation
      type  = "EAC"
    )
  }

  if (!is.null(eac_goethite_col)) {
    res$goethite <- list(
      Q_col = eac_goethite_col,
      alpha = 0.25, # moderate surface area
      k     = 0.04, # h-1: half-time ~17 h
      type  = "EAC"
    )
  }

  if (!is.null(eac_structural_col)) {
    res$structural_Fe <- list(
      Q_col = eac_structural_col,
      alpha = 0.08, # low connectivity (clay interlayer diffusion)
      k     = 0.01, # h-1: very slow exchange
      type  = "EAC"
    )
  }

  ## Bulk EAC fallback (if phase-specific columns unavailable)
  if (!is.null(edc_eac_col) && all(vapply(list(eac_ferrihydrite_col, eac_goethite_col, eac_structural_col), is.null, logical(1)))) {
    res$bulk_EAC <- list(
      Q_col = edc_eac_col,
      alpha = 0.50, # conservative estimate for mixed mineralogy
      k     = 0.20, # h-1
      type  = "EAC"
    )
  }

  ## EDC reservoirs
  if (!is.null(edc_humic_fast_col)) {
    res$humic_fast <- list(
      Q_col = edc_humic_fast_col,
      alpha = 0.80, # fast-exchanging quinone/phenol groups
      k     = 2.00, # h-1
      type  = "EDC"
    )
  }

  if (!is.null(edc_humic_slow_col)) {
    res$humic_slow <- list(
      Q_col = edc_humic_slow_col,
      alpha = 0.10, # condensed aromatic cores, slow access
      k     = 0.02, # h-1
      type  = "EDC"
    )
  }

  ## Bulk EDC fallback
  if (!is.null(edc_edc_col) && is.null(edc_humic_fast_col) && is.null(edc_humic_slow_col)) {
    res$bulk_EDC <- list(
      Q_col = edc_edc_col,
      alpha = 0.45,
      k     = 0.15,
      type  = "EDC"
    )
  }

  res
}
