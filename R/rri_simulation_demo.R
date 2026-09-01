#' Reproducible observable-only HRRI demonstration
#'
#' @description Generates both forcing scenarios with the package simulator,
#' calibrates illustrative targets using a separate baseline simulation, and
#' compares full, 4-soil/3-plant/2-microbe, and single-domain panels with unchanged
#' feature targets and tolerances. No hidden capacity, alpha, k, memory or
#' latent_truth column enters the observed-feature index. This is a software
#' demonstration, not a validation of latent-state recovery or ecological prediction.
#' @param seed Non-negative integer seed.
#' @return Simulation objects, explicit reference specifications and all plotted
#' tables. The capacity-horizon table is an internal equation check using known
#' synthetic parameters; it is deliberately separate from observable scoring.
#' @importFrom stats median aggregate
#' @examples
#' \dontrun{
#'   demo <- rri_simulation_demo(seed = 20260830L)
#'   head(demo$scores)
#' }
#' @export
rri_simulation_demo <- function(seed = 20260830L) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max - 2000L || seed != floor(seed)) stop("Invalid seed.")
  scenarios <- c("flood_drain", "drought_rewet")
  features <- c("Eh", "FeII_mmol_kg", "MnII_mmol_kg", "NO3_mmol_kg", "NH4_mmol_kg",
    "porewater_O2_mmol_L", "SPAD", "FvFm", "ROS_load", "ROL", "log10_mtrA", "log10_amoA_AOA", "log10_nosZ")
  domains <- c(rep("Soil", 6), rep("Physio", 4), rep("Micro", 3))
  tolerances <- c(100, 20, 6, 3, 5, 0.08, 8, 0.10, 0.5, 0.5, 1, 1, 1)
  extract <- function(sim) {
    soil <- sim$soil_data[, features[1:6], drop = FALSE]
    plant <- sim$plant_data[, features[7:10], drop = FALSE]
    micro <- log10(1 + sim$micro_gene_abundance[, c("mtrA", "amoA_AOA", "nosZ")])
    names(micro) <- features[11:13]
    cbind(soil, plant, micro)
  }
  panels <- list(full = features,
    soil4_plant3_micro2 = features[c(1,2,3,5,7,8,9,11,12)],
    soil_only = features[1:6], plant_only = features[7:10], microbe_only = features[11:13])
  outputs <- lapply(seq_along(scenarios), function(j) {
    args <- list(n_plot=2, n_depth=2, n_plant=2, n_time=60, p_micro=12,
      scenario=scenarios[j], n_cycles=1L, disturbance_center=22,
      disturbance_width=0.035, seasonal_amp=0, MNAR_strength=0.15)
    sim <- do.call(simulate_redox_holobiont, c(args, list(seed=seed+j)))
    calibration <- do.call(simulate_redox_holobiont, c(args, list(seed=seed+1000L+j)))
    x <- extract(sim); xc <- extract(calibration)
    pre <- calibration$id$time <= 5
    target <- vapply(xc[pre, , drop=FALSE], stats::median, numeric(1), na.rm=TRUE)
    reference <- data.frame(feature=features, domain=domains, target=unname(target),
      tolerance=tolerances, weight=1, stringsAsFactors=FALSE)
    fitted <- lapply(panels, function(cols) rri_reference_scores(x[,cols,drop=FALSE],
      reference, id=sim$id, min_coverage=0.5))
    score_rows <- do.call(rbind, lapply(names(fitted), function(nm) {
      z <- fitted[[nm]]$row_scores
      data.frame(scenario=scenarios[j], panel=nm, unit_id=z$unit_id, time=z$time,
        RRI=z$RRI, domain_coverage=z$domain_coverage, n_domains=z$n_domains)
    }))
    plotted <- data.frame(scenario=scenarios[j], unit_id=sim$id$unit_id,
      time=sim$id$time, WFPS=sim$forcing$WFPS, Eh=x$Eh, SPAD=x$SPAD, log10_mtrA=x$log10_mtrA)
    # Equation demonstration: each horizon applied to one fixed simulated state.
    row <- which(sim$id$time == 22)[1]
    horizons <- c(0,1,3,6,12,24,48,96)
    reservoirs <- list(EAC=list(Q_col="EAC",alpha="alpha_accept",k="k_accept_h",type="EAC"),
      EDC=list(Q_col="EDC",alpha="alpha_donate",k="k_donate_h",type="EDC"))
    ca <- rri_accessible_capacity(sim$soil_data[rep(row,length(horizons)),], reservoirs,
      tau=horizons,normalise=FALSE)
    capacity <- data.frame(scenario=scenarios[j],tau_h=horizons,
      accessible_EAC=ca$cacc_eac, accessible_EDC=ca$cacc_edc)
    list(simulation=sim, reference=reference, fitted=fitted, scores=score_rows,
      observations=plotted, capacity=capacity,
      calibration_seed=seed+1000L+j, analysis_seed=seed+j)
  })
  names(outputs) <- scenarios
  list(scenarios=outputs, scores=do.call(rbind,lapply(outputs,`[[`,"scores")),
    observations=do.call(rbind,lapply(outputs,`[[`,"observations")),
    capacity=do.call(rbind,lapply(outputs,`[[`,"capacity")),
    metadata=list(seed=seed, panels=panels, calibration_window="days 1-5 in separate simulation",
      tolerance_status="Illustrative declared deviations, not empirically validated thresholds",
      scoring="Fixed target proximity, not universal resilience or estimated latent truth",
      input_domains="DNA represents functional potential, not microbial process activity"))
}

#' Draw the reproducible simulation demonstration
#' @param demo Result of rri_simulation_demo.
#' @param figure observations, coverage or capacity.
#' @return Invisibly returns the plotted data table; draws on the active device.
#' Uses base graphics, so no optional plotting package is required.
#' @importFrom graphics par plot abline mtext lines legend matplot
#' @importFrom stats aggregate
#' @examples
#' \dontrun{
#'   demo <- rri_simulation_demo(seed = 20260830L)
#'   plot_rri_simulation_demo(demo, figure = "observations")
#' }
#' @export
plot_rri_simulation_demo <- function(demo, figure=c("observations","coverage","capacity")) {
  figure <- match.arg(figure); old <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old),add=TRUE)
  scenarios <- names(demo$scenarios)
  if (figure=="observations") {
    d <- demo$observations
    vars <- c("WFPS","Eh","SPAD","log10_mtrA")
    labels <- c("Water-filled pore space (fraction)","Eh (mV)","SPAD index","log10(1 + mtrA copies/g)")
    graphics::par(mfrow=c(4,2),mar=c(3.2,4.2,2.4,1),oma=c(1,0,2,0),mgp=c(2.5,0.7,0))
    for (i in seq_along(vars)) for (sc in scenarios) {
      z <- d[d$scenario==sc,]; a <- stats::aggregate(z[[vars[i]]],list(time=z$time),.rri_mean)
      graphics::plot(a$time,a$x,type="l",lwd=2,col=if(sc==scenarios[1]) "#0072B2" else "#D55E00",
        xlab="Day",ylab=labels[i],ylim=range(d[[vars[i]]],na.rm=TRUE),main=gsub("_"," / ",sc),cex.main=0.9)
      graphics::abline(v=22,lty=3,col="grey60")
    }
    graphics::mtext("Synthetic observables: mean across simulated units",outer=TRUE,font=2,cex=1.1)
  } else if (figure=="coverage") {
    d <- demo$scores; panels <- names(demo$metadata$panels)
    cols <- c("#000000","#0072B2","#D55E00","#009E73","#CC79A7")
    graphics::par(mfrow=c(1,2),mar=c(4.2,4.3,3,1),oma=c(6,0,2,0))
    for (sc in scenarios) {
      graphics::plot(c(1,60),c(0,1),type="n",xlab="Day",ylab="Reference proximity (0-1)",
        main=gsub("_"," / ",sc))
      for (j in seq_along(panels)) {
        z <- d[d$scenario==sc & d$panel==panels[j],]
        a <- stats::aggregate(z$RRI,list(time=z$time),.rri_mean)
        graphics::lines(a$time,a$x,col=cols[j],lwd=2,lty=j)
      }
    }
    graphics::mtext("Unchanged feature scales; domain subsets measure different summaries",outer=TRUE,font=2,cex=1)
    graphics::par(fig=c(0,1,0,0.13),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0));graphics::plot.new()
    graphics::legend("center",legend=c("Full: 6 soil / 4 plant / 3 microbe","Sparse: 4 / 3 / 2",
      "Soil only","Plant only","Microbe only"),col=cols,lty=1:5,lwd=2,bty="n",ncol=3,cex=0.75)
  } else {
    d <- demo$capacity
    graphics::par(mfrow=c(1,2),mar=c(4.2,4.3,3.5,1),oma=c(1,0,2,0))
    for(sc in scenarios) {
      z <- d[d$scenario==sc,]
      graphics::matplot(z$tau_h,z[,c("accessible_EAC","accessible_EDC")],type="l",lty=1:2,
        lwd=2,col=c("#0072B2","#D55E00"),xlab="Access window (h)",
        ylab="Accessible capacity (mmol e-/kg)",main=gsub("_"," / ",sc))
      graphics::legend("bottomright",c("EAC","EDC"),col=c("#0072B2","#D55E00"),lty=1:2,lwd=2,bty="n")
    }
    graphics::mtext("Internal equation check with known synthetic parameters",outer=TRUE,font=2,cex=1.1)
  }
  invisible(d)
}
