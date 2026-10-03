#' Aligned score trajectories and recovery-descriptor availability
#'
#' A two-panel display separating score dynamics from finite-output counts.
#' A finite result is computational availability, not evidence of ecological
#' validity, adequate precision or successful return to baseline.
#' @param res,id Pipeline scores and their observation identifiers.
#' @param rec Recovery table with one row per trajectory in id.
#' @param time_col Name of the numeric time column.
#' @param group_cols Columns identifying trajectories in both id and rec.
#' @param perturb_start,perturb_end Optional exact event endpoints.
#' @param base_size Base font size in points.
#' @return A patchwork object, or a named list of ggplots without patchwork.
#'   Export at approximately 8 by 4.6 inches.
#' @examples
#' id <- data.frame(plot = rep(c("P1", "P2"), each=3), time=rep(1:3,2))
#' res <- list(row_scores=data.frame(RRI=c(.6,.3,.5,.7,.4,.6)))
#' rec <- data.frame(plot=c("P1","P2"), k_recovery=c(.1,NA_real_))
#' plot_rri_recovery_diagnostics(res,id,rec,group_cols="plot")
#' @export
plot_rri_recovery_diagnostics <- function(res,id,rec,time_col="time",
  group_cols=c("plot","depth","plant_id"),perturb_start=NULL,perturb_end=NULL,base_size=9) {
  id <- as.data.frame(id); rec <- as.data.frame(rec)
  if(!all(group_cols %in% names(rec))||!all(group_cols %in% names(id)))
    stop("Both id and rec require all group_cols.",call.=FALSE)
  akey <- unique(.rri_key(id,group_cols)); bkey <- .rri_key(rec,group_cols)
  if(!nrow(rec)||anyDuplicated(bkey)||!setequal(akey,bkey))
    stop("rec must have exactly one row per trajectory in id.",call.=FALSE)
  a <- plot_rri_recovery_map(res,id,time_col=time_col,group_cols=group_cols,
    perturb_start=perturb_start,perturb_end=perturb_end,palette="viridis",
    base_size=base_size,max_groups=length(akey),order_by="key",direction=1,tile_gap=0,event_colour="white")+
    ggplot2::labs(title="a  Aligned score trajectories",subtitle=NULL,caption=NULL,
      x=paste0(time_col," (input units)"),fill="HRRI score (0-1)")+
    theme_ems(base_size)+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),
      axis.line=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_text(size=base_size-2),legend.position="bottom")+
    ggplot2::scale_y_discrete(labels=function(x) gsub("Plant","",gsub(" \\| "," / ",x)))+
    ggplot2::guides(fill=.hrri_vector_colourbar(title="HRRI score (0-1)",title.position="bottom",barwidth=grid::unit(28,"mm"),barheight=grid::unit(2,"mm")))
  metrics <- c("depth_min_frac","tau_lag","k_recovery","t_half","overshoot_frac",
               "H_hysteresis","temporal_asymmetry","incomplete_return_frac")
  labels <- c("Maximum decline","Response lag","Recovery rate","Half-time","Overshoot",
              "Loop area","Deficit asymmetry","Terminal displacement")
  counts <- vapply(metrics,function(m) {
    if(!m %in% names(rec)) return(0L)
    if(!is.numeric(rec[[m]])) stop("Recovery metric columns must be numeric.",call.=FALSE)
    sum(is.finite(rec[[m]]))
  },integer(1))
  d <- data.frame(metric=metrics,label=factor(labels,levels=rev(labels)),n=counts,total=nrow(rec),present=metrics %in% names(rec))
  d$text <- ifelse(d$present,paste0(d$n,"/",d$total),"not supplied")
  b <- ggplot2::ggplot(d,ggplot2::aes(y=.data$label))+
    ggplot2::geom_col(ggplot2::aes(x=.data$total),fill="#edf0f2",width=.62)+
    ggplot2::geom_col(ggplot2::aes(x=.data$n),fill="#167c80",width=.62)+
    ggplot2::geom_text(ggplot2::aes(x=.data$n+.025*.data$total,label=.data$text),hjust=0,size=base_size/3.2)+
    ggplot2::scale_x_continuous(limits=c(0,nrow(rec)*1.3),breaks=pretty(c(0,nrow(rec))),expand=c(0,0))+
    ggplot2::labs(title="b  Descriptor availability",x="Trajectories with finite values",y=NULL)+
    theme_ems(base_size)+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_blank(),axis.line.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
  .hrri_assemble(list(trajectories=a,availability=b),ncol=2,widths=c(1.2,1))
}
