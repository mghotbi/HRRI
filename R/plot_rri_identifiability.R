#' Identifiability, sensitivity and oxygen-accounting illustrations
#'
#' Four analytical panels for a single declared reservoir. They demonstrate
#' non-identifiability from one observation window and local conditioning of
#' two-window inversion; they do not estimate parameters from HRRI scores.
#' @param Q Positive inventory, in electron-equivalent units.
#' @param alpha Accessibility in (0, 1).
#' @param k,tau Positive reference rate and duration in reciprocal units.
#' @param duration_pairs Two-column numeric matrix; each row contains positive
#'   increasing observation durations.
#' @param fe_inventory Positive Fe(II) inventory in mmol per kg.
#' @param base_size Base font size in points.
#' @return A patchwork object, or a named list of ggplots if patchwork is absent.
#'   Export at approximately 7.4 by 6.4 inches.
#' @details The illustration uses C = Q alpha (1 - exp(-k tau)). Rate elasticity
#'   is k tau / expm1(k tau). Two-window sensitivity is the reciprocal absolute
#'   difference of those elasticities, not a confidence interval. Oxygen demand
#'   uses 0.25 mol O2 per mol Fe(II), and 2.25 mol O2 per mol FeS, for oxidation
#'   to Fe(III) and sulfate. These specified endpoints are essential.
#' @examples
#' plot_rri_identifiability()
#' @export
plot_rri_identifiability <- function(Q=100,alpha=.4,k=.15,tau=10,
  duration_pairs=rbind(c(1,40),c(5,20),c(2,4),c(20,40),c(10,12),c(40,80)),
  fe_inventory=50,base_size=9) {
  for(v in list(Q,k,tau,fe_inventory)) if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v<=0)
    stop("Inventories, k and tau must be positive finite scalars.",call.=FALSE)
  if(!is.numeric(alpha)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)
    stop("alpha must lie strictly between zero and one.",call.=FALSE)
  dp <- as.matrix(duration_pairs)
  if(!is.numeric(dp)||ncol(dp)!=2L||!nrow(dp)||any(!is.finite(dp))||any(dp<=0)||any(dp[,1]>=dp[,2]))
    stop("duration_pairs must have two positive increasing columns.",call.=FALSE)
  n <- "#274b6a"; te <- "#167c80"; ru <- "#c6633b"; go <- "#b58130"
  cc <- alpha*(-expm1(-k*tau)); km <- -log1p(-cc)/tau
  kg <- seq(km,max(k*10,km*1.1),length.out=800)
  d <- data.frame(k=kg,alpha=cc/(-expm1(-kg*tau)))
  a <- ggplot2::ggplot(d,ggplot2::aes(x=.data$k,y=.data$alpha))+
    ggplot2::geom_line(colour=n,linewidth=.8)+ggplot2::geom_hline(yintercept=cc,colour="#87939c",linetype=2)+
    ggplot2::annotate("point",x=k,y=alpha,colour=ru,size=2.5)+
    ggplot2::annotate("text",x=max(kg)*.45,y=(1+alpha)/2,label=sprintf("Reference (%.2f, %.2f)",k,alpha),size=2.7)+
    ggplot2::annotate("segment",x=max(kg)*.35,y=(1+alpha)/2-.04,xend=k,yend=alpha,colour="#87939c")+
    ggplot2::annotate("text",x=max(kg)*.65,y=cc+.06,label=sprintf("Asymptote = %.5f",cc),size=2.7)+
    ggplot2::labs(title="a  One window leaves a curve",x=expression(Exchange~rate~k~(time^{-1})),y=expression(Accessibility~alpha))+
    theme_ems(base_size)
  u <- exp(seq(log(.01),log(20),length.out=500)); elast <- function(x) x/expm1(x)
  b <- ggplot2::ggplot(data.frame(u=u,e=elast(u)),ggplot2::aes(x=.data$u,y=.data$e))+
    ggplot2::annotate("rect",xmin=.01,xmax=.1,ymin=0,ymax=Inf,fill=te,alpha=.05)+
    ggplot2::annotate("rect",xmin=3,xmax=20,ymin=0,ymax=Inf,fill=ru,alpha=.05)+
    ggplot2::geom_hline(yintercept=1,linetype=2,colour="#87939c")+
    ggplot2::geom_line(colour=te,linewidth=.8)+ggplot2::scale_x_log10()+
    ggplot2::annotate("text",x=.034,y=.37,label="alpha*k~sensitive",parse=TRUE,size=2.7)+
    ggplot2::annotate("text",x=6,y=.44,label="Near\nsaturation",colour=ru,size=2.7)+
    ggplot2::labs(title="b  Duration changes sensitivity",x=expression(Dimensionless~duration~k*tau),
      y=expression(Rate~elasticity~partialdiff*ln*C[acc]/partialdiff*ln*k))+theme_ems(base_size)
  tie <- fe_inventory*.25/2.25; fes <- seq(0,tie*1.8,length.out=150)
  ledger <- data.frame(fes=fes,fe=fe_inventory*.25,sulfide=2.25*fes)
  c <- ggplot2::ggplot(ledger,ggplot2::aes(x=.data$fes))+
    ggplot2::geom_line(ggplot2::aes(y=.data$fe,colour="Fe(II)"),linewidth=.8)+
    ggplot2::geom_line(ggplot2::aes(y=.data$sulfide,colour="FeS"),linewidth=.8)+
    ggplot2::geom_vline(xintercept=tie,linetype=3,colour="#87939c")+
    ggplot2::annotate("point",x=tie,y=fe_inventory*.25,colour=n,size=2.3)+
    ggplot2::annotate("text",x=tie*1.38,y=fe_inventory*.10,label=sprintf("Equal demand\nat %.2f",tie),size=2.7)+
    ggplot2::scale_colour_manual(values=c("Fe(II)"=n,"FeS"=ru),name=NULL)+
    ggplot2::labs(title="c  Stoichiometric ranking",x=expression(FeS~inventory~(mmol~kg^{-1})),y=expression(O[2]~demand~(mmol~kg^{-1})))+
    theme_ems(base_size)+
    (if(utils::packageVersion("ggplot2") >= "3.5.0")
      ggplot2::theme(legend.position="inside",legend.position.inside=c(.30,.89))
      else ggplot2::theme(legend.position=c(.30,.89)))+
    ggplot2::theme(legend.direction="vertical",legend.background=ggplot2::element_rect(fill="white",colour=NA))
  sensitivity <- 1/abs(elast(k*dp[,1])-elast(k*dp[,2]))
  labs <- paste(dp[,1],dp[,2],sep=", ")
  dd <- data.frame(pair=factor(seq_len(nrow(dp)),levels=rev(seq_len(nrow(dp))),labels=rev(labs)),
    sensitivity=sensitivity, saturated=k*dp[,1]>=1.5)
  dd$label <- ifelse(is.finite(sensitivity),sprintf("%.2f",sensitivity),"Inf")
  # Infinite conditioning is displayed explicitly, never passed as a finite bar.
  dd$display <- ifelse(is.finite(sensitivity),sensitivity,0)
  d <- ggplot2::ggplot(dd,ggplot2::aes(x=.data$display,y=.data$pair))+
    ggplot2::geom_col(ggplot2::aes(fill=.data$saturated),width=.60)+
    ggplot2::geom_text(ggplot2::aes(label=.data$label),hjust=-.15,size=2.7)+
    ggplot2::scale_fill_manual(values=c("FALSE"=te,"TRUE"=go),guide="none")+
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,.20)))+
    ggplot2::labs(title="d  Two-window conditioning",x="Local relative inverse sensitivity",y=expression(Durations~(tau[1]*","~tau[2])))+
    theme_ems(base_size)+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_line(colour="#edf0f2"))
  .hrri_assemble(list(identified_set=a,elasticity=b,oxygen=c,conditioning=d),ncol=2)
}
