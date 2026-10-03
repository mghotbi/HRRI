.hrri_accuracy_paper <- function(acc,panels,score_label,target_label,point_alpha,
                                show_clusters,base_size,ncol,cluster_label) {
  n <- "#274b6a"; te <- "#167c80"; ru <- "#c6633b"; go <- "#b58130"; pu <- "#807195"
  d <- acc$data; d$cluster <- factor(d$cluster)
  cm <- stats::aggregate(cbind(score,target)~cluster,d,mean)
  ag <- stats::setNames(acc$agreement$row_level,acc$agreement$statistic)
  size <- base_size/3.4
  theme <- theme_ems(base_size)+ggplot2::theme(legend.position="top",
    legend.title=ggplot2::element_blank(),legend.text=ggplot2::element_text(size=base_size-2),
    plot.title=ggplot2::element_text(size=base_size+1))
  pts <- if(isTRUE(show_clusters)) ggplot2::geom_point(ggplot2::aes(colour=.data$cluster),size=.55,alpha=point_alpha,show.legend=FALSE) else
    ggplot2::geom_point(colour=te,size=.55,alpha=point_alpha,stroke=0)
  annotation <- function(label) ggplot2::annotate("text",x=-Inf,y=Inf,label=label,hjust=-.03,vjust=1.1,size=size,colour="#233441",lineheight=1.1)
  out <- list()
  if("calibration" %in% panels) {
    co <- if(length(unique(d$target))>1L) stats::coef(stats::lm(score~target,d)) else c(NA_real_,NA_real_)
    p <- ggplot2::ggplot(d,ggplot2::aes(x=.data$target,y=.data$score))+pts+
      ggplot2::geom_abline(slope=1,intercept=0,colour="#73818c",linetype=2,linewidth=.45)
    if(all(is.finite(co))) p <- p+ggplot2::geom_abline(slope=co[2],intercept=co[1],colour=ru,linewidth=.6)
    p <- p+ggplot2::geom_point(data=cm,shape=21,fill="white",colour=n,size=2.2,stroke=.65)+
      annotation(sprintf("r = %.3f   CCC = %.3f\nRMSE = %.3f",ag[["pearson_r"]],ag[["lins_ccc"]],ag[["rmse"]]))+
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.04,.15)))+
      ggplot2::labs(title="a  Association versus agreement",x=target_label,y=score_label)+theme
    out$calibration <- p
  }
  if("agreement" %in% panels) {
    d$avg <- (d$score+d$target)/2; d$diff <- d$score-d$target
    cm$avg <- (cm$score+cm$target)/2; cm$diff <- cm$score-cm$target
    bias <- mean(cm$diff); sd <- stats::sd(cm$diff); loa <- bias+c(-1,1)*1.96*sd
    p <- ggplot2::ggplot(d,ggplot2::aes(x=.data$avg,y=.data$diff))+pts+
      ggplot2::geom_hline(yintercept=bias,colour=ru,linewidth=.5)
    if(is.finite(sd)) p <- p+ggplot2::geom_hline(yintercept=loa,colour=ru,linetype=2,linewidth=.5)
    label <- if(is.finite(sd)) sprintf("%s-mean bias = %.3f\nMean limits [%.3f, %.3f]",cluster_label,bias,loa[1],loa[2]) else "Too few clusters for mean limits"
    p <- p+ggplot2::geom_point(data=cm,shape=21,fill="white",colour=n,size=2.2,stroke=.65)+
      annotation(label)+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.04,.20)))+
      ggplot2::labs(title="b  Differences at two levels",x="Mean of score and target",y="Score minus target")+theme
    out$agreement <- p
  }
  if("precision" %in% panels) {
    dr <- acc$draws
    cl <- if(!is.null(dr$cluster)) dr$cluster[,"r"] else numeric()
    rw <- if(!is.null(dr$naive_r)) dr$naive_r else numeric()
    cl <- cl[is.finite(cl)]; rw <- rw[is.finite(rw)]
    if(!length(cl)||!length(rw)) {
      p <- ggplot2::ggplot()+ggplot2::annotate("text",x=0,y=0,label="Bootstrap draws unavailable\nRun rri_accuracy() with n_boot > 0",size=size)+
        ggplot2::labs(title="c  Conditional bootstrap precision",x=NULL,y=NULL)+theme+
        ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
    } else {
      gr <- range(c(cl,rw)); pad <- max(diff(gr)*.07,.005); gr <- gr+c(-pad,pad)
      pair <- list(cl,rw); src <- c(paste("Whole",tolower(cluster_label)),"Rows (naive)")
      den <- vector("list",2); spikes <- list()
      for(i in 1:2) {
        z <- pair[[i]]
        if(length(z)>1L && stats::sd(z)>0) {
          fit <- stats::density(z,bw=stats::sd(z)*length(z)^(-1/5),from=gr[1],to=gr[2],n=512)
          den[[i]] <- data.frame(r=fit$x,density=fit$y,source=src[i])
        } else spikes[[length(spikes)+1]] <- data.frame(r=z[1],source=src[i])
      }
      dd <- do.call(rbind,den); if(is.null(dd)) dd <- data.frame(r=numeric(),density=numeric(),source=character())
      pal <- stats::setNames(c(te,pu),src)
      a <- (1-acc$conf)/2; widths <- vapply(pair,function(z) diff(stats::quantile(z,c(a,1-a))),numeric(1))
      p <- ggplot2::ggplot(dd,ggplot2::aes(x=.data$r,y=.data$density,colour=.data$source,fill=.data$source))+
        ggplot2::geom_area(alpha=.12,position="identity",colour=NA)+ggplot2::geom_line(linewidth=.65)+
        ggplot2::scale_colour_manual(values=pal,limits=src)+ggplot2::scale_fill_manual(values=pal,limits=src)+
        ggplot2::geom_vline(xintercept=ag[["pearson_r"]],colour=n,linetype=3,linewidth=.45)+
        annotation(sprintf("%.0f%% widths\n%s %.3f\nRows %.3f",100*acc$conf,cluster_label,widths[1],widths[2]))+
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.01,.15)))+
        ggplot2::labs(title="c  Resampling unit and precision",x="Bootstrap Pearson r",y="Density")+theme
      if(length(spikes)) p <- p+ggplot2::geom_point(data=do.call(rbind,spikes),ggplot2::aes(x=.data$r,y=0),inherit.aes=FALSE,colour=te)
    }
    out$precision <- p
  }
  if("error" %in% panels) {
    dc <- acc$decomposition
    lab <- c(squared_bias="Squared\nmean bias",variance_mismatch="Spread\nmismatch",lack_of_correlation="Lack of\ncorrelation")
    dc$label <- factor(unname(lab[dc$component]),levels=rev(unname(lab)))
    zero <- isTRUE(all.equal(attr(dc,"mse"),0,tolerance=0))
    if(zero) dc$percent <- 0
    pal <- stats::setNames(c(n,go,te),unname(lab))
    p <- ggplot2::ggplot(dc,ggplot2::aes(x=.data$percent,y=.data$label,fill=.data$label))+
      ggplot2::geom_col(width=.55)+ggplot2::scale_fill_manual(values=pal,guide="none")+
      ggplot2::geom_text(ggplot2::aes(label=if(zero) "0; share undefined" else sprintf("%.1f%%",.data$percent)),hjust=-.12,size=size)+
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,if(zero) .5 else .25)))+
      ggplot2::labs(title="d  Sources of disagreement",x=if(zero) "Zero MSE; shares undefined" else "Share of mean squared error (%)",y=NULL)+theme+
      ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),panel.grid.major.x=ggplot2::element_line(colour="#edf0f2"))
    out$error <- p
  }
  .hrri_assemble(out,ncol=ncol)
}
