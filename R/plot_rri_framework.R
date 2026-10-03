#' Complementary mechanistic and observational routes
#'
#' A vector schematic of the two HRRI routes. The response curve and biological
#' symbols are illustrative, not measured or simulated observations. Electron
#' accepting and donating capacities must remain directionally separate.
#' @param base_size Base font size in points.
#' @return A ggplot object; export at approximately 7.4 by 4.8 inches.
#' @examples
#' plot_rri_framework()
#' @export
plot_rri_framework <- function(base_size = 10) {
  navy <- "#274b6a"; rust <- "#c6633b"; teal <- "#167c80"; ink <- "#233441"
  p <- ggplot2::ggplot()
  txt <- function(x,y,label,size=8,colour=ink,face="plain",parse=FALSE,hjust=.5)
    ggplot2::annotate("text",x=x,y=y,label=label,size=size*base_size/10/2.845,
      colour=colour,fontface=face,parse=parse,hjust=hjust,lineheight=1.15)
  line <- function(x,y,xend,yend,colour=navy,arrow=FALSE,linetype=1) {
    args <- list(geom="segment",x=x,y=y,xend=xend,yend=yend,colour=colour,
      linewidth=.35,linetype=linetype)
    if(arrow) args$arrow <- grid::arrow(length=grid::unit(1.7,"mm"),type="closed")
    do.call(ggplot2::annotate,args)
  }
  ellipse <- function(x,y,rx,ry,fill,colour=NA,angle=0) {
    t <- seq(0,2*pi,length.out=81); a <- angle*pi/180
    d <- data.frame(x=x+rx*cos(t)*cos(a)-ry*sin(t)*sin(a),
                    y=y+rx*cos(t)*sin(a)+ry*sin(t)*cos(a))
    ggplot2::geom_polygon(data=d,ggplot2::aes(x=.data$x,y=.data$y),
                           fill=fill,colour=colour,linewidth=.4,inherit.aes=FALSE)
  }
  for(i in 1:2) {
    x <- c(.1,5.9)[i]; co <- c(navy,rust)[i]
    p <- p+txt(x,6.15,c("a  MECHANISTIC","b  OBSERVATIONAL")[i],11,co,"bold",hjust=0)+
      txt(x,5.86,c("Declared inventories and exchange parameters","Aligned measurements on the same units")[i],7.6,hjust=0)+
      line(x,5.68,x+4,5.68,co)
  }
  for(i in 1:3) {
    x <- c(.8,2.1,3.4)[i]
    p <- p+ellipse(x,5.07,.39,.39,"#edf4f7",navy)+
      txt(x,5.07,c("Q[i]","alpha[i]","k[i]*','~tau")[i],13,navy,parse=TRUE)+
      txt(x,4.48,c("Inventory","Accessibility","Rate / time")[i],8)+line(x,4.32,x,4.15)
  }
  p <- p+line(.8,4.15,3.4,4.15)+line(2.1,4.15,2.1,3.95,arrow=TRUE)+
    txt(2.1,3.62,"C[acc](tau)==sum(Q[i]*alpha[i]*(1-e^{-k[i]*tau}),i)",11,navy,parse=TRUE)+
    txt(2.1,3.15,"Event-window accessible electron capacity",7.7,face="bold")+
    line(.35,2.99,3.85,2.99,"#acb9c3")+
    txt(2.1,2.65,"D[O[2]]==sum(nu[j]*n[j*','*red],j)",11,navy,parse=TRUE)+
    txt(2.1,2.16,"Separate stoichiometric oxygen demand",8)
  for(i in 1:3) {
    x <- c(6.35,7.65,8.95)[i]
    p <- p+ellipse(x,5.07,.39,.39,"#fff9f3",rust)+
      txt(x,4.48,c("SOIL","PLANT","MICROBIAL")[i],7.7,rust,"bold")+
      line(x,4.30,7.65,4.03,rust)
  }
  for(j in 0:2) p <- p+ellipse(6.35,5.22-j*.14,.25,.075,c("#a6754f","#bb9975","#d8c5a7")[j+1],"white")
  p <- p+line(7.65,4.85,7.65,5.30,"#638d42")+
    ellipse(7.54,5.13,.15,.065,"#82a959",angle=-35)+
    ellipse(7.77,5.23,.15,.065,"#82a959",angle=35)+
    ggplot2::annotate("rect",xmin=7.52,xmax=7.78,ymin=4.79,ymax=4.93,fill="#b38a68")
  for(i in 1:3) p <- p+ellipse(c(8.78,9.05,8.94)[i],c(5.15,5.19,4.94)[i],.14,.05,"white",teal,c(-35,35,55)[i])
  p <- p+txt(7.65,3.79,"ALIGN  -  ORIENT  -  SCALE",9,rust,"bold")+
    txt(7.65,3.46,"Domain coverage and missingness retained",7.7)+
    txt(5,4.11,"NO INVERSION",7.5,face="bold")+
    txt(5,3.68,"Scores do not\nidentify Q, alpha or k",7.3)+
    line(5,2.22,5,3.10,"#aab3b9",linetype=3)
  t <- seq(5.99,9.62,length.out=300)
  d <- data.frame(x=t,y=2.75-.63*exp(-((t-7.05)/.43)^2)+.15*exp(-((t-8.34)/.45)^2))
  p <- p+ggplot2::annotate("rect",xmin=6.7,xmax=7.4,ymin=2.05,ymax=3.05,fill="#b58130",alpha=.12)+
    ggplot2::geom_line(data=d,ggplot2::aes(x=.data$x,y=.data$y),colour=rust,linewidth=.7)+
    line(5.99,2.75,9.62,2.75,"#97a2ab",linetype=2)+
    txt(7.05,3.13,"event",7.7)+txt(7.08,1.90,"decline",7.5)+txt(8.8,2.36,"return / overshoot",7.5)+
    txt(5,1.58,"JOINT INTERPRETATION",9,teal,"bold")+
    ellipse(5,.77,.50,.50,"white",teal)+txt(5,.89,"HRRI",10,teal,"bold")+
    txt(5,.56,"traceable\nworkflow",6.5)+
    ggplot2::annotate("curve",x=2.1,y=1.91,xend=4.48,yend=.90,curvature=.12,
      colour=navy,linewidth=.4,arrow=grid::arrow(length=grid::unit(1.7,"mm"),type="closed"))+
    ggplot2::annotate("curve",x=7.65,y=1.91,xend=5.52,yend=.90,curvature=-.12,
      colour=rust,linewidth=.4,arrow=grid::arrow(length=grid::unit(1.7,"mm"),type="closed"))+
    txt(2.03,.60,"Mechanistic quantities",8.7,navy)+txt(7.93,.60,"Observational diagnostics",8.7,rust)
  p+ggplot2::coord_fixed(xlim=c(-.25,10),ylim=c(.12,6.48),clip="off")+
    ggplot2::theme_void(base_size=base_size)+ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA),plot.margin=ggplot2::margin(8,8,8,8))
}
