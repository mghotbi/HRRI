#' @title Profile of Available HRRI Diagnostics
#'
#' @description
#' Displays available diagnostic summaries labelled Capacity, Connectivity,
#' Kinetics and Memory. The composite RRI is not an axis: it is built from
#' the plant, soil and microbial domains rather than from these four
#' properties, so averaging across it is not defensible. These axes are
#' operational descriptors returned by \code{\link{rri_property_scores}};
#' they are not direct measurements or identified estimates of the theoretical
#' mechanisms bearing the same names.
#'
#' @param props A list returned by \code{\link{rri_property_scores}}, or a
#'   named numeric vector with elements \code{Capacity}, \code{Connectivity},
#'   \code{Kinetics}, \code{Memory} (values in `[0, 1]`).
#' @param rri_value Optional numeric. Composite RRI, reported in the subtitle
#'   for reference. It is not plotted as an axis and does not enter
#'   any average across property axes; no such average is calculated.
#'   Defaults to \code{props$rri_summary} if available.
#' @param group_list Optional named list of property score vectors, one per
#'   group (e.g., per thaw stage or treatment). If supplied, multiple
#'   groups are displayed in separate panels.
#' @param rec Optional recovery table for data-derived coverage notes; use only
#'   the same recovery table used to calculate props. With multiple groups supply
#'   explicit notes instead.
#' @param notes Optional named character vector (one note per property), or a
#'   named list of such vectors matching group_list. Describes score provenance.
#' @param type Display type: separate-axis dot profile (default), or an explicit
#'   radar display for compatibility. Radar area is not a quantitative summary.
#' @param fill_alpha Numeric in `[0, 1]`. Polygon fill transparency.
#' @param colours Character vector of polygon outline/fill colours, recycled
#'   across groups.
#' @param show_values Logical. Annotate each axis tip with the numeric score.
#' @param title Character. Plot title.
#' @param base_size Numeric. Base font size.
#'
#' @details
#' The default profile displays each descriptor separately, without an overall
#' mean or a common favourable direction. Missing values are labelled explicitly.
#' Out-of-range finite scores and infinite values are rejected, not clipped.
#' The optional radar chart uses Cartesian coordinates constructed with \code{ggplot2}; no
#' external radar-chart package is required. Each available axis runs from 0
#' (centre) to 1 (rim). Polygon area has no quantitative meaning, and axes based
#' on different transformations are not necessarily commensurable.
#'
#' \strong{Axis meanings:}
#' \itemize{
#'   \item \strong{Capacity} --- oxidative-oriented soil feature composite; not
#'     accessible capacity unless the caller calculates and supplies it explicitly.
#'   \item \strong{Connectivity} --- association or network-topology descriptor;
#'     not demonstrated electron transfer.
#'   \item \strong{Kinetics} --- cohort- or timescale-relative recovery-speed descriptor;
#'     not a mineral exchange rate.
#'   \item \strong{Memory} --- loop-area and persistent-displacement descriptor;
#'     not an identified causal memory state.
#'   
#' }
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom rlang .data
#'
#' @examples
#' sim <- simulate_redox_holobiont(
#'   n_plot = 3, n_depth = 2, n_plant = 3, n_time = 14,
#'   p_micro = 30, seed = 99
#' )
#'
#' res <- rri_pipeline_st(
#'   ROS_flux     = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data   = sim$micro_data,
#'   id           = sim$id,
#'   reducer      = "per_domain",
#'   scaling      = "pnorm"
#' )
#'
#' rec <- rri_recovery_metrics(
#'   res           = res,
#'   id            = sim$id,
#'   time_col      = "time",
#'   group_cols    = c("plot", "depth", "plant_id"),
#'   perturb_start = 5,
#'   perturb_end   = 8
#' )
#'
#' props <- rri_property_scores(
#'   res       = res,
#'   rec       = rec,
#'   soil_df   = sim$Eh_stability,
#'   eac_col   = "EAC",
#'   edc_col   = "EDC",
#'   humic_col = "dissolved_organic_matter_redox"
#' )
#'
#' plot_rri_properties(props)
#'
#' @export
plot_rri_properties <- function(
  props, rri_value = NULL, group_list = NULL, fill_alpha = 0.20,
  colours = c("#1A3A5C", "#E07B39", "#2E7D32", "#7B3294", "#B2182B"),
  show_values = TRUE, title = "Operational descriptors and their input coverage", base_size = 10,
  type = c("profile", "radar"), rec = NULL, notes = NULL
) {
  type <- match.arg(type)
  axes <- c("Capacity", "Connectivity", "Kinetics", "Memory")
  extract <- function(x) {
    if (is.list(x) && "property_scores" %in% names(x)) x <- x$property_scores
    if (!is.numeric(x) || is.null(names(x)) || anyDuplicated(names(x)))
      stop("Property scores must be a uniquely named numeric vector.", call. = FALSE)
    z <- x[axes]
    if (any(is.infinite(z)) || any(z < 0 | z > 1, na.rm = TRUE))
      stop("Property scores must be in [0, 1] or missing; values are not clipped.", call. = FALSE)
    as.numeric(z)
  }
  if (!length(colours) || anyNA(colours)) stop("Supply at least one colour.", call. = FALSE)
  if (length(fill_alpha) != 1L || !is.finite(fill_alpha) || fill_alpha < 0 || fill_alpha > 1)
    stop("fill_alpha must be in [0, 1].", call. = FALSE)
  if (is.null(group_list)) {
    scores <- list(system = props)
    if (is.null(rri_value) && is.list(props)) rri_value <- props$rri_summary
  } else {
    if (!is.list(group_list) || !length(group_list) || is.null(names(group_list)) ||
        anyNA(names(group_list)) || any(!nzchar(names(group_list))) || anyDuplicated(names(group_list)))
      stop("group_list must be a nonempty list with unique, nonempty names.", call. = FALSE)
    scores <- group_list
  }
  if (!is.null(rri_value) && (!is.numeric(rri_value) || length(rri_value) != 1L ||
      is.infinite(rri_value))) stop("rri_value must be a numeric scalar or NULL.", call. = FALSE)
  d <- do.call(rbind, lapply(seq_along(scores), function(i) {
    data.frame(group = names(scores)[i], axis = axes, value = extract(scores[[i]]))
  }))
  d$group <- factor(d$group, levels = names(scores))
  pal <- stats::setNames(rep(colours, length.out = length(scores)), names(scores))
  subtitle <- "Operational descriptors; no common favourable direction or overall mean"
  if (!is.null(rri_value) && is.finite(rri_value))
    subtitle <- paste0(subtitle, sprintf("\nDomain composite RRI = %.3f (separate quantity)", rri_value))
  caption <- "Missing = unavailable from supplied inputs. Scores depend on their component definitions."
  if (type == "profile") {
    d$axis <- factor(d$axis, levels = rev(axes))
    d$y <- as.numeric(d$axis)
    d$note <- "User-supplied descriptor; interpretation follows its definition"
    for (i in seq_along(scores)) {
      ix <- d$group == names(scores)[i]
      obj <- scores[[i]]
      nt <- stats::setNames(rep("User-supplied descriptor; interpretation follows its definition",4),axes)
      if(is.list(obj) && !is.null(obj$property_table)) {
        tab <- obj$property_table
        nt <- stats::setNames(as.character(tab$method[match(axes,tab$property)]),axes)
        nt["Capacity"] <- "Oxidative-oriented feature score; not electron capacity"
        if(grepl("cross_domain",nt["Connectivity"])) nt["Connectivity"] <- "Cross-domain association; not physical accessibility"
      }
      if(!is.null(rec) && length(scores)==1L) {
        count <- function(n) if(n %in% names(rec)) sum(is.finite(rec[[n]])) else NA_integer_
        nt["Kinetics"] <- sprintf("%s lag values; %s finite rates; inspect component weights",count("tau_lag"),count("k_recovery"))
        nt["Memory"] <- sprintf("%s displacements; %s finite forcing-response loop areas",count("incomplete_return_frac"),count("H_hysteresis"))
      }
      supplied <- if(is.list(notes)) notes[[names(scores)[i]]] else notes
      if(!is.null(supplied)) {
        if(!is.character(supplied)||is.null(names(supplied))||!all(axes %in% names(supplied)))
          stop("notes must name all four properties.",call.=FALSE)
        nt <- supplied[axes]
      }
      d$note[ix] <- unname(nt[as.character(d$axis[ix])])
    }
    d$note[is.na(d$note)] <- "Provenance not supplied"
    d$note <- vapply(d$note,function(s) paste(strwrap(s,width=70),collapse="\n"),character(1))
    d$co <- c(Capacity="#274b6a",Connectivity="#167c80",Kinetics="#b58130",Memory="#c6633b")[as.character(d$axis)]
    if(!missing(colours)) d$co <- rep(colours,length.out=nrow(d))
    ok <- is.finite(d$value)
    z <- d[ok, , drop = FALSE]
    ticks <- d[rep(seq_len(nrow(d)),each=3),,drop=FALSE]
    ticks$x <- rep(c(0,.5,1),nrow(d))
    p <- ggplot2::ggplot(d,ggplot2::aes(x=.data$value,y=.data$y))+
      ggplot2::geom_segment(ggplot2::aes(x=0,xend=1,yend=.data$y),colour="#d9e0e5",linewidth=1.4)+
      ggplot2::geom_point(data=z,ggplot2::aes(colour=.data$co),size=3.2)+
      ggplot2::geom_text(data=d[!ok,,drop=FALSE],ggplot2::aes(x=.5,label="Missing"),colour="#657680",size=base_size/3.3)+
      ggplot2::geom_text(data=ticks,ggplot2::aes(x=.data$x,y=.data$y-.19,label=.data$x),size=base_size/3.6,colour="#657680")+
      ggplot2::geom_text(ggplot2::aes(x=0,y=.data$y-.43,label=.data$note),hjust=0,vjust=.5,size=base_size/3.8,lineheight=1.05,colour="#233441")+
      ggplot2::scale_colour_identity()+
      ggplot2::scale_y_continuous(breaks=1:4,labels=rev(axes),limits=c(.35,4.4),expand=c(0,0))+
      ggplot2::scale_x_continuous(limits=c(-.02,1.04),expand=c(0,0))+
      ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme_void(base_size=base_size)+
      ggplot2::theme(axis.text.y=ggplot2::element_text(face="bold",colour="#233441",size=base_size,
        margin=ggplot2::margin(r=10)))
    if(isTRUE(show_values)) p <- p+ggplot2::geom_text(data=z,
      ggplot2::aes(y=.data$y+.17,label=sprintf("%.3f",.data$value),colour=.data$co),
      size=base_size/3.3,fontface="bold")
  } else {
    d$angle <- rep(seq(0, 2 * pi, length.out = 5)[1:4], length(scores))
    d$x <- d$value * cos(d$angle); d$y <- d$value * sin(d$angle)
    rings <- do.call(rbind, lapply(seq(0.25, 1, 0.25), function(r) {
      a <- seq(0, 2*pi, length.out = 101)
      data.frame(x = r*cos(a), y = r*sin(a), ring = r)
    }))
    labs <- data.frame(axis = axes, x = 1.22*cos(d$angle[1:4]), y = 1.22*sin(d$angle[1:4]))
    complete <- vapply(split(d$value, d$group), function(v) all(is.finite(v)), logical(1))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_path(data = rings, ggplot2::aes(group = .data$ring), colour = "grey85") +
      ggplot2::geom_polygon(data = d[d$group %in% names(complete)[complete], , drop = FALSE],
        ggplot2::aes(group = .data$group, colour = .data$group, fill = .data$group), alpha = fill_alpha) +
      ggplot2::geom_point(data = d[is.finite(d$value), , drop = FALSE],
                           ggplot2::aes(colour = .data$group), size = 2.8) +
      ggplot2::geom_text(data = labs, ggplot2::aes(label = .data$axis), size = base_size / 4) +
      ggplot2::scale_colour_manual(values = pal, guide = "none") +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.5, 1.5)) +
      ggplot2::theme_void(base_size = base_size)
    dl <- d
    dl$x <- 1.02*cos(dl$angle); dl$y <- 1.02*sin(dl$angle)
    dl$label <- ifelse(is.finite(dl$value), sprintf("%.2f", dl$value), "Missing")
    if (!isTRUE(show_values)) dl <- dl[!is.finite(dl$value), , drop = FALSE]
    p <- p + ggplot2::geom_label(data = dl, ggplot2::aes(label = .data$label), size = base_size / 4.5)
    caption <- paste(caption, "Radar area has no quantitative meaning.")
  }
  if (length(scores) > 1L) p <- p + ggplot2::facet_wrap(~group, ncol = 2)
  p + ggplot2::labs(title = title, subtitle = if(type == "radar") subtitle else NULL,
      caption = if(type == "radar") caption else "Separate 0-1 ranges; no common favourable direction or overall mean.") +
    ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA),plot.title = ggplot2::element_text(face = "bold",colour="#233441"),plot.title.position="plot",
      plot.subtitle = ggplot2::element_text(size = base_size * 0.75),
      plot.caption = ggplot2::element_text(size = base_size * 0.65, hjust = 0),
      plot.margin = ggplot2::margin(10, 12, 10, 12))
}
