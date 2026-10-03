
if (!exists("figure_dir", inherits=FALSE)) figure_dir <- file.path(tempdir(),"HRRI_paper_figures")
source(system.file("examples","paper_figure_objects.R",package="HRRI"),local=TRUE)
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
sizes <- rbind(c(7.4,4.8),c(7.4,6.4),c(7.4,6),c(8,4.6),c(7.4,4.5),c(7.4,6.5))
for(i in seq_along(paper_figures)) {
  p <- paper_figures[[i]]
  if(is.list(p) && !inherits(p,c("ggplot","patchwork")))
    stop("Install suggested package patchwork to assemble all paper panels.")
  nm <- names(paper_figures)[i]
  ggplot2::ggsave(file.path(figure_dir,paste0(nm,".pdf")),p,width=sizes[i,1],height=sizes[i,2],device="pdf")
  ggplot2::ggsave(file.path(figure_dir,paste0(nm,".png")),p,width=sizes[i,1],height=sizes[i,2],dpi=600)
  if(requireNamespace("svglite",quietly=TRUE))
    ggplot2::ggsave(file.path(figure_dir,paste0(nm,".svg")),p,width=sizes[i,1],height=sizes[i,2],device=svglite::svglite)
}
if(!requireNamespace("svglite",quietly=TRUE)) message("Install suggested svglite for editable-text SVG export; PDF and PNG were saved.")
utils::write.csv(gal_rec,file.path(figure_dir,"recovery_values_and_fit_status.csv"),row.names=FALSE)
utils::write.csv(props$property_table,file.path(figure_dir,"property_values_and_methods.csv"),row.names=FALSE)
utils::write.csv(acc$agreement,file.path(figure_dir,"agreement_statistics.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(figure_dir,"sessionInfo.txt"))
message("Six paper figures saved to ",normalizePath(figure_dir))
