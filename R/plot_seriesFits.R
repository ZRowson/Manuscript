#' Plot Time-Series Data for a Test Chemical and Active Curve-Fits
#'
#' @author Zachary Rowson \email{Rowson.Zachary@@epa.gov}
#'
# @description
# Plots total movement per 2-minute bin time-series data of a test chemical and
# it's corresponding vehicle control using a Sample Averaged (SA) perspective
# and combines this plot with active curve-fit plots.
#
# @details
# Created: 07/08/2022
# Last edit: 07/08/2022
#
# Sample Averaged (SA) perspective fits to the mean speed at each
# time period for a concentration group.


library(data.table)
library(gabi)

# Functions ---------------------------------------------------------------


seriesFits <- function(pmr0, chemical, tcplOut, list) {

  # Plot SA series
  SAplot <- gabi::plot_tSeries(pmr0, chemical = chemical)

  # Find active endpoints for chemical and curve-fit with plots
  activeEndps <- tcplOut[name==chemical & hitcall>0.8, endp]
  activeRows <- list[activeEndps]
  activeFits <- lapply(activeRows, gabi::concRespCoreZR, do.plot = TRUE, verbose.plot = FALSE)

  # Isolate curve-fit plots, edit, and combine with SA plot and return
  curvePlots <- lapply(activeFits, `[[`, "plot")
  legend <- cowplot::get_legend(curvePlots[[1]])
  curvePlots1 <- lapply(curvePlots, edit_plot, chemical=chemical)

  return( list("SAplot" = SAplot,
              "curves" = curvePlots1,
              "legend" = legend) )
  }

edit_plot <- function(plot, chemical) {
  title <- plot$labels$title
  plot$labels$title <- substr(title, nchar(chemical)+6, nchar(title))
  plot$labels$subtitle <- NULL
  plot <- plot + theme(legend.position="none")
  return(plot)
}


# Apply to Active Chemicals -----------------------------------------------


# Remember to set working directory

tcplOut <- as.data.table(DNT60_tcpl_out)[!endp%in%c("AUC_L","AUC_D","AUC_T")]
activeChems <- tcplOut[hitcall > 0.8, unique(name)]

  chemical <- activeChems[[22]]
  rowList <- rows_n[[chemical]]
  plotList <- seriesFits(pmr0=pmr0, chemical=chemical, tcplOut=tcplOut, list=rowList)

  p1 <- cowplot::plot_grid(plotlist=plotList[["curves"]], nrow=2)
  p2 <- cowplot::plot_grid(p1, plotList[["legend"]], rel_widths = c(3, .4))

  fileName = paste0(chemical, "_SAplotCurves", ".png")
  png(filename = fileName, width = 40, height = 40, unit = "cm", res=300)
  cowplot::plot_grid(plotList[["SAplot"]], p2, nrow=2)
  dev.off()

