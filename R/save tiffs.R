## ---------------------------
##
## Script Name: Save Manuscript Figures as tiffs
##
## Purpose of Script: Save Figures associated with gabi Manuscript as high-resolution tiffs
##
## Author: Zachary Rowson
##
## Date Created: 2022-11-23
## Last Edit: 2023-01-09
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/StephanieProjects/gabi/Manuscript"
##
## ---------------------------
##
## Notes: Runnning this file will require running "R Analyses and Figures for Paper.Rmd"
##        to fill environment with figure objects.
##
##
## ---------------------------


setwd("~/Desktop/StephanieProjects/gabi/Manuscript/Figures/Figures v8")
# Make sure to set proper working directory


# Main Figures ------------------------------------------------------------


# Figure 2: Typical Raw Movement Data and Effect of Transformation procedure and Fluoxetine Exposure on Distributions of Endpoint Data
tiff(filename = "Fluoxetine data and avgS endp.tiff", height = 10, width = 21.5, units = "in", res=300)
FlxPlots
dev.off()

# Figure 3: Count of Active Chemicals by Phase and Endpoint Category
tiff(filename = "Activity by cat and phase.tiff", height = 10, width = 21.5, units = "in", res=300)
cowplot:::plot_grid(activesCatCntPlot, activesPhsCntPlot,
                    labels = "AUTO", label_size = 24,
                    ncol = 2)
dev.off()

# Figure 4: Raw Speed Data and Active Curve-Fits and for Chlorpyrifos
tiff(filename = "Chlorpyrifos Raw Data Plots and Active Endpoints.tiff", height = 10, width = 15, units = "in", res=300)
cowplot::plot_grid(p1_CPF, p2_CPF, nrow=2)
dev.off()

# Figure 5: BMC Confidence Interval Overlap
tiff(filename = "BMC CI overlap.tiff", height = 9, width = 12, units = "in", res=300)
print(BMC_facet)
dev.off()

# Figure 6: BMC Heatmap
tiff(filename = "BMC Heatmap.tiff", height = 6, width = 7, units = "in", res=300)
draw(htlist, merge_legend = FALSE, annotation_legend_list = lgd_list, annotation_legend_side = "top", align_annotation_legend = "heatmap_center")
dev.off()

# Supplemental Figure 7: Developmentally and Neurologically Toxic LOELs
tiff(filename = "LOELs plot.tiff", height = 9, width = 18, units = "in", res=300)
print(compare_LOELs)
dev.off()


# Supplemental Figures ----------------------------------------------------


# Supplemental Figure 2: Fluoxetine Vehicle Control and 4 uM Treatment Group Transformed Endpoint Data
tiff(filename = "endp_raw_facet.tiff", height = 10, width = 15, units = "in", res=300)
print(raw_plot_facet)
dev.off()

# Supplemental Figure 3: Fluoxetine Vehicle Control and 4 uM Treatment Group Transformed Endpoint Data
tiff(filename = "endp_n_facet.tiff", height = 10, width = 15, units = "in", res=300)
print(trans_plot_facet)
dev.off()

# Supplemental Figure 4: Combinations of Endpoint Activity by Endpoint Category or Experimental Phase
tiff(filename = "Combinations of Activity.tiff", height = 9, width = 14, units = "in", res = 300)
cowplot::plot_grid(plot, plot1, rel_widths = c(9,8),
                   labels = "AUTO", label_size = 24)
dev.off()

# Supplemental Figure 5: Activity in Endpoints Other than Average Speed
tiff(filename = "Areas of Additional Activity.tiff", height = 5, width = 7, units = "in", res = 300)
draw(htlist2, merge_legend = FALSE, annotation_legend_list = lgd_list,
     annotation_legend_side = "top", align_annotation_legend = "heatmap_center")

decorate_heatmap_body(paste0("BMC log(","\U03BC","M)"), row_slice = 1, column_slice = 1, {
  grid.rect(unit(1.025,'npc'), unit(0.7,'npc'),
            width = (length(co[[1]]) + length(co[[2]]))/length(co[[1]]) * unit(0.415,'npc') + unit(1,'mm'),
            height = (length(ro[[1]]) + length(ro[[2]]))/length(ro[[1]]) * unit(0.344,'npc') + unit(1,'mm'),
            gp = gpar(lwd=2, lty=1, fill=NA, col="red"), just=c('left','top')
  )
})
decorate_heatmap_body(paste0("BMC log(","\U03BC","M)"), row_slice = 1, column_slice = 1, {
  grid.rect(unit(2.05,'npc'), unit(1,'npc'),
            width = (length(co[[1]]) + length(co[[2]]))/length(co[[1]]) * unit(0.13,'npc') + unit(1,'mm'),
            height = (length(ro[[1]]) + length(ro[[2]]))/length(ro[[1]]) * unit(0.245,'npc') + unit(1,'mm'),
            gp = gpar(lwd=2, lty=1, fill=NA, col="red"), just=c('left','top')
  )
})
dev.off()

# Supplemental Figure 6: Plots for Comparing Stong and Weak Effector Profiles Activity Profiles
tiff(filename = "Stong and Weak Effector.tiff", height = 8, width = 22, units = "in", res = 300)
cowplot::plot_grid(plot.HepEpox, plot.Paraquat, ncol = 2, nrow = 1)
dev.off()

# Supplemental Figure 7: Amphetamine Active Curve-Fits
tiff(filename = "Amphetamine Active Curve-Fits.tiff", height = 8, width = 9, units = "in", res = 300)
cowplot::plot_grid(title, p1, ncol=1, rel_heights = c(0.1,1))
dev.off()

# Supplemental Figure 8: Diazepam Behavior Data and Active Curve-Fits
tiff(filename = "Diazepam Bhv and Active Fits.tiff", height = 14, width = 13, units = "in", res = 300)
cowplot::plot_grid(SAplot_Dzp, p1, nrow=2)
dev.off()

# Supplemental Figure 9: D-Sorbitol Tim-Series and Active Curve-Fits
tiff(filename = "D-sorbitol Bhv and Active Fits.tiff", height = 14, width = 16, units = "in", res = 300)
cowplot::plot_grid(SAplot_Srb, p1_Srb, nrow=2)
dev.off()

# Supplemental Figure 10: Unity Plot to Compare This Study's Median LOELs to Study Median LOELs
tiff(filename = "Median LOEL Unity Plot.tiff", height = 5, width = 5, units = "in", res=300)
print(unity_plot)
dev.off()
