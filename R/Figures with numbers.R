## ---------------------------
##
## Script Name: Figures for Manuscript with Numbers
##
## Purpose of Script: Create figures for manuscript with associated numbers so
##                    that a guide can be created for mentors
##
## Author: Zachary Rowson
##
## Date Created: 2022-06-13
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/gabi/Manuscript"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


library(gabi)
library(data.table)
library(ComplexHeatmap)
library(gridExtra)
library(viridis)
library(grid)


load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/raw data/Padilla_DNT60_lmr0_w_egid.Rdata")
load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")
load("~/Desktop/StephanieProjects/DNT60Analysis/results analysis/data/refer_chem.Rdata")
load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/pipelined data/Padilla_DNT60_tcplfits.rda")
load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/pipelined data/Padilla_DNT60_rows_n.rda")
tcpl_out.dt <- as.data.table(DNT60_tcpl_out)


# Figure 1 - Fluoxetine SA plot -------------------------------------------

plot_SA <- function(chemical) {
              # Units
              unit.t = "min"
              unit.mov = "cm"
              unit.conc = paste0("\U03BC","M")
              prsp = "SA"
              no.A = 10

              # Extract chemical data
              group <- unique(lmr0.egid[cpid == chemical, egid])

              ## Identify movement columns of interest
              t.cols <- grep("vt", names(lmr0.egid), value = TRUE)
              cols <- t.cols[(no.A+1):length(t.cols)]
              A.cols <- t.cols[!(t.cols%in%cols)]

              ## extract data to be plotted, exclude acclimation
              to.fit <- lmr0.egid[cpid==chemical | (wllt=="v" & egid==group), -A.cols, with=FALSE]

              # create appropriate axes titles for plots
              label.y <- "Speed"

              # Format data for plotting

              ## calculate mean and 50% CIs for each vector column by concentration group, excluding concentration
              exclude.A <- t.cols[!(t.cols%in%A.cols)]
              means <- to.fit[, lapply(.SD, function(col) mean(col,na.rm=T)),
                              .SDcols = exclude.A,
                              by = conc]

              ## calculate CI's for transformed values then transform back
              shift <- 1
              logCIs <- to.fit[, lapply(.SD, function(x) log10(x+shift)), .SDcols=exclude.A, by=conc][
                , lapply(.SD, function(x) t.test(x,conf.level=0.50)$conf.int), .SDcols=exclude.A, by=conc]
              CIs <- logCIs[, lapply(.SD, function(x) (10^x)-shift), by=conc][,lapply(.SD, function(col) abs(diff(col))/2), .SDcols=exclude.A, by=conc]

              ## elongate means and CIs data, and join
              means_long <- data.table::melt(means, id.vars = "conc", variable.name = "t", value.name = "mean")
              means_long[, t := sub("vt","",t)]
              CIs_long <- data.table::melt(CIs, id.vars = "conc", variable.name = "t", value.name = "CI")
              CIs_long[, t := sub("vt","",t)]
              stats <- means_long[CIs_long, on = c("conc","t")][, conc := as.factor(conc)]
              stats[, t := as.numeric(t)]

              # create standard error of mean estimates by time period and plot as ribbons or error bars

              # plot time-series data

              ## create title, x- and y-axis titles, legend labels, and legend title
              title <- paste0("Sample Averaged Time-Series for ", chemical)
              title.t <- paste0("Time (",unit.t,")")
              title.mean <- paste0("Mean ", label.y, " (",unit.mov,"/",unit.t,")")
              conc.n <- to.fit[wllq==1, .N, by=.(conc)][order(conc)]
              legend.labels <- paste0(conc.n$conc, ", n=", conc.n$N)

              title.legend <- paste0("Concentration (", unit.conc, ")")

              ## get better colors for plotting
              N <- length(unique(to.fit[,conc]))
              colors <- viridis::viridis(N)

              ## create x-axis breaks and labels
              m <- as.integer(max(means_long[,t]))
              x.breaks <- seq(from=no.A,to=m,by=10)
              x.labels1 <- 2*seq(from=no.A,to=m,by=10)
              x.labels2 <- x.labels1 - 2
              x.labels <- paste(x.labels2, x.labels1, sep="-")

              ## plot
              plot <- ggplot() +
                geom_point(data = stats, aes(x=t, y=mean, color=as.factor(conc))) +
                geom_line(data = stats, aes(x=t, y=mean, color=conc, group=conc)) +
                scale_x_continuous(breaks = x.breaks, labels = x.labels) +
                scale_color_manual(values = colors, labels=legend.labels) +
                geom_ribbon(data = stats,
                            aes(x=t, ymax=mean+CI, ymin=mean-CI, group=conc, fill=conc),
                            alpha = 0.25) +
                geom_rect(aes(xmin=10,xmax=30,ymin=-1.25,ymax=-0.5) ,fill="white", color="black") +
                geom_rect(aes(xmin=30,xmax=50,ymin=-1.25,ymax=-0.5) ,fill="black", color="black") +
                annotate("text", x=c(20,40), y=rep(-0.875,2), color=c("black","white"), label=c("Light Phase","Dark Phase"), size=5) +
                scale_fill_manual(values = colors, labels=legend.labels) +
                labs(title = title, subtitle = "Acclimation Period Excluded: 50% Confidence Bands",
                     x = title.t, y = title.mean, color = title.legend) +
                guides(fill = "none") +
                theme_bw() +
                theme(text = element_text(size = 18))

              return(plot)
}

plot_SA(chemical = "Fluoxetine")


# Figure 2 - BMC heatmap --------------------------------------------------


ac <- c("AUC_L","avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "AUC_D", "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "AUC_T","avgS_T","AUC_r")

## Gather BMC's
bmd <- lapply(tcplfits_n, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  bmd <- fit[["summary"]]$bmd
  if (!(hitcall>0.8 & !is.na(bmd) & bmd != 0)) {
    bmd <- 10000
  }
  bmd
})))
bmd <- do.call("rbind", bmd)
bmd <- log10(bmd[, ac])
bmd <- bmd[,-c(1,9,14)]

# Identify active chemicals
actives <- rownames(bmd)[apply(bmd, 1, function(row) any(row < 4))]

# Isolate data of interest
to.fit <- bmd

# PBDE-47 name is too long
rownames(bmd)[6] <- "PBDE-47"
actives[6] <- "PBDE-47"

# Create a matrix specifying if an up, down, up and down, or down and up arrow should be printed in cells.
to.fit <- lapply(tcplfits_n, function(chm) chm[-c(1,9,14)])
layer.mat <- lapply(to.fit, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  fit_method <- fit[["summary"]]$fit_method
  dir <- sign(fit[["summary"]]$top)
  if (hitcall>0.8 & fit_method=="gnls" & dir==1) {
    layer <- 0
  } else if (hitcall>0.8 & fit_method=="gnls" & dir==-1) {
    layer <- 1
  } else if (hitcall>0.8 & dir==1) {
    layer <- 2
  } else if (hitcall>0.8 & dir==-1) {
    layer <- 3
  } else layer <- 4
})))
layer.mat <- do.call("rbind", layer.mat)
rownames(layer.mat)[6] <- "PBDE-47" # Change PBDE-47 name to match BMD matrix
layer.mat <- layer.mat[actives,]


# Column labels
col_labels <- c(expression("Average Speed in Light"^1), "Average Acceleration in Light", "Average Jerk in Light", expression("Habituation in Light"^3),
                "Startle Acceleration", "Startle Relative to Avg. Speed in Light", "Startle Fold-Change",
                expression("Average Speed in Dark"^1), "Average Acceleration in Dark", "Average Jerk in Dark", expression("Habituation in Dark"^3),
                "Average Speed in Both Phases", expression("AUC in Dark / AUC in Light Ratio"^2)) # Superscripts notate references in poster
# Legends

# Custom heat legend.
f2 = circlize::colorRamp2(seq(min(bmd), max(bmd), length = 8), rev(viridis(8)), space = "sRGB")
heat_lgd = Legend(col_fun = f2,
                  title = paste0("BMC log(","\U03BC","M)"),
                  title_position = "lefttop",
                  legend_width = unit(4,"cm"),
                  direction = "horizontal")
# Column annotation by phase legend.
ann_lgd = Legend(labels = c("Light","Transition","Dark","Light+Dark"),
                 title = "Phase",
                 title_position = "leftcenter",
                 labels_gp = gpar(fontsize=12),
                 title_gp = gpar(fontsize=12),
                 legend_gp = grid::gpar(fill = c("white","grey","black","red")),
                 border = TRUE,
                 nrow = 1,
                 column_gap = unit(5, 'mm'))
# Legend for cell signal arrows.
dir_lgd = Legend(labels = c(paste("\U2191","Gain"),
                            paste("\U2193","Loss"),
                            paste("\U21C5","GainLoss"),
                            paste("\U21F5","LossGain")),
                 title = "Signal Direction",
                 labels_gp = gpar(fontsize=12),
                 title_gp = gpar(fontsize=12),
                 title_position = "leftcenter",
                 nrow = 1,
                 column_gap = unit(0, 'mm')) # Will produce warnings, don't worry.
lgd_list <- list(ann_lgd, dir_lgd)

# Create column annotation indicating the phase of the LMR that is described.
column_ha <- ComplexHeatmap::HeatmapAnnotation(Phase = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                              levels=c("Light","Transition","Dark","Light+Dark")),
                                               border = TRUE,
                                               col = list(Phase=c("Light"="white",
                                                                  "Transition"="grey",
                                                                  "Dark"="black",
                                                                  "Light+Dark"="red")),
                                               annotation_legend_param = list(nrow = 1),
                                               show_annotation_name = FALSE,
                                               show_legend = FALSE)

# Create function to add arrows indicating signal direction in cells.
cell_fun <- function(j, i, x, y, width, height, fill) {
  if (layer.mat[i,j] == 0) {
    grid.text("\U21C5", x, y)
  } else if (layer.mat[i,j] == 1) {
    grid.text("\U21F5", x, y)
  } else if (layer.mat[i,j] == 2) {
    grid.text("\U2191", x, y)
  } else if (layer.mat[i,j] == 3) {
    grid.text("\U2193", x, y)
  }
}


# Create main heat map.
htlist <- ComplexHeatmap::Heatmap(bmd[actives,],

                                  # Specify some parameters for heat legend.
                                  name = paste0("BMC log(","\U03BC","M)"),
                                  col = f2,
                                  border_gp = grid::gpar(col="black",lwd=1),
                                  rect_gp=grid::gpar(col="grey"),
                                  show_heatmap_legend = TRUE,
                                  heatmap_legend_param = list(legend_height = unit(6,"cm"),
                                                              direction = "vertical",
                                                              title_gp = grid::gpar(fontsize=14),
                                                              labels_gp = grid::gpar(fontsize = 14)),

                                  # Heatmap width and height
                                  # width = unit(16,"in"),
                                  # height = unit(18,"in"),

                                  # Column label parameters
                                  column_labels = col_labels, column_names_rot = 45,


                                  # Split columns by phase.
                                  column_split = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                        levels=c("Light","Transition","Dark","Light+Dark")),

                                  # Specify some column parameters.
                                  column_title = NULL,
                                  top_annotation = column_ha,

                                  # Add signal direction arrows.
                                  cell_fun = cell_fun,

                                  # Specify some parameters row dendrogram aesthetics and row labels.
                                  row_title_side = "right",
                                  row_title_rot = 0,
                                  row_split = 5,
                                  row_dend_side = "right",
                                  row_names_side = "left",

                                  # Clustering parameters.
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE,
                                  clustering_distance_rows = "pearson",

                                  # Font sizes
                                  row_names_gp = grid::gpar(fontsize=14),
                                  column_names_gp = grid::gpar(fontsize=14),
                                  row_title_gp = grid::gpar(fontsize=14),
                                  column_title_gp = grid::gpar(fontsize=14)
)

png(filename = "Figures/Fig2 - BMC heatmap.png", width = 22, height = 26, unit = "cm", res=300)
draw(htlist, merge_legend = FALSE, annotation_legend_list = lgd_list,
     annotation_legend_side = "top", align_annotation_legend = "heatmap_center")
dev.off()


# Figure 3 - Diazepam SA plot + endpoints --------------------------------------------


chemical <- "Diazepam"
SAplot_Dzp <- plot_SA(chemical)
row1_Dzp <- rows_n[[chemical]][["avgS_L"]]
row2_Dzp <- rows_n[[chemical]][["RoA_L"]]
row3_Dzp <- rows_n[[chemical]][["avgS_T"]]

fit1_Dzp <- concRespCoreZR(row1_Dzp, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit2_Dzp <- concRespCoreZR(row2_Dzp, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit3_Dzp <- concRespCoreZR(row3_Dzp, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]

p1 <- cowplot::plot_grid(fit1_Dzp, fit2_Dzp, fit3_Dzp, nrow=2)

png(filename = "Figures/Fig3 - Diazepam tSeries + curves.png", width = 40, height = 40, unit = "cm", res=300)
cowplot::plot_grid(SAplot_Dzp, p1, nrow=2)
dev.off()


# Figure 4 - Contrast behavior profiles --------------------------------------------


chemical <- "Paraquat"
SAplot <- plot_SA(chemical)
row1 <- rows_n[[chemical]][["strtlAavg"]]
row2 <- rows_n[[chemical]][["avgS_D"]]
row3 <- rows_n[[chemical]][["hbt_D"]]
row4 <- rows_n[[chemical]][["avgS_T"]]
row5 <- rows_n[[chemical]][["AUC_r"]]

fit1 <- concRespCoreZR(row1, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit2 <- concRespCoreZR(row2, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit3 <- concRespCoreZR(row3, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit4 <- concRespCoreZR(row4, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit5 <- concRespCoreZR(row5, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]

legend <- cowplot::get_legend(fit1)

edit_plot <- function(plot) {
  title <- plot$labels$title
  plot$labels$title <- gsub(paste0(chemical," for "), "", title)
  plot$labels$subtitle <- NULL
  plot <- plot + theme(legend.position="none")
  return(plot)
}

plots <- lapply(list(fit1, fit2, fit3, fit4, fit5), edit_plot)
plots[[4]]$labels$title <- gsub(" experimental", "\nexperimental", plots[[4]]$labels$title)
plots[[1]]$labels$title <- gsub(" versus", "\nversus", plots[[1]]$labels$title)

p1 <- cowplot::plot_grid(plotlist=plots, nrow=2)
p2 <- cowplot::plot_grid(p1, legend, rel_widths = c(3, .4))

png(filename = "Figures/Fig4 - Compare1.png", width = 40, height = 40, unit = "cm", res=300)
cowplot::plot_grid(SAplot, p2, nrow=2)
dev.off()

chemical <- "Heptachlor epoxide"
SAplot <- plot_SA(chemical)
row <- rows_n[[chemical]][["AUC_r"]]

fit <- concRespCoreZR(row, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]

png(filename = "Figures/Fig4 - Compare2.png", width = 40, height = 40, unit = "cm", res=300)
cowplot::plot_grid(SAplot, fit, nrow=2)
dev.off()


# Figure 5 - Chlorpyrifos + endpoints -----------------------------


chemical <- "Chlorpyrifos (ethyl)"
SAplot <- plot_SA(chemical)
row1 <- rows_n[[chemical]][["strtlA"]]
row2 <- rows_n[[chemical]][["strtlAavg"]]
row3 <- rows_n[[chemical]][["avgA_D"]]
row4 <- rows_n[[chemical]][["hbt_D"]]

fit1 <- concRespCoreZR(row1, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit2 <- concRespCoreZR(row2, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit3 <- concRespCoreZR(row3, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit4 <- concRespCoreZR(row4, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]

p1 <- cowplot::plot_grid(fit1, fit2, fit3, fit4, nrow=2)

png(filename = "Figures/Fig5 - CPF.png", width = 40, height = 40, unit = "cm", res=300)
cowplot::plot_grid(SAplot, p1, nrow=2)
dev.off()


# Figure 6 - Areas of Additional Activity ---------------------------------


isolate <- actives[order(actives)]
isolate <- c(isolate[c(1,2,6,8,21)], isolate[c(10,12,13,20,22)])

layer.mat2 <- layer.mat[isolate,]

cell_fun1 <- function(j, i, x, y, width, height, fill) {
  if (layer.mat2[i,j] == 0) {
    grid.text("\U21C5", x, y, gp=gpar(fontsize=12))
  } else if (layer.mat2[i,j] == 1) {
    grid.text("\U21F5", x, y, gp=gpar(fontsize=12))
  } else if (layer.mat2[i,j] == 2) {
    grid.text("\U2191", x, y, gp=gpar(fontsize=12))
  } else if (layer.mat2[i,j] == 3) {
    grid.text("\U2193", x, y, gp=gpar(fontsize=12))
  }
}

htlist2 <- ComplexHeatmap::Heatmap(bmd[isolate,],

                                   # Specify some parameters for heat legend.
                                   name = paste0("BMC log(","\U03BC","M)"),
                                   col = f2,
                                   border_gp = grid::gpar(col="black",lwd=1),
                                   rect_gp=grid::gpar(col="grey"),
                                   show_heatmap_legend = TRUE,
                                   heatmap_legend_param = list(legend_height = unit(3.5,"cm"),
                                                               direction = "vertical",
                                                               title_gp = grid::gpar(fontsize=10),
                                                               labels_gp = grid::gpar(fontsize=10)),

                                   # Heatmap width and height
                                   # width = unit(16,"in"),
                                   # height = unit(18,"in"),

                                   # Column label parameters
                                   column_labels = col_labels, column_names_rot = 45,


                                   # Split columns by phase and row by activity type
                                   column_split = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                         levels=c("Light","Transition","Dark","Light+Dark")),
                                   # row_split = factor(c(rep("1",5), rep("2",5))),

                                   # Specify some column parameters.
                                   column_title = NULL,
                                   top_annotation = column_ha,

                                   # Add signal direction arrows.
                                   cell_fun = cell_fun1,

                                   # Specify some parameters row dendrogram aesthetics and row labels.
                                   row_title = NULL,
                                   row_title_rot = 0,
                                   #row_split = 3,
                                   row_names_side = "left",

                                   # Clustering parameters.
                                   cluster_columns = FALSE,
                                   cluster_rows = FALSE,

                                   # Font sizes
                                   row_names_gp = grid::gpar(fontsize=8),
                                   column_names_gp = grid::gpar(fontsize=10),
                                   row_title_gp = grid::gpar(fontsize=10),
                                   column_title_gp = grid::gpar(fontsize=10)
)

ht <- draw(htlist2, merge_legend = FALSE, annotation_legend_list = lgd_list,
           annotation_legend_side = "top", align_annotation_legend = "heatmap_center")
ro <- row_order(ht)
co <- column_order(ht)

# Highlight startle endpoints and average acceleration in Dark endpoint
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



# Figure 7 - D-Sorbitol Active Endpoints ------------------------------------------------------------------


chemical <- "D-sorbitol"

row1 <- rows_n[[chemical]][["avgA_L"]]
row2 <- rows_n[[chemical]][["strtlA"]]
row3 <- rows_n[[chemical]][["strtlAavg"]]

fit1 <- concRespCoreZR(row1, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit2 <- concRespCoreZR(row2, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]
fit3 <- concRespCoreZR(row3, do.plot = TRUE, verbose.plot = FALSE)[["plot"]]

legend <- cowplot::get_legend(fit1)
title <- cowplot::ggdraw() +
  cowplot::draw_label("D-Sorbitol Exposure: Active Curve-Fits",
                      x = 0, hjust = 0, size=24) +
  theme(plot.margin = margin(0,0,0,7))

edit_plot <- function(plot) {
  title <- plot$labels$title
  plot$labels$title <- gsub(paste0(chemical," for "), "", title)
  plot$labels$subtitle <- NULL
  plot <- plot + theme(legend.position="none")
  return(plot)
}

plots <- lapply(list(fit1, fit2, fit3), edit_plot)
plots[[4]]$labels$title <- gsub(" experimental", "\nexperimental", plots[[4]]$labels$title)
plots[[1]]$labels$title <- gsub(" versus", "\nversus", plots[[1]]$labels$title)

p1 <- cowplot::plot_grid(plotlist=plots, nrow=2)
p2 <- cowplot::plot_grid(p1, legend, rel_widths = c(3, .4))

png(filename = "Figures/Fig7 - DSorbitol.png", width = 35, height = 30, unit = "cm", res=300)
cowplot::plot_grid(title, p2, ncol=1, rel_heights = c(0.1,1))
dev.off()

# Figure 8 - Amphetamine Active Endpoints ------------------------------------------------------------------


