## ---------------------------
##
## Script name: SA LMR Data Plots
##
## Purpose of script: Produce SA LMR data plots for chemicals found active in
##                    gabi analysis. Plots will be used for supplemental files
##                    for the paper.
##
## Author: Zachary Rowson
##
## Date Created: 2022-05-26
##
## Working Directory: "L:/PRIV/NHEERL_Xena/NHEERL_Xena/Zach Rowson/Behav Protocol/Data Analysis"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


library(ggplot2)
library(data.table)
library(magrittr)


unit.t = "min"
unit.mov = "cm"
unit.conc = paste0("\U03BC","M")
prsp = "SA"
no.A = 10

load("~/Desktop/Stephanie Projects/DNT60 Analysis/pipelining/raw data/Padilla_DNT60_pmr0_egid.Rdata")
load("~/Desktop/Stephanie Projects/DNT60 Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")


# Plot for active chemicals  ------------------------------------------------


## Identify active chemicals
chemicals <- as.data.table(DNT60_tcpl_out)[hitcall>0.8, unique(name)]
names(chemicals) <- chemicals

chemicals1 <- c("Bis(tributyltin) Oxide", "Dieldrin", "Heptachlor epoxide",
                "Phenobarbital", "Triethyltin")

plots1 <- lapply(chemicals1, function(name) {

                    group <- unique(data[cpid == name, egid])

                    ## Identify movement columns of interest
                    t.cols <- grep("vt", names(data), value = TRUE)
                    cols <- t.cols[(no.A+1):length(t.cols)]
                    A.cols <- t.cols[!(t.cols%in%cols)]

                    ## extract data to be plotted, exclude acclimation
                    to.fit <- data[cpid==name | (wllt=="v" & egid==group), -A.cols, with=FALSE]

                    # create appropriate axes titles for plots
                    label.y <- "Speed"

                    # format data for plotting

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
                    title <- paste0("Sample Averaged Time-Series for ", name)
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
                              annotate("text", x=c(20,40), y=rep(-0.875,2), color=c("black","white"), label=c("Light Phase","Dark Phase"), size=4) +
                              scale_fill_manual(values = colors, labels=legend.labels) +
                              labs(title = title, subtitle = "Acclimation Period Excluded: 50% Confidence Bands",
                                   x = title.t, y = title.mean, color = title.legend) +
                              guides(fill = "none") +
                              theme_bw() +
                              theme(text = element_text(size = 16))

                    plot
                    }
       )


# Save images -------------------------------------------------------------


dir.create("SA Behavior Profiles")
setwd("SA Behavior Profiles")

save(plots, file="SA Plots of Actives.Rdata")

lapply(chemicals, function(name) {
  png(filename = paste0("SA plot ",name,".png"),
      width=720, height=480)
  print(plots[[name]])
  dev.off()
  })


