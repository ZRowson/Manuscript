library(ggplot2)
library(ggupset)
library(tidyverse)


load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")
tcpl_out.dt <- as.data.table(DNT60_tcpl_out)


# Number active per Chemical ----------------------------------------------


tcpl_out.actives <- tcpl_out.dt[hitcall>0.8 & !endp%in%c("AUC_L","AUC_D","AUC_T","AUC_r")]

tcpl_out.actives[endp%in%c("avgS_L","avgS_D","avgS_T"), cat1 := "Average Speed"]
tcpl_out.actives[endp%in%c("avgA_L","avgA_D","avgJ_L","avgJ_D"), cat1 := "Average Acceleration"]
tcpl_out.actives[endp%in%c("hbt_L","hbt_D"), cat1 := "Habituation"]
tcpl_out.actives[endp%in%c("strtlA","strtlAavg","strtlF"), cat1 := "Startle Response"]

tcpl_out.actives[endp%in%c("avgS_L","avgA_L","avgJ_L","hbt_L"), cat2 := "Light"]
tcpl_out.actives[endp%in%c("avgS_D","avgA_D","avgJ_D","hbt_D"), cat2 := "Dark"]
tcpl_out.actives[endp%in%c("strtlA","strtlAavg","strtlF"), cat2 := "Transition"]
tcpl_out.actives[endp%in%c("avgS_T"), cat2 := "Total"]

actv.chem.cat1 <- tcpl_out.actives[, .(cat1 = factor(unique(cat1),levels=c("Average Speed","Average Acceleration",
                                                                           "Habituation","Startle Response"))),
                                   by=.(name)]
actv.chem.cat2 <- tcpl_out.actives[, .(cat2 = factor(unique(cat2),levels=c("Light","Transition","Dark","Total"))),
                                   by=.(name)]

ggplot(actv.chem.cat1, aes(x=cat1, color=cat1, fill=cat1)) +
  geom_bar() +
  labs(title="Number of Chemicals Active in Endpoint Categories",
       x="Endpoint Categories",
       y="Number of Chemicals Active in Category",
       color="Endpoint Category",
       fill="Endpoint Category") +
  scale_y_continuous(breaks=c(2,4,6,8,10)) +
  scale_x_discrete(labels=c("Average Speed","Average Acceleration","Habituation","Startle Response")) +
  scale_color_manual(values=viridis::viridis(4),
                     labels=c("Average Speed","Average Acceleration","Habituation","Startle Response")) +
  scale_fill_manual(values=viridis::viridis(4),
                    labels=c("Average Speed","Average Acceleration","Habituation","Startle Response"))

ggplot(actv.chem.cat2, aes(x=cat2, color=cat2, fill=cat2)) +
  geom_bar() +
  labs(title="Number of Chemicals Active in Experimental Phases",
       x="Experimental Phase",
       y="Number of Chemicals Active in Phase",
       color="Experimental Phase",
       fill="Experimental Phase") +
  scale_y_continuous(breaks=c(2,4,6,8,10,12)) +
  scale_color_manual(values=viridis::viridis(4)) +
  scale_fill_manual(values=viridis::viridis(4))



# Combinations ------------------------------------------------------------


endp.lists <- split(actv.chem.cat1[,cat1], actv.chem.cat1[,name])
to.plot <- data.table(name = names(endp.lists), list = endp.lists)

phase.lists <- split(actv.chem.cat2[,cat2], actv.chem.cat2[,name])
to.plot1 <- data.table(name = names(phase.lists), list = phase.lists)

plot <- ggplot(to.plot, aes(x=list)) +
          geom_bar() +
          scale_x_upset(reverse = TRUE) +
          labs(x = "Endpoint Category Combinations",
               y = "Number of Chemicals with Combination",
               title = "Combinations of Endpoint Activity by Endpoint Category")

plot1 <- ggplot(to.plot1, aes(x=list)) +
          geom_bar() +
          scale_x_upset(reverse = TRUE) +
          labs(x = "Experimental Phase Combinations",
               y = "Number of Chemicals with Combination",
               title = "Combinations of Endpoint Activity by Experimental Phase")

png(filename = "Figures/Figure - Exp Phase.png", width = 30, height = 20, unit = "cm", res=300)
cowplot::plot_grid(plot, plot1, rel_widths = c(9,8),
                   labels = "AUTO", label_size = 24)
dev.off()
