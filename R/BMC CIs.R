## ---------------------------
##
## Script Name: BMC CI widths and conglomerate potency estimate
##
## Purpose of Script: Evaluate the range of BMCs associated with a chemical and the CIs
##                    associated with individual estimates to make a statement about
##                    the width of CIs and our confidence in the produced BMC as well
##                    our confidence in the conglomerate BMC estimate.
##
## Author: Zachary Rowson
##
## Date Created: 2022-06-30
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


library(data.table)
library(ggplot2)


load("~/Desktop/StephanieProjects/DNT60Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")


# Do confidence intervals for BMCs of active chemicals overlap? -----------


tcpl_out.dt <- as.data.table(DNT60_tcpl_out)

to.plot <- tcpl_out.dt[hitcall>0.8 & endp%in%c("avgS_L","avgA_L","avgJ_L","hbt_L",
                                                 "strtlA","strtlAavg","strtlF",
                                                 "avgS_D","avgA_D","avgJ_D","hbt_D",
                                                 "avgS_T","AUC_r")]
to.plot[, endp := factor(endp, levels=rev(c("avgS_L","avgA_L","avgJ_L","hbt_L",
                                        "strtlA","strtlAavg","strtlF",
                                        "avgS_D","avgA_D","avgJ_D","hbt_D",
                                        "avgS_T","AUC_r")))]
to.plot[name=="Polybrominated diphenyl ether (PBDE)-47", name:="PBDE-47"]
ggplot(to.plot, aes(x=log10(bmd),y=endp)) +
  geom_point() +
  geom_errorbar(aes(xmin=log10(bmdl), xmax=log10(bmdu))) +
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2)) +
  labs(title = "Overlap of BMC Estimates for Chemicals", x = "log(BMC)", y = "Endpoint") +
  facet_wrap(vars(name))
