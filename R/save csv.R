## ---------------------------
##
## Script Name: Save Manuscript Tables as CSV Files
##
## Purpose of Script: Save tables associated with gabi manuscript as csv files
##
## Author: Zachary Rowson
##
## Date Created: 2023-01-09
## Last Edit: 2022-01-09
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


setwd("~/Desktop/StephanieProjects/gabi/Manuscript/Tables/Tables v8")


# Tables ------------------------------------------------------------------


# Table 2
table2 <- compareLOEL_summary[!is.na(median_bhvLOEL_diff)]
write.csv(table2, file = "table2_compareDiffLOEL.csv")

# Table 3
table3 <- devBhvLOEL.diff[!is.na(median_mag_diff), .(cpid = cpid, median_mag_diff = signif(median_mag_diff,digits=3), sd_mag_diff = signif(sd_mag_diff,digits=3))]
write.csv(table3, file = "table3_compareDiffLOELDev.csv")


# Supplemental Tables -----------------------------------------------------


# Supplemental Table 2
suppTable2 <- data.table::copy(sample.endp.stats)
write.csv(suppTable2, file = "suppTable2_sampleStats.csv")

# Supplemental Table 3
suppTable3 <- data.table::copy(sample.endp.stats_n)
write.csv(suppTable3, file = "suppTable3_sampleStats_n.csv")

# Supplemental Table 4
suppTable4 <- bxcx.params.dt[, .(endp,lam.hat,shift)]
write.csv(suppTable4, file = "suppTable4_boxCoxParams.csv")

# Supplemental Table 5
suppTable5 <- count.hits.endp[order(-N)]
write.csv(suppTable5, file = "suppTable5_hitsPerEndpoint.csv")

# Supplemental Table 6
suppTable6 <- count.hits.chm[order(-N)]
write.csv(suppTable6, file = "suppTable6_hitsPerChemical.csv")

# Supplemental Table 7
suppTable7 <- data.table::copy(strgEffectors_formatted)
write.csv(suppTable7, file = "suppTable7_strEffectors.csv")

# Supplemental Table 8
suppTable8 <- med_BMC.dt
write.csv(suppTable8, file = "suppTable8_medBMCs.csv")

# Supplemental Table 9
suppTable9 <- data.table::copy(litReviewSummary1)
write.csv(suppTable9, file = "suppTable9_litReview.csv")

# Supplemental Table 10
suppTable10 <- data.table::copy(litReviewSummary.subset)
write.csv(suppTable10, file = "suppTable10_litReviewSubset.csv")

