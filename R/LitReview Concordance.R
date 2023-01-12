## ---------------------------
##
## Script Name: Literature Review Concordance
##
## Purpose of Script: Evaluate the concordance of gabi results
##                    with results from other studies.
##
## Author: Zachary Rowson
##
## Date Created: 2022-07-22
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/gabi/Manuscript"
##
## ---------------------------
##
## Notes: This comparison will be useful for assessing the results of analysis in terms of other
##        laboratory results. Do chemicals found active tend to be active in other papers and do
##        papers inactive tend to be inactive? Will also help assess the effect that including
##        more endpoints has on detection of chemical activity, are we seeing activity in chemicals
##        widely active in other studies or are we potentially seeing activity previously missed or spurious activity
##        using these endpoints. Also, identifying chemicals widely active in other studies but not active
##        in the gabi analysis may provide insight into the lack of assay sensitivity.
##
##
## ---------------------------


library(data.table)
library(readxl)

# Load tcplFit2 output data and save as data.table

load("~/Desktop/Stephanie Projects/DNT60 Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")
tcplOut.dt <- as.data.table(DNT60_tcpl_out)

# Load literature review table and save sheet of interest as data.table

litReview <- as.data.table( readxl::read_xlsx("tables/copy_Chemical Comparison Tables_Behavior Results and  Literature_FINAL-rev5_2_22-jko.xlsx",
                                              sheet = "gabiComparison") )
View(litReview)


# Format litReview table --------------------------------------------------

# Set column names and row names
newNames <- c( "cpid", sapply(names(litReview), gsub, pattern=" ", replacement=".")[-1] )
setnames(litReview, old = names(litReview), new = newNames)
litReview[1, cpid := "5,5-Diphenylhydantoin"]
litReview[12, cpid := "Bis(tributyltin) Oxide"]
litReview[52, cpid := "Polybrominated diphenyl ether (PBDE)-47"]

# ID chemicals active in net activity metrics
active.1 <- tcplOut.dt[endp%in%c("avgS_L","avgS_D","avgS_T") & hitcall>0.8, unique(name)]

# ID chemicals active in other endpoints
active.2 <- tcplOut.dt[!endp%in%c("avgS_L","avgS_D","avgS_T","AUC_L","AUC_D","AUC_T") & hitcall>0.8,
                       unique(name)]

# Edit litReview columns
litReview[is.na(Active.in.gabi.analysis), Active.in.gabi.analysis := 0]
litReview[cpid %in% active.1, active.in.Net.Activity.metrics := 1]
litReview[Active.in.gabi.analysis==1 & is.na(active.in.Net.Activity.metrics),
          active.in.Net.Activity.metrics := 0]
litReview[cpid %in% active.2, Active.in.other.endpoints := 1]
litReview[Active.in.gabi.analysis==1 & is.na(Active.in.other.endpoints),
          Active.in.other.endpoints := 0]

str(litReview)
litReview[, Active.in.gabi.analysis := as.logical(Active.in.gabi.analysis) ]
litReview[, Number.of.publications.with.activity := as.integer(Number.of.publications.with.activity)]
litReview[, `%.of.papers.with.activity` := as.numeric(`%.of.papers.with.activity`)]

litReview[Number.of.asociated.publications==0,
          `:=` (Number.of.publications.with.activity=NA, `%.of.papers.with.activity`=NA)]


# Questions ---------------------------------------------------------------


# For the 21 chemicals found active, how much concordance are we seeing from other labs?
litReview[Active.in.gabi.analysis==TRUE, mean(`%.of.papers.with.activity`,na.rm=TRUE)]
# On average we are seeing that papers associated with active chemicals are 74.44% active.

# Does this differ for inactive chemicals?
litReview[Active.in.gabi.analysis==FALSE, mean(`%.of.papers.with.activity`,na.rm=TRUE)]
# On average we are seeing that papers associated with active chemicals are 57.57% active.
# So, papers associated with chemicals active in gabi analysis tend to have more associated
# papers that found activity.
# Would it be better to just do some binomial like thing?
litReview[Active.in.gabi.analysis==TRUE, sum(Number.of.publications.with.activity,na.rm=T) / sum(Number.of.asociated.publications)]
litReview[Active.in.gabi.analysis==F, sum(Number.of.publications.with.activity,na.rm=T) / sum(Number.of.asociated.publications)]

# Do the number of percent chemicals without associated papers differ from inactive to
# active chemicals.
n_1 <- litReview[Active.in.gabi.analysis==TRUE & Number.of.asociated.publications==0, .N]
litReview[Active.in.gabi.analysis==TRUE, n_1 / .N]
n_2 <- litReview[Active.in.gabi.analysis==FALSE & Number.of.asociated.publications==0, .N]
litReview[Active.in.gabi.analysis==FALSE, n_2 / .N]
# Active chemicals have no associated publications less often than inactive chemicals
# (32.82% v. 44.75%), though if this difference is statistically significant hasn't been
# established.

# How does active concordance in studies differ between chemicals found active in the
# commonly used net activity metric vs. those found active in other endpoints and those
# active in both categories?
litReview[active.in.Net.Activity.metrics==T & Active.in.other.endpoints==F,
          mean(`%.of.papers.with.activity`,na.rm=TRUE)]
litReview[active.in.Net.Activity.metrics==F & Active.in.other.endpoints==T,
          mean(`%.of.papers.with.activity`,na.rm=TRUE)]
litReview[active.in.Net.Activity.metrics==T & Active.in.other.endpoints==T,
          mean(`%.of.papers.with.activity`,na.rm=TRUE)]
# Chemicals active in net activity metrics but not other endpoints (only one chemical, BPA)
# were active in 83.33% of associated papers. Chemicals active in other endpoints were active
# in 68.75% of associated papers on average. Chemicals active in net activity and other
# other endpoints were active in 80.56% of papers on average.

# To get at Katie's question of whether evaluation of data using these other endpoints results
# in detection of chemical activity where chemical effect may have been missed, look at chemicals
# active in other endpoints only with a large number of associated papers.
litReview[active.in.Net.Activity.metrics==F & Active.in.other.endpoints==T,]
# There are no chemicals with many associated publications but little activity.
# Diethylstilbesterol we see 2 total publications both with no activity.
# There are some chemicals that have no associated publications.
litReview[active.in.Net.Activity.metrics==F & Active.in.other.endpoints==T &
            Number.of.asociated.publications==0,]
# 5,5-Diphenylhydantoin, Heptachlor epoxide, and Triethyltin
# May be that Heptachlor epoxide is assumed to be too teratogenic? Triethyltin as well?
# 5,5-Diphenylhydantoin is interesting.
# Obviously there may be publications that tested these chemicals that were not found.

# Were there any chemicals that were widely tested and widely found active in those other studies?
# Valproate is the most extreme example where it was found active in 6 of 7 associated publications.
# (I think valproate was a near miss in some endpoints). Could the lack of activity in gabi be due
# to it's teratogenic properties? (Isn't is teratogenic?)
# THis is probably amateur , but could valproate be used as an example of what this assay can be used
# for, the screening of chemicals for behavioral effects in the absence of morphological effect.
# Valproate is known to be teratogenic in humans but also neurodevelopmentally toxic. The utility
# of this assay may be identifiying behavioral effects in the absence of moprhological effect which
# is likely more rare. May help explain the lack of sensitivity.

# We also see activity in many publications for lead acetate and acetaminophen. Is the same argument
# possible?

# Inactive chemicals with 3 associated publications with 2 or 3 papers with activity are Haloperidol,
# Nicotine, and Carbamezipine.
