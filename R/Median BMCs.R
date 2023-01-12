## ---------------------------
##
## Script Name: Literature Review Concordance
##
## Purpose of Script: Evaluate the concordance of results from
##                    papers with the gabi results under various
##                    criteria.
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
## Notes:
##
##
## ---------------------------


library(data.table)



# Load data and save as data.table ----------------------------------------


load("~/Desktop/Stephanie Projects/DNT60 Analysis/pipelining/results/Padilla_DNT60_tcpl_out.Rdata")
tcpl_out.dt <- as.data.table(DNT60_tcpl_out)


# For active chemicals calculate 3 median BMC  ----------------------------


