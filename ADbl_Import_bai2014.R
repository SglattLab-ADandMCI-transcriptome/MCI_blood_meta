## Import data from Bai14
## GCH

require(gcrma)
require(openxlsx)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Bai14"

cat("Reading CEL files from bai2014 aka Bai14\n")
files = list.files(paste0(rawlocation,"blood/bai2014/"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"blood/bai2014/"))

Exprs = exprs(RMA)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = gsub("_HG.U133_Plus_2.CEL.gz","",names(Exprs))

covs = read.xlsx(paste0(rawlocation,"blood/bai2014/NIHMS562483-supplement-Supplemental_Table1.xlsx"), startRow = 13)
covs = covs[,c(1,13,4,6,3,1)]
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_tissue = "whole blood"
covs$FACTOR_dx[covs$FACTOR_dx=="NL"] = "CTL"
covs$FACTOR_dx[covs$FACTOR_dx=="AD"] = "AD"
covs$FACTOR_sex[covs$FACTOR_sex=="1"] = "male"
covs$FACTOR_sex[covs$FACTOR_sex=="2"] = "female"
covs$FACTOR_ethnicity[covs$FACTOR_ethnicity=="W"] = "white"
covs$FACTOR_ethnicity[covs$FACTOR_ethnicity=="O"] = "other"
covs$Sample_ID = gsub("_HG-U133_Plus_2.CEL.pimg","",covs$Sample_ID)
covs$Sample_ID = gsub("-",".",covs$Sample_ID)
print(head(covs))

rm(RMA)
