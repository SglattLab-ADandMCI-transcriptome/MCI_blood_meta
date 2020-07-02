## Import data from GSE6613 aka Scherzer07
## GCH

## TODO  #gender, age, ethnicity missing

require(gcrma)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Scherzer07"

cat("Reading CEL files from GSE6613 aka Scherzer07\n")

files = list.files(paste0(rawlocation,"blood/GSE6613/data"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"blood/GSE6613/data/"))

Exprs = exprs(RMA)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = gsub(".CEL.gz","",names(Exprs))


## covariates
covs = fread(paste0(rawlocation,"blood/GSE6613/E-GEOD-6613.sdrf.txt"),data.table=F)
covs = covs[,c(1,4,1,1,1,3)] #gender, age, ethnicity missing
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_sex = NA
covs$FACTOR_age = NA
covs$FACTOR_ethnicity = NA
covs$FACTOR_tissue = "whole blood"
covs$FACTOR_dx[grep("Parkinson",covs$FACTOR_dx)] = "PD"
covs$FACTOR_dx[grep("normal",covs$FACTOR_dx)] = "CTL"
covs$FACTOR_dx[grep("Alzheimer",covs$FACTOR_dx)] = "AD"
covs$FACTOR_dx[!(grepl("PD|CTL|AD",covs$FACTOR_dx))] = "Other"
covs$Sample_ID = gsub("GSE6613","",covs$Sample_ID)
head(covs)

rm(RMA)