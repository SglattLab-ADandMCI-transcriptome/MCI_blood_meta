## Import data from GSE63061 aka AddNeuroMed2
## GCH

require(data.table)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "AddNeuroMed2"

cat("Reading normalized data from GSE63061 aka AddNeuroMed2\n")

## These data were pre asinh'd with limma
Exprs = fread(paste0(rawlocation,"blood/GSE63061/data/GSE63061_normalized.txt.gz"),fill=T,data.table=F)
names(Exprs)[1] = "PROBEID"
## probe/gene translation
## Illumina HumanHT-12 V4.0 expression beadchip
translist = fread(paste0(rawlocation,"blood/GSE63061/GPL10558-50081.txt"),skip=30,data.table=F)
q = length(Exprs$PROBEID)
cat("Replacing probe IDs with gene symbols.\n")
translist = data.frame(translist$Symbol,translist$ID,stringsAsFactors = F)
names(translist) = c("Symbol","ID")
probes = Exprs$PROBEID
for(foo in 1:q){
  if(probes[foo] %in% translist$ID){
    bar = translist$Symbol[probes[foo] == translist$ID]
    bar = strsplit(bar," ")[[1]][1]
    probes[foo] = bar
  }
  cat(foo,"/",q,"\r")
}
Exprs$PROBEID = probes
## already quantile normalized
## lumi does arcsinh
names(Exprs) = gsub("X","",names(Exprs))


## covariates
covs = fread(paste0(rawlocation,"blood/GSE63061/E-GEOD-63061.sdrf.txt"),data.table=F)
covs = covs[,c(4,19,16,5,6,13)]
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")
## This study has a few "oddball" notes for diagnosis.  The updated version of the diagnosis
## is used here, because it is more likely to be accurate in terms of underlying genetics and such,
## but this will have to be noted in the manuscript.
covs$FACTOR_dx[grepl("OTHER",covs$FACTOR_dx)] = "Other"
covs$FACTOR_dx[grepl("borderline MCI",covs$FACTOR_dx)] = "MCI"
covs$FACTOR_dx[grepl("MCI to CTL",covs$FACTOR_dx)] = "CTL"
covs$FACTOR_dx[grepl("CTL to AD",covs$FACTOR_dx)] = "AD"
covs$FACTOR_race[grepl("Western European",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("White_And_Asian",covs$FACTOR_race)] = "other"
covs$FACTOR_race[grepl("Other Caucasian",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("English",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Irish",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Any_Other_White",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Welsh",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Scottish",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("white",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Asian",covs$FACTOR_race)] = "asian"
covs$FACTOR_race[grepl("Indian",covs$FACTOR_race)] = "indian"
covs$FACTOR_race[grepl("Black",covs$FACTOR_race)] = "black"
covs$FACTOR_race[grepl("Scottish",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Unknown",covs$FACTOR_race)] = "unknown"
covs$FACTOR_race[grepl("Other",covs$FACTOR_race)] = "other"
covs$FACTOR_race[grepl("British",covs$FACTOR_race)] = "white"
covs$FACTOR_race[grepl("Carib",covs$FACTOR_race)] = "caribbean"
covs$FACTOR_tissue[grepl("blood",covs$FACTOR_tissue)] = "whole_blood"

rm(translist)
rm(probes)





