## Import data from GSE63060 aka AddNeuroMed1
## GCH

require(data.table)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "AddNeuroMed1"

cat("Reading normalized data from GSE63060 aka AddNeuroMed1\n")

## These data were pre asinh'd with limma
Exprs = fread(paste0(rawlocation,"blood/GSE63060/data/GSE63060_normalized.txt.gz"),fill=T,data.table=F)
##this file has an erroneous extra crlf
Exprs[30881,1] = Exprs[30880,1]
Exprs = Exprs[-30880,]


names(Exprs)[1] = "PROBEID"
## probe/gene translation
## Illumina HumanHT-12 V3.0 expression beadchip
translist = fread(paste0(rawlocation,"blood/GSE63060/GPL6947-13512.txt"),skip=30)
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
covs = fread(paste0(rawlocation,"blood/GSE63060/E-GEOD-63060.sdrf.txt"),data.table=F)
covs = covs[,c(4,17,14,5,6,11)]
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_ethnicity[grepl("Western European",covs$FACTOR_ethnicity)] = "white"
covs$FACTOR_ethnicity[grepl("Other Caucasian",covs$FACTOR_ethnicity)] = "white"
covs$FACTOR_ethnicity[grepl("Unknown",covs$FACTOR_ethnicity)] = "unknown"
covs$FACTOR_tissue[grepl("blood",covs$FACTOR_tissue)] = "whole blood"


#4856076038_D is not in the covs list, excluding.
foo = which(names(Exprs)=="4856076038_D")
Exprs = Exprs[,-foo]

rm(translist)
rm(probes)





