## Import data from GSE18309 aka Chen09
## GCH

# require(gcrma)
require(affy)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Chen09"

cat("Reading CEL files from GSE18309 aka Chen09\n")

files = list.files(paste0(rawlocation,"blood/GSE18309/data"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
# RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"blood/GSE18309/data/"))
# Exprs = exprs(RMA)
batch = ReadAffy(filenames = files,celfile.path = paste0(rawlocation,"blood/GSE18309/data/"))
eset = expresso(batch, 
                bgcorrect.method = "mas",
                normalize.method = "quantiles",
                pmcorrect.method = "mas",
                summary.method = "avgdiff")

Exprs = exprs(eset)
Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = gsub(".CEL.gz","",names(Exprs))

## probe/gene translation
## [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
translist = fread(paste0(rawlocation,"blood/GSE18309/GPL570-55999.txt"),skip=16)
q = length(Exprs$PROBEID)
cat("Replacing probe IDs with gene symbols.\n")
for(foo in 1:q){
  if(Exprs$PROBEID[foo] %in% translist$ID){
    bar = translist$`Gene Symbol`[Exprs$PROBEID[foo] == translist$ID]
    bar = strsplit(bar," ")[[1]][1]
    Exprs$PROBEID[foo] = bar
  }
  cat(foo,"/",q,"\r")
}


## covariates
covs = fread(paste0(rawlocation,"blood/GSE18309/E-GEOD-18309.sdrf.txt"),data.table=F)
covs = covs[,c(1,7,1,1,4,6)] #age and gender missing
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")
covs$FACTOR_age = NA
covs$FACTOR_sex = NA
covs$FACTOR_race = "asian"
covs$FACTOR_tissue = "whole blood"
covs$FACTOR_dx[covs$FACTOR_dx=="normal"] = "CTL"
covs$FACTOR_dx[covs$FACTOR_dx=="Alzheimers disease"] = "AD"
covs$FACTOR_dx[covs$FACTOR_dx=="mild cognitive impairment"] = "MCI"

# rm(RMA)
rm(batch)
rm(eset)
rm(translist)


