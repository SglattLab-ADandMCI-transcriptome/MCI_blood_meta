## Import data from Bai14
## GCH

# require(gcrma)
require(affy)
require(openxlsx)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Bai14"

cat("Reading CEL files from bai2014 aka Bai14\n")
files = list.files(paste0(rawlocation,"blood/bai2014/"),".CEL",full.names = F)
# ## RMA is similar to log2 and normalize = T will quantile normalize
# RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"blood/bai2014/"))
# Exprs = exprs(RMA)
batch = ReadAffy(filenames = files,celfile.path = paste0(rawlocation,"blood/bai2014/"))
eset = expresso(batch, 
                 bgcorrect.method = "mas",
                 normalize.method = "quantiles",
                 pmcorrect.method = "mas",
                 summary.method = "avgdiff")

Exprs = exprs(eset)
Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = gsub("_HG.U133_Plus_2.CEL.gz","",names(Exprs))

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


covs = read.xlsx(paste0(rawlocation,"blood/bai2014/NIHMS562483-supplement-Supplemental_Table1.xlsx"), startRow = 13)
covs = covs[,c(1,13,4,6,3,1)]
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")
covs$FACTOR_tissue = "whole blood"
covs$FACTOR_dx[covs$FACTOR_dx=="NL"] = "CTL"
covs$FACTOR_dx[covs$FACTOR_dx=="AD"] = "AD"
covs$FACTOR_sex[covs$FACTOR_sex=="1"] = "male"
covs$FACTOR_sex[covs$FACTOR_sex=="2"] = "female"
covs$FACTOR_race[covs$FACTOR_race=="W"] = "white"
covs$FACTOR_race[covs$FACTOR_race=="O"] = "other"
covs$Sample_ID = gsub("_HG-U133_Plus_2.CEL.pimg","",covs$Sample_ID)
covs$Sample_ID = gsub("-",".",covs$Sample_ID)
print(head(covs))

# rm(RMA)
rm(batch)
rm(eset)
rm(translist)
