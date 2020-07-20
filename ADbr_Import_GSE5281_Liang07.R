## Import data from GSE5281, Liang07
## GCH

require(data.table)
# require(gcrma)
require(affy)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Liang07"

cat("Reading CEL files from GSE5281 aka Liang07\n")
files = list.files(paste0(rawlocation,"GSE5281/"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
# RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"GSE5281/"))
# Exprs = exprs(RMA)
batch = ReadAffy(filenames = files,celfile.path = paste0(rawlocation,"GSE5281/"))
eset = expresso(batch, 
                bgcorrect.method = "mas",
                normalize.method = "quantiles",
                pmcorrect.method = "mas",
                summary.method = "avgdiff")

Exprs = exprs(eset)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = unlist(lapply(strsplit(names(Exprs),"\\."), `[[`,1))

refs = fread("references/GPL570-55999.txt",data.table=F,skip=16)
rnames = character()
for(i in 1:length(Exprs$PROBEID)){
  rnames[i] = refs$`Gene Symbol`[which(refs$ID == Exprs$PROBEID[i])]
}
rnames = sub(" //.*","",rnames)
Exprs$PROBEID = rnames

rawcovs = readxl::read_excel(paste0(rawlocation,"GSE5281/GSE5281_sample_characteristics.xls"))

covs = data.frame(Sample_ID = rawcovs$`GEO Accession:`,
                  FACTOR_tissue = rawcovs$`Organ Region:`,
                  FACTOR_dx = rawcovs$`Disease State:`,
                  FACTOR_age = rawcovs$`Age:`,
                  FACTOR_sex = rawcovs$`Sex:`,
                  FACTOR_ethnicity = rawcovs$`Ethnicity:`)

covs$FACTOR_age = unlist(lapply(strsplit(as.character(covs$FACTOR_age)," "), `[[`,1))
print(covs$FACTOR_tissue)
covs$FACTOR_tissue = factor(covs$FACTOR_tissue)
levels(covs$FACTOR_tissue) = c("ERC","HIP","TCX","PSC","PVC","FCX")
print(covs$FACTOR_tissue)
print(covs$FACTOR_dx)
covs$FACTOR_dx = factor(covs$FACTOR_dx)
levels(covs$FACTOR_dx) = c("AD","CTL")
print(covs$FACTOR_dx)
print(covs$FACTOR_ethnicity)
covs$FACTOR_ethnicity = factor(covs$FACTOR_ethnicity)
levels(covs$FACTOR_ethnicity) = c("caucasian","unknown")
print(covs$FACTOR_ethnicity)
print(head(covs))

# rm(RMA)
rm(rawcovs)
rm(batch)
rm(eset)
rm(translist)
