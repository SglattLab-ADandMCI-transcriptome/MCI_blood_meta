## Import data from GSE28146 aka Blalock11
## GCH

require(data.table)
# require(gcrma)
require(affy)
require(GEOquery)
require(forcats)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Blalock11"

cat("Reading CEL files from GSE28146 aka Blalock11\n")
files = list.files(paste0(rawlocation,"GSE28146/"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
# RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"GSE28146/"))
# Exprs = exprs(RMA)
batch = ReadAffy(filenames = files,celfile.path = paste0(rawlocation,"GSE28146/"))
eset = expresso(batch, 
                bgcorrect.method = "mas",
                normalize.method = "quantiles",
                pmcorrect.method = "mas",
                summary.method = "avgdiff")

Exprs = exprs(eset)
Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = unlist(lapply(strsplit(names(Exprs),"\\."), `[[`,1))

rawcovs = getGEO(filename = paste0(rawlocation,"GSE28146/GSE28146_series_matrix.txt.gz"))

Exprs$PROBEID = unlist(lapply(strsplit(rawcovs@featureData@data$`Gene Symbol`," "), `[`, 1))

covs = data.frame(Sample_ID = rawcovs@phenoData@data$geo_accession,
                  FACTOR_tissue = "HIP",
                  FACTOR_dx = rawcovs@phenoData@data$characteristics_ch1.2,
                  FACTOR_age = rawcovs@phenoData@data$characteristics_ch1.1,
                  FACTOR_sex = rawcovs@phenoData@data$characteristics_ch1,
                  FACTOR_ethnicity = "unknown")

covs$FACTOR_age = unlist(lapply(strsplit(as.character(covs$FACTOR_age)," "), `[[`,2))
print(covs$FACTOR_dx)
covs$FACTOR_dx = fct_collapse(covs$FACTOR_dx,
                              CTL = c("disease status: Control"),
                              AD = c("disease status: Incipient",
                                     "disease status: Moderate",
                                     "disease status: Severe"))
print(covs$FACTOR_dx)
print(covs$FACTOR_sex)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
levels(covs$FACTOR_sex) = c("female","male")
print(covs$FACTOR_sex)
print(head(covs))

# rm(RMA)
rm(batch)
rm(eset)
rm(translist)
rm(rawcovs)
