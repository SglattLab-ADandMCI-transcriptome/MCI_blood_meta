## Import data from GSE48350 aka Berchtold14
## GCH

require(data.table)
# require(gcrma)
require(affy)
require(GEOquery)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Berchtold14"

cat("Reading CEL files from GSE48350 aka Berchtold14\n")
files = list.files(paste0(rawlocation,"GSE48350/"),".CEL",full.names = F)
## RMA is similar to log2 and normalize = T will quantile normalize
# RMA = justGCRMA(filenames = files,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"GSE48350/"))
# Exprs = exprs(RMA)
batch = ReadAffy(filenames = files,celfile.path = paste0(rawlocation,"GSE48350/"))
eset = expresso(batch, 
                bgcorrect.method = "mas",
                normalize.method = "quantiles",
                pmcorrect.method = "mas",
                summary.method = "avgdiff")

Exprs = exprs(eset)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = unlist(lapply(strsplit(names(Exprs),"\\.|_"), `[[`,1))

rawcovs = getGEO(filename = paste0(rawlocation,"GSE48350/GSE48350_series_matrix.txt.gz"))

Exprs$PROBEID = unlist(lapply(strsplit(rawcovs@featureData@data$`Gene Symbol`," "), `[`, 1))

covs = data.frame(Sample_ID = rawcovs@phenoData@data$geo_accession,
                  FACTOR_tissue = rawcovs@phenoData@data$`brain region:ch1`,
                  FACTOR_dx = "temp",
                  FACTOR_age = rawcovs@phenoData@data$`age (yrs):ch1`,
                  FACTOR_sex = rawcovs@phenoData@data$`gender:ch1`,
                  FACTOR_ethnicity = "unknown")

dx = paste(rawcovs@phenoData@data$`individual:ch1`,rawcovs@phenoData@data$source_name_ch1)
dxAD = grepl("(, AA)|(_AD)",dx)
dxCTL = grepl(", C",dx)
dx[dxAD] = "AD"
dx[dxCTL] = "CTL"
covs$FACTOR_dx = dx
print(covs$FACTOR_dx)
print(covs$FACTOR_tissue)
covs$FACTOR_tissue = factor(covs$FACTOR_tissue)
levels(covs$FACTOR_tissue) = c("ERC","HIP","PCG","PCG","FCX")
print(covs$FACTOR_tissue)
print(head(covs))

# rm(RMA)
rm(rawcovs)
rm(batch)
rm(eset)
rm(translist)
