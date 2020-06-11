## Import data from GSE36980 aka Hokama13
## GCH

require(data.table)
require(oligo)
require(GEOquery)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Hokama13"

cat("Reading CEL files from GSE36980 aka Hokama13\n")
files = list.files(paste0(rawlocation,"GSE36980/"),".CEL",full.names = T)
## RMA is similar to log2 and normalize = T will quantile normalize
batch = read.celfiles(filenames = files)
RMA = oligo::rma(batch,normalize = T)

Exprs = exprs(RMA)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
names(Exprs) = unlist(lapply(strsplit(names(Exprs),"\\."), `[[`,1))

rawcovs = getGEO(filename = paste0(rawlocation,"GSE36980/GSE36980_series_matrix.txt.gz"))

Exprs$PROBEID = unlist(lapply(strsplit(rawcovs@featureData@data$gene_assignment," "), `[`, 3))

covs = data.frame(Sample_ID = rawcovs@phenoData@data$geo_accession,
                  FACTOR_tissue = rawcovs@phenoData@data$`tissue:ch1`,
                  FACTOR_dx = rawcovs@phenoData@data$title,
                  FACTOR_age = rawcovs@phenoData@data$`age:ch1`,
                  FACTOR_sex = rawcovs@phenoData@data$`Sex:ch1`,
                  FACTOR_ethnicity = "Japanese")

print(covs$FACTOR_dx)
dx = as.character(covs$FACTOR_dx)
dxAD = grepl("^AD_",dx)
dxCTL = grepl("non-AD_",dx)
dx[dxAD] = "AD"
dx[dxCTL] = "CTL"
covs$FACTOR_dx = dx
print(covs$FACTOR_dx)
print(covs$FACTOR_tissue)
covs$FACTOR_tissue = factor(covs$FACTOR_tissue)
levels(covs$FACTOR_tissue) = c("FCX","HIP","TCX")
print(covs$FACTOR_tissue)
print(covs$FACTOR_sex)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
levels(covs$FACTOR_sex) = c("female","male")
print(covs$FACTOR_sex)
print(head(covs))

rm(RMA)
rm(rawcovs)
rm(batch)
