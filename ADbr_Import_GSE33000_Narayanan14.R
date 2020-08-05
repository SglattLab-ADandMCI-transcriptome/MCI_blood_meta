## Import data from GSE33000, aka Narayanan14
## GCH

## TODO so negative in log2

require(data.table)
require(GEOquery)
require(limma)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Narayanan14"

cat("Reading files from GSE33000 aka Narayanan14\n")
rawcovs = getGEO(filename = paste0(rawlocation,"GSE33000/GSE33000_series_matrix.txt.gz"))
data = fread(paste0(rawlocation,"GSE33000/GSE33000_raw_data.txt.gz"),data.table =F)

Exprs = data[grep("Pool",names(data),invert=T)]
foo = strsplit(Exprs$Gene,",")
foo = lapply(foo,`[`,1)
foo = unlist(foo)
Exprs = Exprs[,-c(1:3)]
Exprs = normalizeQuantiles(Exprs)
# Exprs = log2(Exprs)
Exprs = asinh(Exprs)
Exprs = data.frame(PROBEID = foo, Exprs)

##all subjects are channel 2
rawcovs = rawcovs@phenoData@data

covs = data.frame(Sample_ID = rawcovs$geo_accession,
                  FACTOR_tissue = "PFC",
                  FACTOR_dx = rawcovs$`disease status:ch2`,
                  FACTOR_age = rawcovs$`age:ch2`,
                  FACTOR_sex = rawcovs$`gender:ch2`,
                  FACTOR_race = "caucasian")

covs$FACTOR_dx = factor(covs$FACTOR_dx)
levels(covs$FACTOR_dx) = c("AD","HD","CTL")
covs$FACTOR_age = sub(" yrs","",covs$FACTOR_age)

names(Exprs) = c("PROBEID",as.character(covs$Sample_ID))

head(covs)

rm(data)
rm(rawcovs)
