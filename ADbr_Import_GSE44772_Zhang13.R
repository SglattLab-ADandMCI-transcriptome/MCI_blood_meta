## Import data from GSE44772, aka Zhang13
## GCH

require(data.table)
require(GEOquery)
require(limma)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Zhang13"

cat("Reading files from GSE44772 aka Zhang13\n")
rawcovs = getGEO(filename = paste0(rawlocation,"GSE44772/GSE44772_series_matrix.txt.gz"))
dataC = fread(paste0(rawlocation,"GSE44772/GSE44768_originalDataFile_CR.txt.gz"),data.table =F)
dataF = fread(paste0(rawlocation,"GSE44772/GSE44770_originalDataFile_PFC.txt.gz"),data.table =F)

foo = strsplit(dataC$gene,",")
foo = lapply(foo,`[`,1)
foo = unlist(foo)

ExprsC = dataC[grep("Pool",names(dataC),invert=T)]
ExprsC = ExprsC[,-c(1:3)]
ExprsC = normalizeQuantiles(ExprsC)
ExprsC = log2(ExprsC)

ExprsF = dataF[grep("Pool",names(dataF),invert=T)]
ExprsF = ExprsF[,-c(1:3)]
ExprsF = normalizeQuantiles(ExprsF)
ExprsF = log2(ExprsF)

## TODO are the pfc subjects the same as from Narayanan14?
# Exprs = data.frame(PROBEID = foo, ExprsC, ExprsF)
Exprs = data.frame(PROBEID = foo, ExprsC)

##all subjects are channel 2
rawcovs = rawcovs@phenoData@data

covs = data.frame(Sample_ID = rawcovs$geo_accession,
                  FACTOR_tissue = rawcovs$source_name_ch2,
                  FACTOR_dx = rawcovs$characteristics_ch2,
                  FACTOR_age = rawcovs$characteristics_ch2.1,
                  FACTOR_sex = rawcovs$characteristics_ch2.2,
                  FACTOR_ethnicity = "caucasian")

covs$FACTOR_dx = factor(covs$FACTOR_dx)
print(covs$FACTOR_dx)
covs$FACTOR_dx = factor(covs$FACTOR_dx)
levels(covs$FACTOR_dx) = c("AD","CTL")
print(covs$FACTOR_dx)
foo = strsplit(as.character(covs$FACTOR_tissue),"_")
foo = lapply(foo,`[`,2)
foo = unlist(foo)
covs$FACTOR_tissue = foo
covs$FACTOR_age = sub("age: ","",covs$FACTOR_age)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
print(covs$FACTOR_sex)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
levels(covs$FACTOR_sex) = c("female","male")
print(covs$FACTOR_sex)

## rename exprs to gse
foo = sub("X","",names(Exprs))
for(i in 2:length(foo)){
  foo[i] = rawcovs$geo_accession[which(rawcovs$title == foo[i])]
}
names(Exprs) = foo

foo = which(covs$Sample_ID %in% names(Exprs))
covs = covs[foo,]

print(head(covs))

rm(dataC)
rm(dataF)
rm(ExprsC)
rm(ExprsF)
rm(rawcovs)

