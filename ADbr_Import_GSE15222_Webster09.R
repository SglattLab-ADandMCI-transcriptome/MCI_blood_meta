## Import data from GSE15222 aka Webster09
## GCH

require(data.table)
require(GEOquery)
require(biomaRt)
require(limma)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Webster09"

cat("Reading text file from GSE15222 aka Webster09\n")
data = getGEO(filename = paste0(rawlocation,"GSE15222/GSE15222_series_matrix.txt.gz"))
rawcovs = fread(paste0(rawlocation,"GSE15222/samples.covar"))

Exprs = exprs(data)
## TODO i'm still not sure about this guy.
## These data already background subtracted
## asinh approximates ln(x) + a constant for "large x"
# Exprs = normalizeQuantiles(Exprs)
# Exprs = log2(Exprs)
Exprs = asinh(Exprs)  ## TODO final decision on arcsinh
Exprs = Exprs/log(2)
Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)

# Column 1    Group identifier
# Column 2    Individual identifier
# Column 3    Diagnosis (1=unaffected, 2=affected)
# Column 4    Age
# Column 5    APOE
# Column 6    Region (1=frontal, 2=parietal, 3=temporal, 4=cerebellar)
# Column 7    PMI
# Column 8    Site
# Column 9    Hybridization Date

covs = data.frame(Sample_ID = data@phenoData@data$geo_accession,
                  FACTOR_tissue = data@phenoData@data$source_name_ch1,
                  FACTOR_dx = data@phenoData@data$source_name_ch1,
                  FACTOR_age = data@phenoData@data$characteristics_ch1,
                  FACTOR_sex = data@phenoData@data$characteristics_ch1,
                  FACTOR_ethnicity = "Caucasian")

tissues = factor(levels = c(1:4))
for(i in 1:nrow(rawcovs)){
  ID = paste0(rawcovs[i,1],"-",rawcovs[i,2])
  index = which(data@phenoData@data$title == ID)
  tissues[index] = rawcovs[i,6]
}
levels(tissues) = c("FCX", "PCX", "TCX", "CER")
covs$FACTOR_tissue = tissues

## information from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2667989/
## all noncaucasians were excluded

print(covs$FACTOR_dx)
covs$FACTOR_dx[grep("Alzheimer",covs$FACTOR_dx)] = "AD"
covs$FACTOR_dx[grep("normal",covs$FACTOR_dx)] = "CTL"
print(covs$FACTOR_dx)
print(covs$FACTOR_age)
covs$FACTOR_age = sub("^.*,AGE_","",covs$FACTOR_age)
covs$FACTOR_age = sub(",PMI.*","",covs$FACTOR_age)
# Fixes age data by copying from another location
for(i in 1:length(covs$FACTOR_age)){
  if(covs$FACTOR_age[i] == "tissue: Alzheimer's Disease frozen cortical tissue") {
    foo = data@phenoData@data$characteristics_ch1.2[i]
    covs$FACTOR_age[i] = substr(foo,6,20)
  }
}
print(covs$FACTOR_age)
print(covs$FACTOR_sex)
covs$FACTOR_sex = sub("GENDER_","",covs$FACTOR_sex)
covs$FACTOR_sex = sub(",AGE.*","",covs$FACTOR_sex)
# Fixes gender data by copying from another location
for(i in 1:length(covs$FACTOR_sex)){
  if(covs$FACTOR_sex[i] == "tissue: Alzheimer's Disease frozen cortical tissue") {
    foo = data@phenoData@data$characteristics_ch1.1[i]
    foo = tolower(substr(foo, 9,14))
    covs$FACTOR_sex[i] = as.character(foo)
  }
}
print(covs$FACTOR_sex)
print(head(covs))

# Change probe IDs to more familiar symbols
IDs =  data@featureData@data$GB_ACC
IDs = sub("\\..*","",IDs)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
foo = getBM(attributes=c("refseq_mrna","hgnc_symbol"),filters="refseq_mrna",values=IDs,mart=ensembl)
genes = character()
for(i in 1:length(IDs)){
  thing = foo$hgnc_symbol[which(foo$refseq_mrna == IDs[i])]
  if(length(thing)>0) genes[i] = thing[1] else genes[i] = NA
}

Exprs$PROBEID = genes

rm(data)
rm(foo)
