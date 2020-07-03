## Import data from syn4552659 aka MSBB_ATP
## GCH

require(data.table)
require(gcrma)
require(biomaRt)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MSBB_ATP"

rawcovs = fread(paste0(rawlocation,"syn4552659/AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv"))
keys = fread(paste0(rawlocation,"syn4552659/AMP-AD_MSBB_MSSM_AffymetrixArray_sample_key.txt"))

nums = substr(keys$individualIdentifier,17,30)
keys$Sample_ID = paste0("MSBB",nums,keys$REGION)

for(i in 1:length(keys$Sample_ID)){
  cat(i,"\r")
  keys$FACTOR_sex[i] = rawcovs$Sex[which(rawcovs$individualIdentifier == keys$individualIdentifier[i])]
  keys$FACTOR_age[i] = rawcovs$Age[which(rawcovs$individualIdentifier == keys$individualIdentifier[i])]
  keys$FACTOR_ethnicity[i] = rawcovs$Race[which(rawcovs$individualIdentifier == keys$individualIdentifier[i])]
  keys$FACTOR_dx[i] = rawcovs$CDR[which(rawcovs$individualIdentifier == keys$individualIdentifier[i])]
}

covs = data.frame(Sample_ID = keys$fileName,
                  FACTOR_tissue = keys$REGION,
                  FACTOR_dx = keys$FACTOR_dx,
                  FACTOR_age = keys$FACTOR_age,
                  FACTOR_sex = keys$FACTOR_sex,
                  FACTOR_ethnicity = keys$FACTOR_ethnicity)

print(covs$FACTOR_dx)
covs$FACTOR_dx[covs$FACTOR_dx >= 1] = "AD"
covs$FACTOR_dx[covs$FACTOR_dx < 1] = "CTL"
print(covs$FACTOR_dx)

print(covs$FACTOR_ethnicity)
covs$FACTOR_ethnicity = factor(covs$FACTOR_ethnicity)
levels(covs$FACTOR_ethnicity) = c("asian","black","hispanic","white")
print(covs$FACTOR_ethnicity)
print(covs$FACTOR_sex)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
levels(covs$FACTOR_sex) = c("female","male")
print(covs$FACTOR_sex)
print(covs$FACTOR_tissue)
covs$FACTOR_tissue = factor(covs$FACTOR_tissue)
levels(covs$FACTOR_tissue) = c("other","other","other","TCX",
                               "other","other","other","other",
                               "other","other","other","PFC",
                               "other","other","other","other",
                               "HIP","other","other")
print(covs$FACTOR_tissue)
print(head(covs))

cat("Reading CEL files from syn4552659 aka MSBB_ATP\n")
afiles = paste0(grep("U133A",covs$Sample_ID,value=T))
bfiles = paste0(grep("U133B",covs$Sample_ID,value=T))
## combo chips
## RMA is similar to log2 and normalize = T will quantile normalize
aRMA = justGCRMA(filenames = afiles,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"syn4552659/"))
aExprs = exprs(aRMA)
aExprs = data.frame(PROBEID = rownames(aExprs), aExprs)
## right now the 49804 b file is corrupted, which kills the whole thing. excluding.
## relevant to a small portion of batch effect correction.
bfiles = bfiles[-grep("AMP-AD_MSBB_MSSM_AffymetrixHG-U133B_49804hg133b11.cel",bfiles,ignore.case = T)]
bRMA = justGCRMA(filenames = bfiles,normalize = T,optimize.by = "memory",celfile.path = paste0(rawlocation,"syn4552659"))
bExprs = exprs(bRMA)
bExprs = data.frame(PROBEID = rownames(bExprs), bExprs)

anames = sub("^.*U133A_","MSBB",names(aExprs))
anames = sub("hg133.*","",anames)
bnames = sub("^.*U133B_","MSBB",names(bExprs))
bnames = sub("hg133.*","",bnames)
lonea = which(!(anames %in% bnames))
loneb = which(!(bnames %in% anames))
aExprs = aExprs[,-lonea]
bExprs = bExprs[,-loneb]
anames = anames[-lonea]
bnames = bnames[-loneb]
names(aExprs) = anames
names(bExprs) = bnames
all((anames %in% bnames))
all((bnames %in% anames))
Exprs = rbind(aExprs,bExprs)


IDs =  Exprs$PROBEID
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
## listattributes for affymetrix a and b
foo0 = getBM(attributes=c("affy_hg_u133a","hgnc_symbol"), filters="affy_hg_u133a", values=IDs,mart=ensembl)
foo1 = getBM(attributes=c("affy_hg_u133a_2","hgnc_symbol"), filters="affy_hg_u133a_2", values=IDs,mart=ensembl)
foo2 = getBM(attributes=c("affy_hg_u133b","hgnc_symbol"), filters="affy_hg_u133b", values=IDs,mart=ensembl)
names(foo0) = c("probe","hgnc_symbol")
names(foo1) = c("probe","hgnc_symbol")
names(foo2) = c("probe","hgnc_symbol")
foo = rbind(foo0,foo1,foo2)
genes = character()
for(i in 1:length(IDs)){
  thing = foo$hgnc_symbol[which(foo$probe == as.character(IDs[i]))]
  if(length(thing)>0) genes[i] = thing[1] else genes[i] = NA
}
Exprs$PROBEID=genes

covs$Sample_ID = sub("^.*U133A_","MSBB",covs$Sample_ID)
covs$Sample_ID = sub("^.*U133B_","MSBB",covs$Sample_ID)
covs$Sample_ID = sub("hg133.*","",covs$Sample_ID)

rm(aRMA)
rm(bRMA)
rm(rawcovs)
rm(aExprs)
rm(bExprs)
