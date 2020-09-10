## Import data from Shigemizu20
## GCH

require(data.table)
require(openxlsx)
require(biomaRt)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Shigemizu20"

cat("Reading normalized data from Shigemizu20 aka Shigemizu20\n")

data = fread(paste0(rawlocation,"blood/Shigemizu20/cpm.txt"), data.table=F)
rawcovs = read.xlsx(paste0(rawlocation,"blood/Shigemizu20/13195_2020_654_MOESM1_ESM.xlsx"))


names(rawcovs) = rawcovs[1,]
rawcovs = rawcovs[-1,]

covs = data.frame(Sample_ID = rawcovs$ID,
                  FACTOR_tissue = "whole_blood",
                  FACTOR_dx = rawcovs$phenotype,
                  FACTOR_age = rawcovs$age,
                  FACTOR_sex = rawcovs$sex,
                  FACTOR_race = "asian")  ## as per email from Dr. Shigemizu

covs$FACTOR_dx[covs$FACTOR_dx == "CN"] = "CTL"
covs$FACTOR_sex[covs$FACTOR_sex == "M"] = "male"
covs$FACTOR_sex[covs$FACTOR_sex == "F"] = "female"

##add ID prefix
newID = gsub("^","S20_",covs$Sample_ID)
covs$Sample_ID = newID


newID = gsub("^","S20_",names(data)[-1])
genes = data$gene
Exprs = data[,-1]
allzero = rowSums(Exprs)
allzero = which(allzero == 0)
Exprs = Exprs[-allzero,]
genes = genes[-allzero]
Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)
names(Exprs) = newID
Exprs = data.frame(PROBEID = genes, Exprs) 


## convert to hgnc symbols
cat("\nPerforming HGNC symbol conversion.")
IDs =  Exprs$PROBEID
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genetable = getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), values=IDs,mart=ensembl)
genes = character()
for(i in 1:length(IDs)){
  ID = strsplit(IDs[i],"\\.")[[1]][1]
  thing = genetable$hgnc_symbol[which(genetable$ensembl_gene_id == ID)]
  if(length(thing)>0) genes[i] = thing[1] else genes[i] = NA
}
cat("\n",sum(is.na(genes)),"of",length(genes),"genes retired or otherwise missing.\n")
Exprs$PROBEID=genes


rm(data)
rm(rawcovs)