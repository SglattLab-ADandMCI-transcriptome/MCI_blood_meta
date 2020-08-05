## Import data from syn5550404 aka MayoRNAseq
## GCH

## we had decided 90 or above would be shifted to 90

require(data.table)
require(limma)
require(biomaRt)
require(edgeR)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MayoRNAseq"

CBEdata = fread(paste0(rawlocation,"syn5550404/MayoRNAseq_RNAseq_CBE_transcriptCounts.tsv"), data.table=F)
TCXdata = fread(paste0(rawlocation,"syn5550404/MayoRNAseq_RNAseq_TCX_transcriptCounts.tsv"), data.table=F)

CBEcovs = fread(paste0(rawlocation,"syn5550404/MayoRNAseq_RNAseq_CBE_covariates.csv"), data.table=F)
TCXcovs = fread(paste0(rawlocation,"syn5550404/MayoRNAseq_RNAseq_TCX_covariates.csv"), data.table=F)

genes = TCXdata$ensembl_id
CBEnames = names(CBEdata)[-1]
TCXnames = names(TCXdata)[-1]

CBEExprs = cpm(CBEdata[-1],log=F)
TCXExprs = cpm(TCXdata[-1],log=F)
CBEExprs = normalizeQuantiles(CBEExprs)
TCXExprs = normalizeQuantiles(TCXExprs)
CBEExprs = asinh(CBEExprs)
TCXExprs = asinh(TCXExprs)

Exprs = cbind(CBEExprs,TCXExprs)
Exprs = data.frame(PROBEID = genes, Exprs)
names(Exprs) = c("PROBEID",CBEnames,TCXnames)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
foo = getBM(attributes=c("ensembl_transcript_id","hgnc_symbol"),filters="ensembl_transcript_id",values=Exprs$PROBEID,mart=ensembl)
genes = character()
for(i in 1:length(Exprs$PROBEID)){
  thing = foo$hgnc_symbol[which(foo$ensembl_transcript_id == Exprs$PROBEID[i])]
  if(length(thing)>0) genes[i] = thing[1] else genes[i] = NA
}

Exprs$PROBEID = genes

covs = data.frame(Sample_ID = c(CBEcovs$SampleID,TCXcovs$ID),
                  FACTOR_tissue = c(rep("CBE",nrow(CBEcovs)),rep("TCX",nrow(TCXcovs))),
                  FACTOR_dx = c(CBEcovs$Diagnosis,TCXcovs$Diagnosis),
                  FACTOR_age = c(CBEcovs$AgeAtDeath,TCXcovs$AgeAtDeath),
                  FACTOR_sex = c(CBEcovs$Sex,TCXcovs$Gender),
                  FACTOR_race = "unknown")

covs$FACTOR_dx = factor(covs$FACTOR_dx)
levels(covs$FACTOR_dx) = c("AD","CTL","Pathologic Aging","PSP")
levels(covs$FACTOR_sex) = c("female","male")
foo = as.character(covs$FACTOR_age)
foo[which(foo == "90_or_above")] = "90"
covs$FACTOR_age = foo
head(covs)

rm(CBEdata)
rm(CBEcovs)
rm(TCXdata)
rm(TCXcovs)
