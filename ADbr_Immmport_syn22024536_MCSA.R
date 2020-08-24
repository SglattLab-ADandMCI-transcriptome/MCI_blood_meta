## Import data from syn22024536 aka MCSA
## GCH

stop("TODO")
## TODO this is not at all done or ready

require(data.table)
require(limma)
require(biomaRt)
require(edgeR)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MCSA"

data = fread(paste0(rawlocation,"syn22024536/mcsa_rnaseq_rawcount.txt"), data.table=F)

covs = fread(paste0(rawlocation,"syn22024536/MCSA_individual_human_metadata.csv"), data.table=F)

genes = data$ensembl_id
CBEnames = names(CBEdata)[-1]


Exprs = cpm(CBEdata[-1],log=F)
Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)

Exprs = data.frame(PROBEID = genes, Exprs)
names(Exprs) = c("PROBEID",CBEnames)

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

rm(data)
rm(covs)
