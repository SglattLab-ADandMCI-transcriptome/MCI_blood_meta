## Raw probe level mRNA expression data were exported
## from GenomeStudio software (Illumina Inc.) for preprocessing
## with background correction, variance stabilizing transformation,
## quantile normalization and probe filtering using the lumi package of BioConductor.
## https://www.synapse.org/#!Synapse:syn17062532

##do not forget we used ComBat

## For import of syn17062532 aka MayoeGWAS

require(data.table)
require(sva)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MayoeGWAS"

## TODO take out subjects who are in the pilot RNAseq study? no

CBEdata = fread(paste0(rawlocation,"syn17062532/MayoEGWAS_arrayExpression_CBE.csv"), data.table=F)
TCXdata = fread(paste0(rawlocation,"syn17062532/MayoEGWAS_arrayExpression_TCX.csv"), data.table=F)

CBEdata$IID = paste0(CBEdata$IID,"C")
TCXdata$IID = paste0(TCXdata$IID,"T")

CBEcovs = fread(paste0(rawlocation,"syn17062532/MayoEGWAS_arrayExpression_CBE_covariates.csv"), data.table=F)
TCXcovs = fread(paste0(rawlocation,"syn17062532/MayoEGWAS_arrayExpression_TCX_covariates.csv"), data.table=F)

CBEbatch = numeric()
CBEbatch[which(CBEcovs$plate0==1)] = 0
CBEbatch[which(CBEcovs$plate1==1)] = 1
CBEbatch[which(CBEcovs$plate2==1)] = 2
CBEbatch[which(CBEcovs$plate3==1)] = 3
CBEbatch[which(CBEcovs$plate4==1)] = 4
CBEbatch[which(is.na(CBEbatch))] = 6
TCXbatch = numeric()
TCXbatch[which(TCXcovs$plate2==1)] = 2
TCXbatch[which(TCXcovs$plate3==1)] = 3
TCXbatch[which(TCXcovs$plate4==1)] = 4
TCXbatch[which(TCXcovs$plate5==1)] = 5
TCXbatch[which(is.na(TCXbatch))] = 6

##run ComBat
prb = names(CBEdata)[-c(1,2)]
subj = c(CBEdata$IID,TCXdata$IID)
CBEdata = ComBat(t(CBEdata[,-c(1,2)]),CBEbatch)
TCXdata = ComBat(t(TCXdata[,-c(1,2)]),TCXbatch)

Exprs = rbind(t(CBEdata),t(TCXdata))
# subj = Exprs$IID
# prb = names(Exprs)[-c(1,2)]
# Exprs = Exprs[,-c(1,2)]
Exprs = t(Exprs)
Exprs = log2(Exprs)
Exprs = data.frame(PROBEID = prb, Exprs)
names(Exprs) = c("PROBEID",subj)

refs = fread("references/GPL8432-11703.txt",data.table=F,skip=29)
rnames = character()
for(i in 1:length(Exprs$PROBEID)){
  rnames[i] = refs$ILMN_Gene[which(refs$ID == Exprs$PROBEID[i])]
}
Exprs$PROBEID = rnames

CBEcovs$IID = paste0(CBEcovs$IID,"C")
TCXcovs$IID = paste0(TCXcovs$IID,"T")

covs = data.frame(Sample_ID = c(CBEcovs$IID,TCXcovs$IID),
                  FACTOR_tissue = c(rep("CBE",nrow(CBEcovs)),rep("TCX",nrow(TCXcovs))),
                  FACTOR_dx = c(CBEcovs$Dxn,TCXcovs$Dxn),
                  FACTOR_age = c(CBEcovs$Age,TCXcovs$Age),
                  FACTOR_sex = c(CBEcovs$Sex,TCXcovs$Sex),
                  FACTOR_ethnicity = "unknown")

## Diagnosis (AD=1, non-AD = 0)
covs$FACTOR_dx[covs$FACTOR_dx == 0] = "CTL"
covs$FACTOR_dx[covs$FACTOR_dx == 1] = "AD"
## female = 1 (350 total)
covs$FACTOR_sex[covs$FACTOR_sex == 0] = "male"
covs$FACTOR_sex[covs$FACTOR_sex == 1] = "female"
head(covs)



rm(CBEdata)
rm(CBEcovs)
rm(TCXdata)
rm(TCXcovs)
rm(refs)
