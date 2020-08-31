## Import data from syn22024536 aka MCSA
## GCH

require(data.table)
require(edgeR)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MCSA"

data = fread(paste0(rawlocation,"blood/syn22024536/mcsa_rnaseq_rawcount.txt"), data.table=F)
rawcovs = fread(paste0(rawlocation,"blood/syn22024536/MCSA_individual_human_metadata.csv"), data.table=F)
rawages = fread(paste0(rawlocation,"blood/syn22024536/MCSA_biospecimen_metadata.csv"), data.table=F)

## this is reading in a line of NA at the bottom, 64254
genes = data$GeneName[-64254]

Exprs = data[-64254,-c(1:7)]
names = names(Exprs)
Exprs = cpm(Exprs,log=F)
Exprs = normalizeQuantiles(Exprs)

## TODO ask Chunling/Glatt about deleting all with all 0 == 16987

Exprs = asinh(Exprs)

Exprs = data.frame(PROBEID = genes, Exprs)
names(Exprs) = c("PROBEID",names)

covs = data.frame(Sample_ID = rawcovs$individualID,
                  FACTOR_tissue = "whole_blood",
                  FACTOR_dx = rawcovs$diagnosis,
                  FACTOR_age = NA,
                  FACTOR_sex = rawcovs$sex,
                  FACTOR_race = rawcovs$race)

ages = data.frame(Sample_ID = rawages$specimenID,
                  number = rawages$individualID,
                  age = rawages$samplingAge)

order = numeric()
for(i in 1:nrow(covs)){
  order[i] = which(ages$number == covs$Sample_ID[i])
}
ages = ages[order,]
all(ages$number == covs$Sample_ID)  ##TRUE

covs$Sample_ID = ages$Sample_ID
covs$FACTOR_age = ages$age

covs$FACTOR_race[covs$FACTOR_race == "Caucasian"] = "white"
covs$FACTOR_race = tolower(covs$FACTOR_race)

covs$FACTOR_dx[grep("Alzh",covs$FACTOR_dx)] = "AD"
covs$FACTOR_dx[grep("normal",covs$FACTOR_dx)] = "CTL"
covs$FACTOR_dx[grep("mild",covs$FACTOR_dx)] = "MCI"


head(covs)

rm(data)
rm(rawcovs)
rm(rawages)