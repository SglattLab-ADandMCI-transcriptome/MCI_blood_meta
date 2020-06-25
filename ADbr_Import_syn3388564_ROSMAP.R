## Import of syn3388564 aka ROSMAP RNAseq
## GCH

require(data.table)
require(dplyr)
require(biomaRt)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "ROSMAP_RNAseq"

cat("\nReading files from syn3388564 aka ROSMAP RNAseq")

# This has already been normalized, so only requires log2 and combination/formatting
# They also used ComBat
data1 = fread(paste0(rawlocation,"syn3388564/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv"),data.table=F)
data2 = fread(paste0(rawlocation,"syn3388564/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv"),data.table=F)
genes = data1$gene_id
all(genes == data1$gene_id)
all(genes == data2$gene_id)

data1 = data1[,-c(1,2)]
data2 = data2[,-c(1,2)]

## add 1 so log doesn't produce -inf
data1 = log2(data1+1)
data2 = log2(data2+1)

Exprs = data.frame(PROBEID = genes,data1,data2)

##Sequencing samples that were sequenced more than once include the sequencing
##batch number as the last 2 digits of the id. For example, sample 123_456789_02
##is sample "123_456789", from batch 2. The sample IDs are otherwise randomly-
##assigned identifiers. Some samples were re-sequenced due to poor quality or
##low output, they may have "redo" in the sample identifier.

foo = strsplit(names(Exprs)[-1],"_")
foo = data.frame(foo)
bar = unlist(foo[1,])
baz = which(duplicated(bar))
baz = bar[baz]
baz = unique(baz)
bax = bar %in% baz
foo[,bax]

names(Exprs)[grep("redo",names(Exprs))]
grep("4_140501",names(Exprs))

## The above tests show we just need to remove a few duplicates.
foo = foo[-3,]
specimens = character()
for(i in 1:ncol(foo)){
  specimens[i] = paste0(foo[1,i],"_",foo[2,i])
}
specimens = gsub("X","",specimens)
specimens = unique(specimens)

remove = which(duplicated(bar,fromLast=T)) +1
Exprs = Exprs[,-remove]
cat("\n",length(remove),"duplicate/redundant runs removed")

## Convert to individual IDs and attach covariates
rawkey = fread(paste0(rawlocation,"syn3388564/ROSMAP_biospecimen_metadata.csv"),data.table=F)
rawcovs = fread(paste0(rawlocation,"syn3388564/ROSMAP_Clinical_2019-05_v3.csv"),data.table=F)

all(specimens %in% rawkey$specimenID)
key = data.frame(specimen = specimens,individual = NA)
for(specimen in key$specimen){
  key$individual[which(key$specimen == specimen)] = rawkey$individualID[grep(specimen,rawkey$specimenID)]
}
newnames = c("PROBEID",key$individual)
names(Exprs) = newnames

## two specimens are missing individual IDS
remove = which(newnames == "")
Exprs = Exprs[,-remove]
cat("\n",length(remove),"entries removed for missing Individual IDs")

covs = data.frame(Sample_ID = names(Exprs)[-1],
                  FACTOR_tissue = "PFC",
                  FACTOR_dx = NA,
                  FACTOR_age = NA,
                  FACTOR_sex = NA,
                  FACTOR_ethnicity = NA)

old = 0
for(i in 1:nrow(covs)){
  age = rawcovs$age_death[which(rawcovs$individualID == covs$Sample_ID[i])]
  if(age == "90+"){
    age = 90
    old = old+1
  }
  covs$FACTOR_age[i] = age
}
cat("\n",old,"ages listed as 90+ converted to 90")

for(i in 1:nrow(covs)){
  sex = rawcovs$msex[which(rawcovs$individualID == covs$Sample_ID[i])]
  if(sex){
    sex = "male"
  } else {
    sex = "female"
  }
  covs$FACTOR_sex[i] = sex
}

# Value Coding
# 1 NCI: No cognitive impairment (No impaired domains)
# 2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
# 3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
# 4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)
# 5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
# 6 Other dementia: Other primary cause of dementia

for(i in 1:nrow(covs)){
  dx = rawcovs$dcfdx_lv[which(rawcovs$individualID == covs$Sample_ID[i])]
  covs$FACTOR_dx[i] = case_when(
    dx == 1 ~ "CTL",
    dx == 2 ~ "MCI",
    dx == 3 ~ "MCI+",
    dx == 4 ~ "AD",
    dx == 5 ~ "AD+",
    dx == 6 ~ "Other",
  )
}
sum(covs$FACTOR_dx == "CTL")
sum(covs$FACTOR_dx == "MCI")
sum(covs$FACTOR_dx == "MCI+")
sum(covs$FACTOR_dx == "AD")
sum(covs$FACTOR_dx == "AD+")
sum(covs$FACTOR_dx == "Other")

# Value Coding
# 1 White
# 2 Black or African American
# 3 American Indian or Alaska Native
# 4 Native Hawaiian or Other Pacific Islander
# 5 Asian
# 6 Other
# 7 Unknown
# Are you of Spanish/Hispanic/Latino origin?
#   value coding
# 1 Yes
# 2 No
esp = 0
wesp = 0
for(i in 1:nrow(covs)){
  ethnicity = rawcovs$race[which(rawcovs$individualID == covs$Sample_ID[i])]
  hispanic = rawcovs$spanish[which(rawcovs$individualID == covs$Sample_ID[i])]
  covs$FACTOR_ethnicity[i] = case_when(
    hispanic == 1 ~ "hispanic",
    ethnicity == 1 ~ "caucasian",
    ethnicity == 2 ~ "black",
    ethnicity == 3 ~ "nativeamerican",
    ethnicity == 4 ~ "islander",
    ethnicity == 5 ~ "asian",
    ethnicity == 6 ~ "other",
    TRUE ~ "unknown"
  )
  if(hispanic == 1){
    esp = esp+1
    if(ethnicity==1){
      wesp = wesp+1
    }
  }
}
cat("\n",esp,"total hispanic and",wesp,"white hispanic (all coded as hispanic)\n")

print(head(covs))

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
cat("\n",sum(is.na(genes)),"of",length(genes),"genes retired or otherwise missing")
Exprs$PROBEID=genes

rm(genes)
rm(data1)
rm(data2)
rm(rawkey)
rm(rawcovs)
rm(key)
rm(foo)
rm(ensembl)
rm(genetable)
rm(IDs)