## Import data from syn22024536 aka MCSA
## GCH

require(data.table)
require(edgeR)
require(sva)
require(biomaRt)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "MCSA"

cat("Reading raw count data from syn22024536 aka MCSA\n")

# data = fread(paste0(rawlocation,"blood/syn22024536/mcsa_rnaseq_rawcount.txt"), data.table=F)
data = fread(paste0(rawlocation,"blood/syn22024536/MCSA_Data_20201210/S422_RawCounts_60715genes.txt"), data.table=F)
rawcovs = fread(paste0(rawlocation,"blood/syn22024536/MCSA_individual_human_metadata.csv"), data.table=F)
rawages = fread(paste0(rawlocation,"blood/syn22024536/MCSA_biospecimen_metadata.csv"), data.table=F)
binfo = fread(paste0(rawlocation, "blood/syn22024536/MCSA_Data_20201210/S422_BatchInfo.txt"))

## original file has line of NA at the bottom, 64254
# genes = data$GeneName[-64254]
# Exprs = data[-64254,-c(1:7)]
genes = data$ID
Exprs = data[,-1]
names = names(Exprs)
Exprs = cpm(Exprs,log=F)
Exprs = data.frame(Exprs)

cat("Filtering genes with less than", filterCPM, "CPM in",
    filterPercent*100, "percent or more of subjects\n")
filtered = logical()
for(i in 1:nrow(Exprs)){
  cat("\r",i,"of",nrow(Exprs),"genes")
  filtered[i] = FALSE
  quux = Exprs[i,]<filterCPM
  if(sum(quux) >= filterPercent*ncol(Exprs)){
    filtered[i] = TRUE
  }
}
cat(sum(filtered),"genes filtered.\n")
filtered = which(filtered)
foo = sum(Exprs[filtered,])
bar = sum(Exprs)
cat(foo,"total CPM filtered, or",foo/bar,"of all")
if(length(filtered) > 0){
  Exprs = Exprs[-filtered,]
  genes = genes[-filtered]
}


##COMBAT
cat("COMBAT for known batch effects.\n")
print("all(names(Exprs) %in% binfo$sName)")
print(all(names(Exprs) %in% binfo$sName))
foo = which(binfo$sName %in% names(Exprs))
binfo = binfo[foo,]
all(names(Exprs) == binfo$sName)
batch = factor(binfo$run)

Exprs = ComBat(Exprs,batch)


Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)

Exprs = data.frame(PROBEID = genes, Exprs)
names = gsub("S","",names)
names = gsub("$","_PAX",names)
names(Exprs) = c("PROBEID",names)


## probeid from ensg to gene names
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
bm = getBM(attributes=c('ensembl_gene_id_version', 'hgnc_symbol'), 
      filters = 'ensembl_gene_id_version', 
      values = Exprs$PROBEID, 
      mart = mart)
Exprs$PROBEID[which(!Exprs$PROBEID %in% bm$ensembl_gene_id_version)]
newnames = Exprs$PROBEID
for(q in 1:length(newnames)){
  temp = bm$hgnc_symbol[which(bm$ensembl_gene_id_version == newnames[q])]
  if(length(temp) >= 1 && temp != ""){
    newnames[q] = temp[1]
  }
}
Exprs$PROBEID = newnames


covs = data.frame(Sample_ID = rawcovs$individualID,
                  FACTOR_tissue = "whole_blood",
                  FACTOR_dx = rawcovs$diagnosis,
                  FACTOR_age = NA,
                  FACTOR_sex = rawcovs$sex,
                  FACTOR_race = rawcovs$race)

ages = data.frame(Sample_ID = rawages$specimenID,
                  number = rawages$individualID,
                  age = rawages$samplingAge)

ages$age = gsub("\\+","", ages$age)

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

