## Import data from syn22024498 aka ROSMAP
## GCH

require(data.table)
require(edgeR)
require(sva)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "ROSMAP"

cat("Reading raw count data from syn22024498 aka ROSMAP\n")

data = fread(paste0(rawlocation,"blood/syn22024498/S598_ExpMatrix_9904genes.txt"), data.table=F)
rawcovs = fread(paste0(rawlocation,"blood/syn22024498/ROSMAP_clinical.csv"), data.table=F)
rawages = fread(paste0(rawlocation,"blood/syn22024498/ROSMAP_biospecimen_metadata.csv"), data.table=F)
binfo = fread(paste0(rawlocation, "blood/syn22024498/S615_BatchInfo.txt"))

## this is reading in a line of NA at the bottom, 64254
genes = data$ID
names = names(data)
Exprs = data[,-1]  ## CPM wants rows as genes
Exprs = cpm(Exprs,log=F)
Exprs = data.frame(Exprs)

cat("Filtering genes with less than", filterCPM, "CPM in",
    filterPercent*100, "percent or more of subjects\n")
filtered = logical()
for(i in 1:nrow(Exprs)){
  filtered[i] = FALSE
  quux = Exprs[i,]<filterCPM
  if(sum(quux) >= filterPercent*ncol(Exprs)){
    filtered[i] = TRUE
  }
}
filtered = which(filtered)
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
batch = factor(binfo$bInfo[])

Exprs = ComBat(Exprs,batch)

Exprs = normalizeQuantiles(Exprs)
Exprs = asinh(Exprs)

Exprs = data.frame(PROBEID = genes, Exprs)

Exprs = Exprs[,-which(names(Exprs) == "Sample_615")] ## this is a control sample


covs = data.frame(Sample_ID = rawcovs$individualID,
                  FACTOR_tissue = "whole_blood",
                  FACTOR_dx = rawcovs$dcfdx_lv,
                  FACTOR_age = NA,
                  FACTOR_sex = rawcovs$msex,
                  FACTOR_race = rawcovs$race)

rawages = rawages[which(rawages$tissue == "blood"),]

ages = data.frame(Sample_ID = rawages$specimenID,
                  number = rawages$individualID,
                  age = rawages$samplingAge)

ages$age = gsub("\\+","", ages$age)

foo = which(ages$Sample_ID %in% names(Exprs))  ##get sampleIDs/indivIDs we have
ages = ages[foo,]  #trim ages

all(ages$number %in% rawcovs$individualID)
which(!(ages$number %in% covs$Sample_ID))

toremove = numeric()
for(name in unique(ages$number)){
  # name = unique(ages$number)[10]
  things = which(ages$number == name)
  if(length(things) > 1){
    foo = which(ages$age[things] == max(ages$age[things]))
    cat("max",ages$age[things]," -------- ",ages$age[things[foo]],"\n")
    bad = things[-foo]
    toremove = c(toremove,bad)
  }
}
ages = ages[-toremove,]

foo = which(covs$Sample_ID %in% ages$number)
covs = covs[foo,]

lineup = numeric()
for(i in 1:length(covs$Sample_ID)){
  lineup[i] = which(ages$number == covs$Sample_ID[i])
}
ages = ages[lineup,]

covs$FACTOR_age = as.numeric(ages$age)

covs$FACTOR_sex[which(covs$FACTOR_sex == 1)] = "male"
covs$FACTOR_sex[which(covs$FACTOR_sex == 0)] = "female"

# 1 White
# 2 Black or African American
# 3 American Indian or Alaska Native
# 4 Native Hawaiian or Other Pacific Islander
# 5 Asian
# 6 Other
# 7 Unknown
covs$FACTOR_race[which(covs$FACTOR_race == 1)] = "white"
covs$FACTOR_race[which(covs$FACTOR_race == 2)] = "black"
covs$FACTOR_race[which(covs$FACTOR_race == 3)] = "nativeamerican"
covs$FACTOR_race[which(covs$FACTOR_race == 4)] = "pacificislander"
covs$FACTOR_race[which(covs$FACTOR_race == 5)] = "asian"
covs$FACTOR_race[which(covs$FACTOR_race == 6)] = "other"
covs$FACTOR_race[which(covs$FACTOR_race == 7)] = "unknown"

# 1 NCI: No cognitive impairment (No impaired domains)
# 2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
# 3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
# 4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)
# 5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
# 6 Other dementia: Other primary cause of dementia
covs$FACTOR_dx[which(covs$FACTOR_dx == 1)] = "CTL"
covs$FACTOR_dx[which(covs$FACTOR_dx == 2)] = "MCI"
covs$FACTOR_dx[which(covs$FACTOR_dx == 3)] = "MCIplus"
covs$FACTOR_dx[which(covs$FACTOR_dx == 4)] = "AD"
covs$FACTOR_dx[which(covs$FACTOR_dx == 5)] = "ADplus"
covs$FACTOR_dx[which(covs$FACTOR_dx == 6)] = "Other"

summary(covs$FACTOR_age)
table(covs$FACTOR_race)
table(covs$FACTOR_dx)
table(covs$FACTOR_sex)
table(covs$FACTOR_sex,covs$FACTOR_dx)

head(covs)

samplestosave = numeric()
foo = names(Exprs)
for(i in 1:length(ages$Sample_ID)){
  samplestosave[i] = which(foo == ages$Sample_ID[i])
}
all(foo[samplestosave] == ages$Sample_ID)
genes = Exprs$PROBEID
Exprs = Exprs[,c(1,samplestosave)]
names(Exprs) = c("PROBEID",ages$number)


rm(data)
rm(rawcovs)
rm(rawages)