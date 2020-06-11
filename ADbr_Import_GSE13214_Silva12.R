## Import data from GSE13214 aka Silva12
## GCH

require(data.table)
require(GEOquery)
require(limma)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Silva12"

cat("Reading SOFT file from GSE13214 aka Silva12\n")
data = getGEO(filename = paste0(rawlocation,"GSE13214/GSE13214_family.soft.gz"))

chans = integer()
rawcovs = character()
for(i in 1:length(data@gsms)){
  test = c(data@gsms[[i]]@header$characteristics_ch1, data@gsms[[i]]@header$characteristics_ch2)
  chans[i] = grep("braak",test)
  rawcovs[i] = paste(data@gsms[[i]]@header$title,test[chans[i]])
}
Exprs = list()
for(i in 1:length(data@gsms)){
  chan = c(chans[i]*2,chans[i]*2+1)
  name = names(data@gsms)[i]
  foo = data@gsms[[i]]@dataTable@table[,chan]
  bar = foo[,1]-foo[,2]
  bar[which(bar <0)] = 0
  Exprs[[i]] = bar
  names(Exprs)[i] = name
}
Exprs = ldply(Exprs)
Exprs = Exprs[,-1]
Exprs = t(Exprs)
Exprs = normalizeQuantiles(Exprs)
Exprs = log2(Exprs)
Exprs = data.frame(PROBEID = data@gpls$GPL1930@dataTable@table$GenName,Exprs)
names(Exprs) = c("PROBEID",names(data@gsms))

covs = data.frame(Sample_ID = names(data@gsms),
                  FACTOR_tissue = "unknown",
                  FACTOR_dx = "other",
                  FACTOR_age = "unknown",
                  FACTOR_sex = "unknown",
                  FACTOR_ethnicity = "unknown",
                  stringsAsFactors = F)

covs$FACTOR_tissue[grep("hippocampus",rawcovs)] = "HIP"
covs$FACTOR_tissue[grep("frontal_| cortex",rawcovs)] = "FCX"
covs$FACTOR_dx[grep("Alzheimer",rawcovs)] = "AD"
covs$FACTOR_dx[grep("normal_brain",rawcovs)] = "CTL"
covs$FACTOR_age = sub("^.*Age: ","",rawcovs)
covs$FACTOR_age = sub(" years","",covs$FACTOR_age)
covs$FACTOR_sex[grep("male",rawcovs)] = "male"
covs$FACTOR_sex[grep("female",rawcovs)] = "female"
print(head(covs))

## This dataset actually has double replicates for each sample
n = nrow(covs) /2
n = c(1:n) *2+1
## per Glatt, average them
## will do geometric because it's easier
Exprs1 = Exprs[,-c(1,n)]
Exprs2 = Exprs[,n]
foo = Exprs$PROBEID
Exprs = (Exprs1 + Exprs2) /2
Exprs = data.frame(PROBEID = foo, Exprs)
n = n-1
covs = covs[-n,]

rm(data)
rm(rawcovs)
rm(Exprs1)
rm(Exprs2)