## Import data from E-MTAB-6094 aka Leandro14
## GCH

require(data.table)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Leandro14"

cat("Reading data from E-MTAB-6094 aka Leandro14\n")

files = list.files(paste0(rawlocation,"/blood/E-MTAB-6094/"), full.names = T, pattern = "US09523754", recursive=T)
Exprs = list()
for(foo in 1:length(files)){
  cat("\nReading file:",files[[foo]])
  NONread = read.table(files[[foo]], skip=9, h=T, sep="\t", fill = T)
  subj = sub(".*_ *(.*?) *.txt.*", "\\1", files[[foo]])
  NONread = data.frame(NONread$GeneName, NONread$gMedianSignal)
  names(NONread)[2] = subj
  NONread = as.data.frame(t(NONread),stringsAsFactors = F)
  names(NONread) = NONread[1,]
  Exprs[[foo]] = NONread[2,]
  names(Exprs)[foo] = subj
}
foo = ldply(Exprs,as.numeric)
names(foo) = c("id",names(Exprs[[1]]))

norm = normalizeQuantiles(t(foo[,-1]))
# trans = log(norm,2)
trans = asinh(norm)

Exprs = data.frame(names(foo)[-1],trans, stringsAsFactors = F)
names(Exprs) = c("PROBEID",toupper(foo$id))


## covariates
covs = fread(paste0(rawlocation,"/blood/E-MTAB-6094/E-MTAB-6094.sdrf.txt"),data.table = F)
covs = fread(paste0(rawlocation,"blood/E-MTAB-6094/extractedcovs.txt"),data.table=F)
covs = covs[,c(1,1,3,2,7,1)]
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_dx = gsub("P.*","AD",covs$FACTOR_dx)
covs$FACTOR_dx = gsub("C.*","CTL",covs$FACTOR_dx)
## https://en.wikipedia.org/wiki/Pardo_Brazilians
covs$FACTOR_ethnicity = gsub("brown","other",covs$FACTOR_ethnicity)
covs$FACTOR_sex[covs$FACTOR_sex == "F"] = "female"
covs$FACTOR_sex[covs$FACTOR_sex == "M"] = "male"
covs$FACTOR_tissue = "whole blood"

covs = fread(paste0(rawlocation,"blood/E-MTAB-6094/extractedcovs.txt"),data.table=F)


rm(norm)
rm(trans)
rm(NONread)
rm(foo)
