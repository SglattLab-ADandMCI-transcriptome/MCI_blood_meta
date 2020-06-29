## Import data from GSE4229 aka Maes09
## GCH

## TODO some ages known, but not all. decide things.

require(data.table)
require(readtext)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Maes09"

cat("Reading data from GSE4229 aka Maes09\n")

files = list.files(paste0(rawlocation,"/blood/GSE4229"), full.names = T, pattern = "sample_table", recursive=T)
foo = fread(files[1])
baz = as.data.frame(matrix(nrow=nrow(foo),ncol=0))
for(thing in files){
  foo = fread(thing)
  bar = foo[,2,drop=F]
  colnames(bar)=thing
  baz = cbind(baz,bar)
}
rownames(baz)=rownames(foo)
foo = unlist(lapply(strsplit(files,"/|_"), `[[`, 4))
colnames(baz)=foo
baz = normalizeQuantiles(baz)
baz = baz+1 ##log(0) is frown
baz = log(baz,2)
Exprs = data.frame(PROBEID = rownames(baz), baz)

subjects = gsub("_sample_table.txt","",files)
subjects = gsub(paste0(rawlocation,"/blood/GSE4229/data/"),"",subjects)
names(Exprs) = c("PROBEID",subjects)

## probe/gene translation
## NIA MGC, Mammalian Genome Collection
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL1211&id=7513&db=GeoDb_blob92
translist = fread(paste0(rawlocation,"/blood/GSE4229/data/array/GPL1211_GEO_Table.html"),skip=25,nrows=9600)
Exprs$PROBEID = translist$`GENE_SYMBOL</strong><strong>`


## covariates
covs = fread(paste0(rawlocation,"/blood/GSE4229/E-GEOD-4229.sdrf.txt"),data.table=F)
covs = covs[,c(1,6,6,1,1,7)] #age, ethnicity missing
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_age = NA
covs$FACTOR_ethnicity = NA
covs$FACTOR_dx[grepl("control",covs$FACTOR_dx)] = "CTL"
covs$FACTOR_dx[grepl("Alzheimer's",covs$FACTOR_dx)] = "AD"
covs$FACTOR_sex[grepl("Male",covs$FACTOR_sex)] = "male"
covs$FACTOR_sex[grepl("Female",covs$FACTOR_sex)] = "female"
covs$FACTOR_tissue[grepl("mononuclear cell",covs$FACTOR_tissue)] = "whole blood"
covs$Sample_ID = gsub(" 1","",covs$Sample_ID)

ages = readtext(paste0(rawlocation,"/blood/GSE4229/1-s2.0-S0197458006002983-mmc2.doc"))
## from pmid 16979800
ages = ages$text
ages = unlist(strsplit(ages,"\\|"))
ages = ages[-grep("^\\s+$",ages)]
ages = data.frame(matrix(ages[-c(1:11,110)],ncol=7,byrow=T), stringsAsFactors = F)
ages$dx = c(rep("NEC",14),rep("AD",14))
names(ages) = c("Subject","Gender","Age","Education","MMSE","ApoE","Total RNA","dx")
ages$Subject = paste(ages$dx,trimws(ages$Subject),trimws(ages$Gender),sep="")
rownames(ages) = ages$Subject
ages$Age = as.numeric(ages$Age)

## rename from gse to IDs
bar = fread(paste0(rawlocation,"/blood/GSE4229/maestitles.txt"), header=F,data.table=F)
names(bar) = c("FACTOR_sampleID","FACTOR_ID")
subjects = names(Exprs)[-1]
for(i in 1:length(subjects)) {
  subjects[i] = bar$FACTOR_ID[bar$FACTOR_sampleID == subjects[i]]
}
names(Exprs) = c("PROBEID",subjects)
subjects = covs$Sample_ID
for(i in 1:length(subjects)) {
  subjects[i] = bar$FACTOR_ID[bar$FACTOR_sampleID == subjects[i]]
}
covs$Sample_ID = subjects

## and finally attach the ages
for(i in ages$Subject){
  covs$FACTOR_age[covs$Sample_ID == i] = ages$Age[ages$Subject == i]
}

