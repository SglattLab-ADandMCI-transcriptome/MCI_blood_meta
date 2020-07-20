## Import data from GSE85426 aka Samsudin16
## GCH

## TODO ethnicity question

require(data.table)
require(GEOquery)
require(plyr)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Samsudin16"

cat("Reading normalized data from GSE85426 aka Samsudin16\n")

## recall that 2 subjects have been removed due to missing covs
files = list.files(paste0(rawlocation,"blood/GSE85426/raw microarray"), full.names =  T, pattern = '.txt', recursive =T)
Exprs = list()
for(q in 1:length(files)){
  bar = fread(files[[q]], h=T,skip=9,data.table=F)
  Exprs[[q]] = bar$gMedianSignal
  name = gsub(".+/","",files[[q]])
  name = gsub(" \\(.+","",name)
  names(Exprs)[q] = name
}
Exprs = ldply(Exprs)

trans = normalizeQuantiles(t(Exprs[,-1]))
# trans_norm = log2(trans)
trans_norm = asinh(trans)
IDs = Exprs$.id
Exprs = data.frame(PROBEID = bar$ProbeName,trans_norm)
names(Exprs) = c("PROBEID",IDs)

## probe to gene names
## Agilent-028004 SurePrint G3 Human GE 8x60K Microarray
translist = fread(paste0(rawlocation,"blood/GSE85426/GPL14550-9757.txt"),skip=10,stringsAsFactors = F)
q = length(Exprs$PROBEID)
cat("Replacing probe IDs with gene symbols.\n")
translist = data.frame(translist$GENE_SYMBOL,translist$ID,stringsAsFactors = F)
names(translist) = c("Symbol","ID")
probes = Exprs$PROBEID
probes = unlist(lapply(strsplit(probes,",|\\|"), `[[`, 1))
for(foo in 1:q){
  if(probes[foo] %in% translist$ID){
    bar = translist$Symbol[probes[foo] == translist$ID]
    bar = strsplit(bar," ")[[1]][1]
    if(!is.na(bar)){  ## I don't know why it pulls NA out of that, but here we go.
      probes[foo] = bar
    }
  }
  cat(foo,"/",q,"\r")
}
Exprs$PROBEID = probes


## covariates
## note 2 samples are lost because we do not have the raw form
covs = fread(paste0(rawlocation,"blood/GSE85426/E-GEOD-85426.sdrf.txt"),data.table=F)
covs = covs[,c(1,33,13,5,1,10)] #ethnicity missing
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
covs$FACTOR_ethnicity = NA
covs$FACTOR_dx[grepl("control",covs$FACTOR_dx)] = "CTL"
covs$FACTOR_dx[grepl("Alzheimer",covs$FACTOR_dx)] = "AD"
covs$Sample_ID = gsub(" 1","",covs$Sample_ID)

## These are extracted from the SOFT file available on GSE85426
foo = names(Exprs)
foo = unlist(lapply(strsplit(foo,".t|\\s|\\("), `[[`, 1))
names(Exprs) = foo

gse = getGEO(filename = paste0(rawlocation,"blood/GSE85426/GSE85426_family.soft.gz"))
foo = GSMList(gse)
foo = lapply(foo, Meta)
foo = unlist(foo)
foo = foo[grepl("description2",names(foo))]

translist = data.frame(foo,names(foo))
names(translist) = c("name","id")
translist$id = gsub(".description2","",translist$id)
IDs = names(Exprs)
IDs = gsub("(\\(| |\\.).*","",IDs)
thing = character()
for(i in 1:length(IDs)){
  thing = translist$id[translist$name == IDs[i]]
  if(length(thing)>0){
    IDs[i] = thing
  }
}
names(Exprs) = IDs

covs = covs[covs$Sample_ID %in% IDs,]

rm(gse)
rm(trans_norm)
rm(trans)
rm(translist)
