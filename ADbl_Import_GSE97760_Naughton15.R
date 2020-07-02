## Import data from GSE97760 aka Naugton15
## GCH

## TODO age and ethnicity search

require(data.table)
require(plyr)
require(GEOquery)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Naugton15"

cat("Reading data files from GSE97760 aka Naugton15\n")

files = list.files(paste0(rawlocation,"blood/GSE97760/GSE97760_RAW"), full.names =  T, pattern = '.txt', recursive =T)
Exprs = list()
for(q in 1:length(files)){
  cat("\nReading",files[q])
  bar = fread(files[[q]], h=T,skip=9,data.table=F)
  Exprs[[q]] = bar$gMedianSignal
  name = gsub(".+/","",files[[q]])
  name = gsub(" \\(.+","",name)
  names(Exprs)[q] = name
}
Exprs = ldply(Exprs)

trans = normalizeQuantiles(t(Exprs[,-1]))
trans_norm = log2(trans)
IDs = Exprs$.id
Exprs = data.frame(PROBEID = bar$ProbeName,trans_norm)
IDs = gsub("_.*","",IDs)
names(Exprs) = c("PROBEID",IDs)

## probe to gene names
## Agilent-039494 SurePrint G3 Human GE v2 8x60K Microarray 039381
translist = fread(paste0(rawlocation,"blood/GSE97760/GPL16699-15607.txt"),skip=10,stringsAsFactors = F)
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


## covs TODO
gse = getGEO("GSE97760")
covs = data.frame(
    Sample_ID = gse$GSE97760_series_matrix.txt.gz@phenoData@data$geo_accession,
    FACTOR_dx = gse$GSE97760_series_matrix.txt.gz@phenoData@data$`disease:ch1`,
    FACTOR_sex = gse$GSE97760_series_matrix.txt.gz@phenoData@data$`gender:ch1`,
    FACTOR_age = NA,
    FACTOR_ethnicity = NA,
    FACTOR_tissue = gse$GSE97760_series_matrix.txt.gz@phenoData@data$`tissue:ch1`
)
## names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_ethnicity","FACTOR_tissue")
          
covs$FACTOR_sex = tolower(covs$FACTOR_sex)
covs$FACTOR_dx = gsub("healthy","CTL",covs$FACTOR_dx)
covs$FACTOR_dx = gsub("advanced Alzheimer's disease","AD",covs$FACTOR_dx)


rm(gse)
rm(trans)
rm(trans_norm)
rm(translist)
rm(bar)
rm(probes)
