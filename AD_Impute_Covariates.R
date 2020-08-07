## impute missing covariates, but only if we have 70% of them
## GCH

## TODO all of this, probably only for a few ethnicities

require(smcfcs)
require(data.table)


datafolder = "./data_for_analysis/"

covfiles = list.files(datafolder,"_SampleFactors_allstudies.txt")

covtis = sub("_SampleFactors_allstudies.txt","",covfiles)

tissues = covtis

for(tissue in tissues){
  message("Ethnicities: ",tissue)
  datExprCovs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"), data.table=F, stringsAsFactors = F)
  print(table(datExprCovs$FACTOR_race))
  for(study in unique(datExprCovs$FACTOR_studyID)){
    message(study)
    print(table(datExprCovs$FACTOR_race[datExprCovs$FACTOR_studyID == study]))
  }
}



# datExprCovs$FACTOR_race[grep("unknown|other",datExprCovs$FACTOR_race)] = NA
