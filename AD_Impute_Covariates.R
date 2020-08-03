## impute missing covariates, but only if we have 70% of them
## GCH

require(smcfcs)
require(data.table)


datafolder = "./data_for_analysis/"

covfiles = list.files(datafolder,"_SampleFactors_allstudies.txt")

covtis = sub("_SampleFactors_allstudies.txt","",covfiles)

tissues = covtis

for(tissue in tissues){
  message("Ethnicities: ",tissue)
  datExprCovs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"), data.table=F, stringsAsFactors = F)
  print(table(datExprCovs$FACTOR_ethnicity))
  for(study in unique(datExprCovs$FACTOR_studyID)){
    message(study)
    print(table(datExprCovs$FACTOR_ethnicity[datExprCovs$FACTOR_studyID == study]))
  }
}



# datExprCovs$FACTOR_ethnicity[grep("unknown|other",datExprCovs$FACTOR_ethnicity)] = NA