## validation test split-off
## for use in final validation, not even k-fold stuff

set.seed(13210)

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"

require(data.table)

data = fread("./data_for_analysis/whole_blood_ScaledWithFactors_OutliersRemoved_allstudies.txt", data.table = F)

data$FACTOR_age = as.numeric(sub("\\+","",data$FACTOR_age))

data = data[grep(paste0(controllabel,"|",caselabel,"$"),data$FACTOR_dx),]
study_id = unique(data$FACTOR_studyID)
for(study in study_id){
  foo = sum(data$FACTOR_dx[data$FACTOR_studyID == study]==caselabel)
  if(foo==0){
    data = data[-which(data$FACTOR_studyID == study),]
  }
}


## create split index
## TODO do this smarter
index = sample(1:5,size = nrow(data),replace=T)  ##temporary random assignment strategy


## split studes and write files
validation = data[which(index == 1),]
training = data[which(index != 1),]

fwrite(validation,"./data_for_analysis/whole_blood_ScaledWithFactors_validation.txt",row.names=F)
fwrite(training,"./data_for_analysis/whole_blood_ScaledWithFactors_training.txt",row.names=F)

sink("./classifier_results/set_split_tables.txt")

cat("\nTraining set:")
print(table(training$FACTOR_dx))
print(table(training$FACTOR_studyID))
print(table(training$FACTOR_dx,training$FACTOR_studyID))

cat("\nValidation set:")
print(table(validation$FACTOR_dx))
print(table(validation$FACTOR_studyID))
print(table(validation$FACTOR_dx,validation$FACTOR_studyID))

sink()
