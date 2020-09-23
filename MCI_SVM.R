setwd("~/PsychGENe/brain/")
## SVM for MCI blood

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")

require(data.table)
require(pROC)
# require(dataPreparation)
require(e1071)
# require(ggplot2)
require(caret)

scalefiles = list.files("./data_for_analysis/","_GeneExpression_allstudies.txt")

tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)
tissues = "whole_blood"

tissue = tissues[1]
for (tissue in tissues){
  # expr = fread(paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"),data.table=F)
  # covs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"),data.table=F)
  # data = data.frame(covs,expr)
  # cat("all(data$FACTOR_sampleID == data$FACTOR_sampleID.1)")
  # print(all(data$FACTOR_sampleID == data$FACTOR_sampleID.1))
  # data = data[,-which(names(data) == "FACTOR_sampleID.1")]
 
  data = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_OutliersRemoved_allstudies.txt"),data.table=F)
  
  data$FACTOR_age = as.numeric(sub("\\+","",data$FACTOR_age))
  
  data = data[grep(paste0(controllabel,"|",caselabel,"$"),data$FACTOR_dx),]
  study_id = unique(data$FACTOR_studyID)
  for(study in study_id){
    foo = sum(data$FACTOR_dx[data$FACTOR_studyID == study]==caselabel)
    if(foo==0){
      data = data[-which(data$FACTOR_studyID == study),]
    }
  }
  cat("\nStudies: ",unique(data$FACTOR_studyID))
  
  PredListNames = paste('FACTOR_', c("dx", "sex", "age","ethnicity","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  
  
  ## SVM
  ## TODO don't need to scale?
  ## TODO cuda and big grid
  # foo = paste(names(predictors),collapse = " + ")
  # form = as.formula(paste0("~ ", foo))
  x = data[,-grep("FACTOR_",names(data))]
  x = x[,which(colSums(is.na(x))<1)]
  x = data.frame(scale(x))
  dx = factor(predictors$FACTOR_dx)
  
  ##split training and test
  trainindex = sample(1:nrow(x), nrow(x)*.7)
  train = x[trainindex,]
  test = x[-trainindex,]
  traindx = dx[trainindex]
  testdx = dx[-trainindex]
  
  machine = svm(traindx ~.,data = train)
  
  ## test
  result = predict(machine, test)
  
  cmatrix = confusionMatrix(result,testdx,positive="MCI")
  print(cmatrix)
  sink(file = paste("./classifier_results/",tissue,"_SVM_cmatrix_",analysislabel,".txt", sep = ""))
  print(cmatrix)
  sink()
  
  ## kfold cross validation?
  
  pr = prcomp(x)
  pca = data.frame(pr$x)
  png(file = paste("./QCplots/",tissue,"_PCAtest_",analysislabel,".png", sep = ""),
      res = 300, units = "in", height = 8, width = 8)
  plot(pca$PC1 ~ pca$PC2)
  dev.off()
}