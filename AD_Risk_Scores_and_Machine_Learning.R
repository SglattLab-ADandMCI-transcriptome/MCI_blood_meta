setwd("~/PsychGENe/brain/")

##risk scores and machine learning for AD brains

require(data.table)
require(limma)

scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")

tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

for (tissue in tissues){
  data = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_OutliersRemoved_allstudies.txt"),data.table=F)
  data = data[grep("CTL|AD",data$FACTOR_dx),]
  data$FACTOR_dx = factor(data$FACTOR_dx, levels = c("CTL","AD"))
  
  ## create risk score of odds ratios
  # Model matrix (basic)
  PredListNames = paste('FACTOR_', c("dx", "sex", "age","ethnicity","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  design = model.matrix( ~ ., predictors)
  y = t(data[,-grep("FACTOR_",names(data))])
  y = y[which(rowSums(is.na(y))<1),]
  lowvar = apply(y,1,var)
  fit = lmFit(y, design)
  dif = eBayes(fit)
  tops = topTable(dif)
  #lmFit, topTable
  # top differentially expressed genes
  # odds ratios
  # testing/ROC for numbers of genes (20,50,100)
  
  ## create svm or something model
  # leave one out
  # ROCs
  
  ## random forest also?
}
