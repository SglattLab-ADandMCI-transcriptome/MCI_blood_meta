setwd("~/PsychGENe/MCI_blood_meta/")
## SVM for MCI blood

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")

numfolds = 5

# make the new folder
if(!dir.exists("./classifier_results/")){ dir.create("./classifier_results/")}

require(data.table)
require(limma)
require(pROC)
require(e1071)
require(caret)

scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")

tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

tissue = tissues[1]
for (tissue in tissues){
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

  PredListNames = paste('FACTOR_', c("dx", "sex", "age","race","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  
  ## SVM
  ## TODO cuda and big grid
  x = data[,-grep("FACTOR_",names(data))]
  x = x[-which(rowSums(is.na(predictors)) > 0),]
  predictors = predictors[-which(rowSums(is.na(predictors)) > 0),]
  x = x[,which(colSums(is.na(x))<1)]  ##TODO jon re missingness tolerance
  dx = factor(predictors$FACTOR_dx)

  cat("\nStudies: ",unique(predictors$FACTOR_studyID,"\n"))
  
  ## residuals to take out known covariates
  ## TODO can't use SVs but can use cell abundance?
  design = model.matrix( ~ ., predictors[,-which(names(predictors)=="FACTOR_dx")])
  y=t(x)
  resfit = lmFit(y, design)
  residual = resid(resfit,y)
  residual = data.frame(t(residual))
  
    
  ## create folds
  foldindex = sample(1:numfolds,size = nrow(residual),replace=T)

  genes = rownames(residual)
  
  sink(file = paste("./classifier_results/",tissue,"_SVM_cmatrix_",analysislabel,".txt", sep = ""))
  sink()
  
  machines = list()
  results = list()
  eachdx = list()
  cmatrices = list()
  for(fold in 1:numfolds){
    cat("\nPerforming SVM for fold:",fold,"\n")
    ##split training and test
    trainindex = which(foldindex != fold)
    train = residual[trainindex,]
    test = residual[-trainindex,]
    traindx = dx[trainindex]
    testdx = dx[-trainindex]
    eachdx[[fold]] = testdx
    
    machine = svm(traindx ~.,data = train,probability = T)
    machines[[fold]] = machine
    
    ## test
    result = predict(machine, test, decision.values = T, probability = T)
    results[[fold]] = result
    
    cmatrix = confusionMatrix(result,testdx,positive="MCI")
    cmatrices[[fold]] = cmatrix 
    print(cmatrix)
    sink(file = paste("./classifier_results/",tissue,"_SVM_cmatrix_",analysislabel,".txt", sep = ""), append = TRUE)
    cat("SVM fold",fold)
    print(cmatrix)
    sink()
  
    pr = prcomp(x)
    pca = data.frame(pr$x)
    png(file = paste("./QCplots/",tissue,"_SVM_PCAtest_",analysislabel,"_fold_",fold,".png", sep = ""),
        res = 300, units = "in", height = 8, width = 8)
    plot(pca$PC1 ~ pca$PC2)
    dev.off()
    
    ## ROC I suppose
    probs = attributes(result)$probabilities
    probs = probs[,1]
    riskroc = roc(response = testdx, predictor = probs)
    png(file = paste("./QCplots/",tissue,"_SVM_",analysislabel,"_ROC_fold_",fold,".png", sep = ""),
        res = 300, units = "in", height = 8, width = 8)
    plot(riskroc,main=paste0("SVM probabilities ROC, test set fold ",fold))
    text(.8,0,labels = paste0("AUC = ",riskroc$auc))
    dev.off()
  }
  ##save fold info and compare
  ##threshold for max sens/spec
}