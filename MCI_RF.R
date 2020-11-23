setwd("~/PsychGENe/MCI_blood_meta/")
## RF for MCI blood
set.seed(13210)

## TODO cuda and big grid ????
## todo make not svm

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
require(plyr)
require(foreach)
require(doParallel)

registerDoParallel(cores = 8)

scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")
tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

tissue = tissues[1]
for (tissue in tissues){
  data = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_training.txt"),data.table=F)
  
  data$FACTOR_age = as.numeric(sub("\\+","",data$FACTOR_age))
  
  data = data[grep(paste0(controllabel,"|",caselabel,"$"),data$FACTOR_dx),]
  study_id = unique(data$FACTOR_studyID)
  for(study in study_id){
    foo = sum(data$FACTOR_dx[data$FACTOR_studyID == study]==caselabel)
    if(foo==0){
      data = data[-which(data$FACTOR_studyID == study),]
    }
  }

  PredListNames = paste('FACTOR_', c("sampleID","dx", "sex", "age","race","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  
  x = data[,-grep("FACTOR_",names(data))]
  badpreds = which(colSums(is.na(predictors)) > 0)
  if(length(badpreds) > 0){
    predictors = predictors[,-badpreds]
  }
  badgenes = which(colSums(is.na(x)) > 1)
  if(length(badgenes) > 0){
    x = x[,-badgenes]
  }
  dx = factor(predictors$FACTOR_dx)
  cat("\nStudies: ",unique(predictors$FACTOR_studyID,"\n"))
  
  ## cell abundance PCs
  ##Generate abundances PCs
  cat("Generating abundances principal components\n")
  abfiles = list.files("./deconvolution","abundances", full.names = T)
  ablist = list()
  for(file in 1:length(abfiles)){
    ablist[[file]] = fread(abfiles[file], data.table=F)
  }
  abundances = ldply(ablist)
  abnames = abundances$V1
  abundances = abundances[,-1]
  row.names(abundances) = abnames
  abpc = prcomp(abundances)
  abpc3 = data.frame(abpc$x[,c(1:3)])
  names(abpc3) = c("cells PC1", "cells PC2", "cells PC3")
  
  ##include first 3 principal components of blood cell abundances
  abindex = numeric()
  subjects = predictors$FACTOR_sampleID
  for(j in 1:length(subjects)){
    abindex[j] = which(abnames == subjects[j])
  }
  print("all(subjects == abnames[abindex])")
  print(all(subjects == abnames[abindex]))
  abtemp = abpc3[abindex,]
  predictors = data.frame(predictors,abtemp)
  dx = factor(predictors$FACTOR_dx)

  ## residuals to take out known covariates
  design = model.matrix( ~ ., predictors[,-which(names(predictors) %in% c("FACTOR_dx","FACTOR_sampleID"))])
  y=t(x)
  resfit = lmFit(y, design)
  residual = resid(resfit,y)
  residual = data.frame(t(residual))
  genes = names(residual)
  
    
  
  
  ###TODO OR do we somewhere here let caret do the grid search and folding
  ## TODO we need some ranking and selection here?
  ## fsstats in exprso or RFE in caret
  ## default gamma  (default: 1/(data dimension))
  ## default cost = 1
  rfectl = rfeControl(functions = caretFuncs,
                      method = "repeatedcv",
                      # repeats = 5,
                      repeats = 2,
                      allowParallel = TRUE,
                      verbose = TRUE)

  profile = rfe(x = residual,
                y = dx,
                method = "svmRadial",
                sizes = c(10,20),
                rfeControl = rfectl)
  
  
  tctrl = trainControl(method = 'repeatedcv',
                       number = 5,
                       repeats = 5,
                       search = "grid",
                       allowParallel = TRUE,
                       # savePredictions = "all",
                       verboseIter = TRUE)
  
  grid = expand.grid(C = c(1:2), ## maybe 10^-3:3?
                     sigma = c(1)*(1/ncol(residual))) ## maybe .1 to 10?
  
  machines = train(x = residual,
  # machines = train(x = subset, ## TODO put back to residual or RFE stuff
                   y = dx,
                   method = "svmRadial",
                   # classProbs = TRUE,
                   trControl = tctrl,
                   tuneGrid = grid,
                   verboseIter = TRUE)
  
  
  
  # ###### SVM
  # ## create folds
  # foldindex = sample(1:numfolds,size = nrow(residual),replace=T)
  # 
  # sink(file = paste("./classifier_results/",tissue,"_SVM_cmatrix_",analysislabel,".txt", sep = ""))
  # sink()
  # 
  # machines = list()
  # results = list()
  # eachdx = list()
  # cmatrices = list()
  # for(fold in 1:numfolds){
  #   cat("\nPerforming SVM for fold:",fold,"\n")
  #   ##split training and test
  #   trainindex = which(foldindex != fold)
  #   train = residual[trainindex,]
  #   test = residual[-trainindex,]
  #   traindx = dx[trainindex]
  #   testdx = dx[-trainindex]
  #   eachdx[[fold]] = testdx
  #   
  #   machine = svm(traindx ~.,data = train,
  #                 probability = T, scale = F)
  #   machines[[fold]] = machine
  #   
  #   ## run on training data
  #   trainresult = predict(machine, train, decision.values = T, probability = T)
  #   
  #   ## test
  #   result = predict(machine, test, decision.values = T, probability = T)
  #   results[[fold]] = result
  #   
  #   cmatrix = confusionMatrix(result,testdx,positive="MCI")
  #   cmatrices[[fold]] = cmatrix 
  #   print(cmatrix)
  #   sink(file = paste("./classifier_results/",tissue,"_SVM_cmatrix_",analysislabel,".txt", sep = ""), append = TRUE)
  #   cat("SVM fold",fold)
  #   print(cmatrix)
  #   sink()
  # 
  #   pr = prcomp(x)
  #   pca = data.frame(pr$x)
  #   png(file = paste("./QCplots/",tissue,"_SVM_PCAtest_",analysislabel,"_fold_",fold,".png", sep = ""),
  #       res = 300, units = "in", height = 8, width = 8)
  #   plot(pca$PC1 ~ pca$PC2)
  #   dev.off()
  #   
  #   ## ROC for test and train for this fold
  #   probs = attributes(result)$probabilities
  #   probs = probs[,1]
  #   png(file = paste("./QCplots/",tissue,"_SVM_",analysislabel,"_test_ROC_fold",fold,".png", sep = ""),
  #       res = 300, units = "in", height = 8, width = 8)
  #   riskroc = roc(response = testdx, predictor = probs)
  #   plot(riskroc,main=paste0("SVM test set ROC, fold ",fold))
  #   text(.8,0,labels = paste0("AUC = ",riskroc$auc))
  #   dev.off()
  # 
  #   trainprobs = attributes(trainresult)$probabilities
  #   trainprobs = trainprobs[,1]
  #   png(file = paste("./QCplots/",tissue,"_SVM_",analysislabel,"_train_ROC_fold",fold,".png", sep = ""),
  #       res = 300, units = "in", height = 8, width = 8)
  #   riskroc = roc(response = traindx, predictor = trainprobs)
  #   plot(riskroc,main=paste0("SVM training set ROC, fold ",fold))
  #   text(.8,0,labels = paste0("AUC = ",riskroc$auc))
  #   dev.off()
  # }
  # ##save fold info
  # ## get average auc test/train
  # ## get max sens/spec
  # ## draw test combined roc curves and write average auc on it
  # ## draw train combined roc curves and write average auc on it
  # ##threshold for max sens/spec
  
  trainpreds = predict.train(machines)
  confusionMatrix(as.factor(dx),trainpreds)
  trainprobs = predict.train(machines, type="prob")

  
  ## TODO do final training with best parameters
  ## roc/auc for full training set
  
  ##TODO compare to validation set
  ## run on the validation set
  ## make roc curve with average auc on it
  ## max sens/spec
}