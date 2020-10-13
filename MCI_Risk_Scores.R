setwd("~/PsychGENe/MCI_blood_meta/")
##risk scores for MCI blood

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
# require(metafor)

scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")
tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

tissue = tissues[1]
for (tissue in tissues){
  ## FACTOR_etc in the first few columns, genes in the rest, subjects as rows
  data = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_training.txt"),data.table=F)
  
  ## TODO PCs of blood regress

  data$FACTOR_age = as.numeric(sub("\\+","",data$FACTOR_age))
  
  data = data[grep(paste0(controllabel,"|",caselabel,"$"),data$FACTOR_dx),]
  study_id = unique(data$FACTOR_studyID)
  for(study in study_id){
    foo = sum(data$FACTOR_dx[data$FACTOR_studyID == study]==caselabel)
    if(foo==0){
      data = data[-which(data$FACTOR_studyID == study),]
    }
  }
  ## create risk score of odds ratios
  # Model matrix (basic)
  PredListNames = paste('FACTOR_', c("dx","sex", "age","race","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  design = model.matrix( ~ ., predictors)
  y = data.frame(t(data[,-grep("FACTOR_",names(data))]))
  names(y) = data$FACTOR_sampleID
  y = y[,-which(rowSums(is.na(predictors)) > 0)]
  predictors = predictors[-which(rowSums(is.na(predictors)) > 0),]
  y = y[which(rowSums(is.na(y))<1),]
  lowvar = apply(y,1,var)
  fit = lmFit(y, design)
  dif = eBayes(fit)
  tops = topTable(dif)
  
  cat("\nStudies: ",unique(data$FACTOR_studyID),"\n")
  
  ## residuals to take out known covariates
  ## TODO can't use SVs but can use cell abundance?
  design = model.matrix( ~ ., predictors[,-which(names(predictors)=="FACTOR_dx")])
  resfit = lmFit(y, design)
  residual = data.frame(resid(resfit,y))
  names(residual) = names(y)
  
  # residual = y ##for testing only TODO delete this
  
  
  ## k-fold PTRS
  # load ptrs
  source("~/PsychGENe/ptrs/ptrs.R")
  
  ## create folds
  foldindex = sample(1:numfolds,size = ncol(residual),replace=T)
  genes = rownames(residual)
  write(foldindex, file = paste("./classifier_results/",tissue,"_ptrs_indices_",analysislabel,".txt", sep = ""))
  
  fold = 1
  allcoefs = list()
  fits = list()
  for(fold in 1:numfolds){
    ## construct test and train tables
    train = residual[,which(foldindex != fold)]
    test = residual[,which(foldindex == fold)]
    minipredictors = predictors[which(foldindex != fold),]
    testdx = predictors$FACTOR_dx[which(foldindex == fold)]
    train = as.data.frame(t(train))
    
    
    ##miniature mega-analysis

    ## linear model
    cat("\nFitting linear model for fold",fold,"\n")
    minidx = factor(minipredictors$FACTOR_dx)
    miniFit = lm(as.matrix(train) ~ -1 + minidx)
    miniSummary = summary(miniFit)
    miniCoefs = lapply(miniSummary, function(x) x$coefficients)
    miniCoefs = ldply(miniCoefs)
    miniCoefs = miniCoefs[c(1:(nrow(miniCoefs)/2))*2,]
    miniCoefs$.id = gsub("Response ","",miniCoefs$.id)
    weights = data.frame(GeneSymbol = miniCoefs$.id,
                         P = miniCoefs$`Pr(>|t|)`,
                         score = miniCoefs$Estimate)
    
    # ## logit model
    # cat("\nFitting logit model for fold",fold,"\n")
    # minidx = ifelse(minipredictors$FACTOR_dx == "CTL", 0, 1)
    # miniFit = glm(minidx ~ as.matrix(train), 
    #                   family=binomial(link='logit'),
    #                   control = list(maxit = 100))
    # miniSummary = summary(miniFit)
    # miniCoefs = data.frame(miniSummary$coefficients)
    # rownames(miniCoefs) = gsub("as\\.matrix\\(train\\)","", rownames(miniCoefs))
    # miniCoefs = miniCoefs[-which(rownames(miniCoefs) == "(Intercept)"),]
    # weights = data.frame(GeneSymbol = rownames(miniCoefs),
    #                      P = miniCoefs$Pr...z..,
    #                      score = miniCoefs$Estimate)


    fits[[fold]] = miniFit
    fwrite(weights,paste0("./classifier_results/",tissue,"_",analysislabel,"_",fold,"_weights.txt"),sep="\t")

    # load weights
    weights = load_weights(file=paste0("./classifier_results/",tissue,"_",analysislabel,"_",fold,"_weights.txt"),
                           gene_col = "GeneSymbol",
                           weight_col = "score",
                           p_col = "P")
    
    datExpr = data.frame(t(data[-grep("FACTOR_",names(data))]))
    genes = names(data)[-grep("FACTOR_",names(data))]
    names(datExpr) = data$FACTOR_sampleID
    row.names(datExpr) = genes
    
    # run ptrs on data
    p_thres = c(0.9, 0.5, 0.1, 0.05, 0.01, 0.005, 0.002, 0.001)
    score_df = ptrs(dat = test, weight_table = weights, p_thres = p_thres)
    
    badscore = (colSums(is.na(score_df))>0)
    score_df = score_df[,!badscore,drop=F]
    p_thres = p_thres[!badscore]
    
    trainscore_df = ptrs(dat = t(train), weight_table = weights, p_thres = p_thres)
    
    # plot score density
    for(i in 1:ncol(score_df)){
      png(file = paste("./QCplots/",tissue,"_ptrs_",analysislabel,"_density_fold_",fold,"_",i,".png", sep = ""),
          res = 300, units = "in", height = 8, width = 8)
      plot(density(score_df[,i]),main = paste0("Density for risk scores, fold ",fold,", p < ",p_thres[i]))
      dev.off()
    }
  
    # contrast cases and controls on ptrs
    # ## linear model
    # fit = lm(as.matrix(score_df) ~ testdx)
    # coefs = summary(fit)
    # coefs = lapply(coefs, function(x) x$coefficients)
    # coefs = ldply(coefs)
    # coefs = coefs[c(1:(nrow(coefs)/2))*2,]
    # coefs = data.frame(p_thres,coefs)
    # print(coefs)
    
    ## logit model
    cat("\nFitting logit model on ptrs for fold",fold,"\n")
    testdxb = ifelse(testdx == "CTL", 0, 1)
    riskFit = glm(testdxb ~ as.matrix(score_df),
                      family=binomial(link='logit'),
                      control = list(maxit = 100))
    riskSummary = summary(riskFit)
    riskCoefs = data.frame(riskSummary$coefficients)
    print(riskCoefs)

    ## ROC for test and train
    for(i in 1:ncol(score_df)){
      png(file = paste("./QCplots/",tissue,"_ptrs_",analysislabel,"_test_ROC_fold",fold,"_",i,".png", sep = ""),
          res = 300, units = "in", height = 8, width = 8)
      riskroc = roc(response = testdx, predictor = score_df[,i])
      plot(riskroc,main=paste0("Risk score test set ROC, fold ",fold,", p < ",p_thres[i]))
      text(.8,0,labels = paste0("AUC = ",riskroc$auc))
      dev.off()
    }
    
    for(i in 1:ncol(score_df)){
      png(file = paste("./QCplots/",tissue,"_ptrs_",analysislabel,"_train_ROC_fold",fold,"_",i,".png", sep = ""),
          res = 300, units = "in", height = 8, width = 8)
      riskroc = roc(response = minidx, predictor = trainscore_df[,i])
      plot(riskroc,main=paste0("Risk score training set ROC, fold ",fold,", p < ",p_thres[i]))
      text(.8,0,labels = paste0("AUC = ",riskroc$auc))
      dev.off()
    }
  }
  ## TODO save fits and compare
  ## save AUCs
  ## save coefs
  
  
  ## TODO compare to validation set
  
  
  
  #lmFit, topTable
  # top differentially expressed genes
  # odds ratios
  # testing/ROC for numbers of genes (20,50,100)
  
  # leave one out
  # ROCs
  
  ## random forest also?
  
  ## neural net
}
