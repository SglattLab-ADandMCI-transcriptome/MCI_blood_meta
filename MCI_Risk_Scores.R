setwd("~/PsychGENe/MCI_blood_meta/")
##risk scores for MCI blood

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")

# make the new folder
if(!dir.exists("./classifier_results/")){ dir.create("./classifier_results/")}

require(data.table)
require(limma)
require(pROC)
require(e1071)

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
  

  ## TODO account for study and other covs?
  ## resid function on a fit w/o dx predictor
  ## output = resid(lm(GeneExprMatrix ~ Age + Sex, data = FactorData))
  ## create risk score of odds ratios
  # Model matrix (basic)
  PredListNames = paste('FACTOR_', c("dx", "sex", "age","ethnicity","studyID"), collapse= "|", sep="")
  predictors = data[,grep(PredListNames, colnames(data))]
  predictors = data.frame(predictors)
  design = model.matrix( ~ ., predictors)
  y = t(data[,-grep("FACTOR_",names(data))])
  y = y[,-which(rowSums(is.na(predictors)) > 0)]
  y = y[which(rowSums(is.na(y))<1),]
  lowvar = apply(y,1,var)
  fit = lmFit(y, design)
  dif = eBayes(fit)
  tops = topTable(dif)
  
  
  ## TODO use lmfit to get out the other predictors?
  
  
  ## create weight file
  sigs = fread(paste0("./meta_analysis/",tissue,"_",analysislabel,"_meta_significant_arcsinh_and_pvals.csv"))
  weights = data.frame(GeneSymbol = sigs$GeneSymbol,
                       P = sigs$P,
                       arcsinh = sigs$arcsinh)
  fwrite(weights,paste0("./meta_analysis/",tissue,"_",analysislabel,"_meta_weights.txt"),sep="\t")
  
  
  # load ptrs
  source("~/PsychGENe/ptrs/ptrs.R")
  
  # load weights
  weights = load_weights(file=paste0("./meta_analysis/",tissue,"_",analysislabel,"_meta_weights.txt"),
                         gene_col = "GeneSymbol",
                         weight_col = "arcsinh",
                         p_col = "P")
  
  datExpr = data.frame(t(data[-grep("FACTOR_",names(data))]))
  genes = names(data)[-grep("FACTOR_",names(data))]
  names(datExpr) = data$FACTOR_sampleID
  row.names(datExpr) = genes
  
  # run ptrs on data
  p_thres = c(0.9, 0.5, 0.1, 0.05, 0.01, 0.005, 0.002, 0.001)
  score_df = ptrs(dat = datExpr, weight_table = weights, p_thres = p_thres)
  
  badscore = (colSums(is.na(score_df))>0)
  score_df = score_df[,!badscore,drop=F]
  p_thres = p_thres[!badscore]
  
  # plot score density
  for(i in 1:ncol(score_df)){
    png(file = paste("./QCplots/",tissue,"_ptrs_",analysislabel,"_density_",i,".png", sep = ""),
        res = 300, units = "in", height = 8, width = 8)
    plot(density(score_df[,i]),main = paste0("Density for risk scores, p < ",p_thres[i])) # column named ptrs_1 is equal to p_thres = 0.5
    dev.off()
    }

  # contrast cases and controls on ptrs
  fit = lm(as.matrix(score_df) ~ data$FACTOR_dx)
  coefs = summary(fit)
  # coefs = lapply(coefs, function(x) broom::tidy(x$coefficients))
  coefs = lapply(coefs, function(x) x$coefficients)
  coefs = ldply(coefs)
  # coefs = coefs[grepl("dx", coefs$.rownames), ]
  coefs = coefs[c(1:(nrow(coefs)/2))*2,]
  coefs = data.frame(p_thres,coefs)
  print(coefs)
  
  ## ROC I suppose
  for(i in 1:ncol(score_df)){
    png(file = paste("./QCplots/",tissue,"_ptrs_",analysislabel,"_ROC_",i,".png", sep = ""),
        res = 300, units = "in", height = 8, width = 8)
    riskroc = roc(response = predictors$FACTOR_dx, predictor = score_df[,i])
    plot(riskroc,main=paste0("Risk score ROC, p < ",p_thres[i]))
    text(.4,.2,labels = paste0("AUC = ",riskroc$auc))
    dev.off()
  }
  
  
  
  
  
  #lmFit, topTable
  # top differentially expressed genes
  # odds ratios
  # testing/ROC for numbers of genes (20,50,100)
  
  ## TODO account for study and other covs

  
  
  
  # leave one out
  # ROCs
  
  ## random forest also?
  
  ## neural net
}
