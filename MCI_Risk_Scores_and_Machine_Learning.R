setwd("~/PsychGENe/brain/")

##risk scores and machine learning for MCI blood

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")

require(data.table)
require(limma)

scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")

tissues = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)
tissues = "whole_blood"

for (tissue in tissues){
  data = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_OutliersRemoved_allstudies.txt"),data.table=F)
  data = data[grep(paste0(controllabel,"|",caselabel,"$"),data$FACTOR_dx),]
  

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
  
  # run ptrs on random data
  score_df = ptrs(dat = datExpr, weight_table = weights, p_thres = c(0.5, 0.1, 0.05, 0.005, 0.002, 0.001))
  
  # plot score density
  plot(density(score_df$ptrs_1)) # column named ptrs_1 is equal to p_thres = 0.5
  
  score_df = score_df[,!(colSums(is.na(score_df))>0),drop=F]

  # contrast cases and controls on ptrs
  fit = lm(as.matrix(score_df) ~ data$FACTOR_dx)
  coefs = summary(fit)
  print(coefs$coefficients)
  # coefs = lapply(coefs, function(x) broom::tidy(x$coefficients))
  # coefs = ldply(coefs)
  # coefs = coefs[grepl("dx", coefs$.rownames), ]
  # print(coefs)
  
  
  
  
  
  
  
  
  #lmFit, topTable
  # top differentially expressed genes
  # odds ratios
  # testing/ROC for numbers of genes (20,50,100)
  
  ## create svm or something model
  # leave one out
  # ROCs
  
  ## random forest also?
}
