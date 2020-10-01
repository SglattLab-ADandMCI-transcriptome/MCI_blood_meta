setwd("~/PsychGENe/MCI_blood_meta/")
##cell abundance deconvolution for MCI blood

analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")

require(data.table)
require(plyr)
require(ggplot2)

source("~/psychgene/CIBERSORT/CIBERSORT.R")

# make the new folder
if(!dir.exists("./deconvolution/")){ dir.create("./deconvolution/")}

scalefiles = list.files("./data_for_analysis/","_GeneExpression_allstudies.txt")

tissues = sub("_GeneExpression_allstudies.txt","",scalefiles)
tissue = tissues[1]

cat("\nEstimating cell types with CIBERSORT\n")

for (tissue in tissues){
  message("Tissue: ",tissue)
  data = fread(paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"),data.table=F)
  covs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"),data.table=F)
  covs$FACTOR_age = gsub("\\+","",covs$FACTOR_age)
  studies = unique(covs$FACTOR_studyID)
  study = studies[1]
  
  wb_coef_set = list()
  for(study in studies){
    message("Study: ",study)
    smallIndex = which(covs$FACTOR_studyID == study)
    smallCovs = covs[smallIndex,]
    smallData = data[smallIndex,]
    badGenes = which(colSums(is.na(smallData)) > 0)
    smallData = smallData[,-badGenes]
    
    # Tab-delimited tabular input format with no double quotations and no missing entries.
    # HUGO gene symbols in column 1; Mixture labels in row 1
    # Data should be in non-log space. Note: if maximum expression value is <50; 
    # CIBERSORT will assume that data are in log space, and will anti-log 
    # all expression values by 2x.
    
    genes = names(smallData)
    genes[1] = "GeneSymbol"
    subjects = smallData$FACTOR_sampleID
    smallData = smallData[,-1]
    smallData = sinh(smallData)
    smallData = data.frame(subjects,smallData)
    smallData = data.frame(genes,t(smallData))
    
    fwrite(smallData,file=paste0("./deconvolution/temp.txt"),sep="\t",
           row.names = FALSE, col.names=FALSE)
    
    abundances = CIBERSORT(sig_matrix = paste0("~/psychgene/CIBERSORT/LM22.txt"),
                           mixture_file = paste0("./deconvolution/temp.txt"),
                           QN=FALSE, perm = 100)
    
    abundances = data.frame(abundances)
    fwrite(abundances,file = paste0("./deconvolution/abundances_",study,".txt"),
           row.names=TRUE, col.names=TRUE)

    wb.coef = abundances[,-c(23:25)]
    
    cat("Generating abundance plot.\n")
    png(paste("./QCplots/",tissue,"_deconvolution_",study,".png", sep =""), res=300,units="in",height=4,width=6)
    boxplot(wb.coef, las = 2, outline = TRUE, col = 'lightblue',
            pars = list(),
            main = paste0(study," deconvolution. Tissue: ", tissue))
    dev.off()

    cat("Generating linear model of abundances vs known covariates.\n")
    wb.coef = wb.coef[which(smallCovs$FACTOR_dx %in% c(caselabel,controllabel)),]
    predictors = smallCovs[which(smallCovs$FACTOR_dx %in% c(caselabel,controllabel)),]
    predictors = predictors[,which(names(predictors) %in% covariateslist)]
    predictors$FACTOR_age = as.numeric(predictors$FACTOR_age)
    predictors$FACTOR_dx = factor(predictors$FACTOR_dx)
    badPreds = colSums(is.na(predictors)) > 0
    if(any(badPreds)) predictors = predictors[,-which(badPreds)]
    badPreds = apply(predictors,2,unique)
    badPreds = ldply(badPreds, length)
    badPreds = which(badPreds$V1 == 1)
    if(any(badPreds)) predictors = predictors[,-badPreds,drop=F]
    
    design =  model.matrix( ~ .,data = predictors)
    # design = design[,-grep("CTL",dimnames(design)[[2]])]

    coefs = data.frame("NA")
    fit = try(lm(as.matrix(wb.coef)  ~ design))
    if(class(fit) != "try-error") {
      sumry = summary(fit)
      coefs = ldply(sumry, function(x) x$coefficients)
      coefs = data.frame(factor = c(rep(dimnames(design)[[2]],ncol(wb.coef))),coefs)
      coefs = coefs[!is.na(coefs$`Pr...t..`), ]
      coefs$FDR = p.adjust(coefs$`Pr...t..`, 'fdr')
      coefs$bonf = p.adjust(coefs$`Pr...t..`, 'bonferroni')
      coefs$.id = gsub("Response ", "", coefs$.id)
      wb_coef_set[[which(studies == study)]] = data.frame(StudyID = study, coefs)
    } else {
      wb_coef_set[[which(studies == study)]] = data.frame()
    }
    
    fwrite(coefs,file = paste("./deconvolution/",tissue,"_deconvolution_lm_",study,".csv", sep=""), sep=",")
    
    
    ######### this part is half finished/half broken and i'm not sure what i want to do with it
    ######### but I don't just want to delete it yet
    # cat("\nGenerating surrogate variables and confusion matrix")
    # svobj = NULL
    # foo = exprs(exprs)
    # svdf = data.frame(NULL)
    # bar = try(num.sv(foo,mod, method = 'leek'))
    # if(class(bar) != 'try-error'){
    #   svobj = sva(foo,mod, n.sv = bar)
    #   svdf = as.data.frame(svobj$sv)
    # }
    # 
    # if(ncol(svdf) > 0){
    #   colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
    #   predictors = data.frame(predictors, svdf)
    # 
    #   wb.coef = wb.coef[,colSums(wb.coef)>0]
    #   cors = cor(data.frame(svdf, wb.coef),use='pairwise.complete.obs')
    # 
    #   png(paste("./QCplots/",tissue,"_deconvolution_SVDF",paste("_", study_id[[i]]),".png", sep =""), res=300,units="in",height=6,width=6)
    #   corrplot::corrplot(cors, tl.col = 'black', number.cex = .5,
    #                      main = paste(study_id[[i]],"deconvolution. Tissue:", datExprTissue$FACTOR_tissue[which(datExprTissue$FACTOR_studyID == study_id[[i]])[1]]),
    #                      mar=c(0,0,3,0),
    #                      order='hclust',method='color', addCoef.col = "black",
    #                      addrect = 3, rect.lwd = 5,
    #                      rect.col = 'black', outline = T,
    #                      tl.srt = 45, tl.cex = 0.75)
    #   dev.off()
    # }
  }
  
  
  
  
  
  
  
  # difference between case-controls for cell types
  # names(wb_coef_set) = study_id
  if(length(wb_coef_set) > 0){
    cat("\nGenerating estimates for abundance difference between cases/controls\n")
    wb_coef_all = ldply(wb_coef_set)
    wb_coef_all = wb_coef_all[which(wb_coef_all$factor == "FACTOR_dxMCI"),]
    
    png(paste0("./QCplots/",tissue,"_deconvolution_diffMeans.png"),res=300,units="in",height=5,width=6.5)
    g = ggplot(wb_coef_all, aes(x = .id, col = StudyID, y = Estimate)) + 
      geom_point(pch = 1, position=position_dodge(0.9)) +
      theme_classic() + 
      geom_vline(xintercept = seq(from = 1.5, length.out = length(unique(wb_coef_all$.id))-1), linetype = 3, lwd = 0.3 ,col = 'navy') +
      theme(axis.text.x = element_text(hjust = 1, angle = 90), legend.text = element_text(size = 5), legend.position = 'bottom') +
      geom_hline(aes(yintercept = 0), col  ='navy', lwd = 0.3, linetype=2) +
      xlab(NULL) + 
      ylab("Estimated difference in concentration") +
      scale_color_discrete(NULL) + 
      geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error), width = 0.1, position=position_dodge(0.9))
    print(g)
    dev.off()
    
    wb_coef_all = split(wb_coef_all, wb_coef_all$.id)

    
    
    
    
    
    
    
    
    
    
    
    
    cellMeta = list()
    for(x in 1:length(wb_coef_all)){
      temp = wb_coef_all[[x]]; # subset by cell type
      if (length(temp$StudyID) < 3){
        cat("\nskipping",names(wb_coef_all)[x])
        next
      }
      weighted_meta = metafor::rma(yi =  temp$Estimate, sei = temp$Std..Error, weighted = TRUE)
      meta_beta = weighted_meta$b
      meta_se = weighted_meta$se
      meta_pval = weighted_meta$pval
      cellMeta[[x]] = data.frame(CellType = unique(temp$.id), Beta = meta_beta, SE = meta_se, P = meta_pval)
    }
    cellMeta = ldply(cellMeta)
    cellMeta$FDR = p.adjust(cellMeta$P, 'fdr')
    cellMeta$bonf = p.adjust(cellMeta$P, 'bonferroni')
    cellMeta
    write.csv(cellMeta, file = paste0("deconvolution/",tissue,"_deconvolution_meta-sumstats.csv",sep=""),row.names = FALSE)
  }
}
