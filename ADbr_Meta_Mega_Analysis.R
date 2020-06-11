setwd("~/PsychGENe/brain/")

## meta analysis for brain data
## GCH w/JH and WB

# load these packages (install if needed)
require(plyr)
require(ggplot2)
require(data.table)
library(limma)
require(BRETIGEA)
# require(dtangle)
# require(dtangle.data)
## https://umich.app.box.com/v/dtangledatapkg
require(AnnotationDbi)
require(sva)
require(hgu133a.db)
require(metafor)

metafolder = "./meta_analysis"
megafolder = "./mega_analysis"
forestfolder = "./forestplots"

if(!dir.exists(metafolder)){
  dir.create(metafolder)
}

if(!dir.exists(forestfolder)){
  dir.create(forestfolder)
}

datfiles = list.files("./data_for_analysis/","_GeneExpression_allstudies.txt")
covfiles = list.files("./data_for_analysis/","_SampleFactors_allstudies.txt")
scalefiles = list.files("./data_for_analysis/","_ScaledWithFactors_OutliersRemoved_allstudies.txt")

dattis = sub("_GeneExpression_allstudies.txt","",datfiles)
covtis = sub("_SampleFactors_allstudies.txt","",covfiles)
scaletis = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

if(all(dattis == covtis, covtis == scaletis)){
  tissues = dattis
} else {
  stop("Tissue files are... not right.")
}

for(tissue in tissues){
  ## DGE per tissue per study, non-scaled data
  datExprTissue = fread(paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"), data.table=F, stringsAsFactors = F)
  datExprCovs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"), data.table=F, stringsAsFactors = F)
  datExprCovs$FACTOR_age = as.numeric(sub("\\+","",datExprCovs$FACTOR_age))
  genes = datExprTissue$SYMBOL
  datExprTissue = data.frame(datExprCovs,t(datExprTissue[,-1]))
  names(datExprTissue) = c(names(datExprCovs),genes)
  datExprTissue = datExprTissue[grep("CTL|AD",datExprTissue$FACTOR_dx),]
  datExprTissue$FACTOR_dx = factor(datExprTissue$FACTOR_dx, levels = c("CTL","AD"))
  
  study_id = unique(datExprTissue$FACTOR_studyID)
  save_results = list()
  wb_coef_set = list()
  for( i in 1:length(study_id)){
    
    cat("\nStudy-wise differential expression analysis:",tissue,i,"~",study_id[[i]])
    expr.tmp = datExprTissue[datExprTissue$FACTOR_studyID %in% study_id[[i]],]
    
    N = table(expr.tmp$FACTOR_dx)
    Nca = N[[2]]
    Nco = N[[1]]
    cat("\n   Detected", Nca,"cases and",Nco,"controls")
    ## TODO need info on what we are doing here if we throw anything out
    
    # IF SAMPLE SIZE IS TOO SMALL, SKIP TO NEXT STUDY
    if(Nco < 3 | Nca < 3){
      cat("\nSkipping for small sample size.")
      next
    }
    
    # calculate missingness column-wise
    colmisrate = colSums(is.na(expr.tmp[,!grepl("FACTOR_",colnames(expr.tmp))]))/nrow(expr.tmp) # proportion of missingness
    colexclude = colmisrate[colmisrate > 0.0001] # missingness filter 
    
    expr.tmp = expr.tmp[,!colnames(expr.tmp) %in% names(colexclude)]
    
    # statistical analysis (via linear model)
    y = expr.tmp[,colnames(expr.tmp) %in% names(colmisrate)]
    x = expr.tmp[,!colnames(expr.tmp) %in% names(colmisrate)]
    
    # identify columns in y that are non-numeric
    class_check = lapply(y,class)
    class_check = unlist(class_check)
    wrong_class = which(class_check != "numeric")
    
    # add wrong class to demographics
    if(length(wrong_class)>0){
      x = data.table(x, y[,colnames(y) %in% names(wrong_class)]);
      y = y[,!colnames(y) %in% names(wrong_class)]}
    
    # remove columns with missingness ( > 50%)
    colmis = colSums(is.na(x)) >= nrow(x)*.20
    colmis = which(colmis == TRUE)
    if(length(colmis) > 0){x = x[,!colnames(x) %in% names(colmis)]}
    
    # remove rows with some missingness
    rowmis = rowSums(is.na(x))
    rowmis = which(rowmis > 0)
    if(length(rowmis) > 0){x = x[-rowmis]; y = y[-rowmis]}
    
    # Model matrix (basic)
    PredListNames = paste('FACTOR_', c("dx", "sex", "age","ethnicity"), collapse= "|", sep="")
    predictors = x[,grep(PredListNames, colnames(x))]
    predictors = data.frame(predictors)
    
    # set level of disease groups 
    predictors$FACTOR_dx = as.factor(predictors$FACTOR_dx)
    predictors$FACTOR_dx = relevel(predictors$FACTOR_dx, ref = 'CTL')
    
    counts = lapply(predictors,unique)
    count_class = unlist(lapply(counts,length))
    count_class = count_class[count_class <= 1,drop=F]
    
    if(length(count_class) > 0){predictors = predictors[,!(colnames(predictors) %in% names(count_class))]}
    
    N = table(predictors$FACTOR_dx)
    Nca = N[[2]]
    Nco = N[[1]]
    cat("\n   Kept", Nca,"cases and",Nco,"controls")
    ##TODO do i need subject ids in here somewhere
    
    { # deconvolution SECTION
      ## TODO fix this for brainz instead of leukocytes lol
      # estimate leukocyte abundance
      message("\nEstimating cell types with deconvolution...")
      ##TODO evaulate if this is producing sensible/correct results
      
      exprs = as.data.frame(y) # extract normalized gene expression intensities for subjects
  
      # ##LM22
      # foo = which(names(exprs) %in% colnames(dtangle.data::newman_pbmc$data$log))
      # exprs = exprs[,foo]
      # foo = which(colnames(dtangle.data::newman_pbmc$data$log) %in% names(exprs))
      # ps = dtangle.data::newman_pbmc$annotation$pure_samples
      # refs = dtangle.data::newman_pbmc$data$log[,foo]
      # foo = sapply(strsplit(names(ps)," "),"[",1)
  
      # ## GSE28492
      # refs = fread("./references/brain_deconvolution/deconvolutionSamples.txt", data.table=F)
      # rownames(refs) = refs$V1
      # refs = refs[,-1]
      # foo = which(colSums(is.na(refs)) > 0)
      # refs = refs[,-foo]
      # foo = which(names(exprs) %in% colnames(refs))
      # exprs = exprs[,foo]
      # foo = numeric()
      # for(asdf in 1:length(names(exprs))){
      #   foo[asdf] = which(names(refs) == names(exprs)[asdf])
      # }
      # refs = refs[,foo]
      # ps = as.list(c(1:nrow(refs)))
      # names(ps) = rownames(refs)
      # foo = sapply(strsplit(names(ps),"\\."),"[",1)
      # 
      # ## Collapse cell types
      # baz = list()
      # for(x in foo){
      #   baz[[x]] = c(unlist(ps[foo == x]))
      # }
      # ps = baz
      # 
      # ## dtangle(row=sample col=gene, uses top10% of genes, data_type adjusts gamma fuzzfactor, )
      # ## TODO data_type... by studyid? i guess just change if rnaseq.
      # data_type = "microarray-gene"
      # 
      # res = dtangle(exprs, references = refs,
      #               pure_samples = ps,
      #               data_type=data_type,
      #               n_markers = .2,
      #               marker_method = "ratio")
      
      exprs = data.frame(t(exprs))
      names(exprs) = rownames(y)
      rownames(exprs) = names(y)
      res = brainCells(exprs)
      
      # wb.coef = res$estimates
      # bar = which(rownames(wb.coef) %in% rownames(refs))
      # if(length(bar) > 0) wb.coef = wb.coef[-bar,]
      
      wb.coef = res
      
      png(paste("./QCplots/",tissue,"_deconvolution_",study_id[[i]],".png", sep =""), res=300,units="in",height=4,width=6)
      boxplot(wb.coef, las = 2, outline = TRUE, col = 'lightblue', main = paste(study_id[[i]],"deconvolution. Tissue:", datExprTissue$FACTOR_tissue[which(datExprTissue$FACTOR_studyID == study_id[[i]])[1]]))
      dev.off()
      
      ## PCA reduction of peripheral leukocyte proportions
      design =  model.matrix( ~ -1 + ., data = predictors)
      design = design[,!colnames(design) %in% c("FACTOR_dxAD")]
      
      coefs = data.frame("NA")
      fit = try(lm(as.matrix(wb.coef)  ~ design))
      if(class(fit) != "try-error") {
        fit = summary(fit)
        coefs = lapply(fit, function(x) broom::tidy(x$coefficients))
        coefs = ldply(coefs)
        coefs = coefs[grepl("FACTOR_dx", coefs$.rownames),]
        coefs = coefs[!is.na(coefs$Pr...t..), ]
        coefs$FDR = p.adjust(coefs$Pr...t.., 'fdr')
        coefs$.id = gsub("Response ", "", coefs$.id)
        wb_coef_set[[i]] = data.frame(StudyID = study_id[[i]], coefs)
      } else {
        wb_coef_set[[i]] = data.frame()
      }
      
      fwrite(coefs,file = paste(metafolder,"/",tissue,"_deconvolution_",study_id[[i]],".csv", sep=""), sep=",")

    }
    
    # remove genes with low variance
    var_filter = lapply(y, sd)
    low_var_filter = var_filter[is.na(var_filter) | var_filter < .001]
    
    if(length(low_var_filter) > 0){y = y[,!colnames(y) %in% names(low_var_filter)]}
    
    ## Surrogate variable analysis - default method 
    exprs = ExpressionSet(as.matrix(t(y)))
    # predictors = predictors[,colnames(predictors) %in% c("FACTOR_Age", "FACTOR_dx", "FACTOR_Sex", "FACTOR_Batch")ALSE]
    mod = model.matrix(~ ., data= predictors) # model with known factors and covariates
    mod0 = model.matrix(~1,data=predictors) # intercept only model
    
    svobj = NULL
    foo = exprs(exprs)
    
    ##TODO wb.coef is called outside the if here...
    ##TODO add some cat() so you know what's happening
    
    svdf = data.frame(NULL)
    bar = try(min(c(num.sv(foo,mod, method = 'leek'), num.sv(foo,mod, method = 'be'))))
    if(class(bar) != 'try-error'){
      svobj = sva(foo,mod, n.sv = bar)
      svdf = as.data.frame(svobj$sv)
    }
      
    if(ncol(svdf) > 0){
      colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
      predictors = data.frame(predictors, svdf)
      
      wb.coef = wb.coef[,colSums(wb.coef)>0]
      cors = cor(data.frame(svdf, wb.coef),use='pairwise.complete.obs')
      
      png(paste("./QCplots/",tissue,"_deconvolution_SVDF",paste("_", study_id[[i]]),".png", sep =""), res=300,units="in",height=6,width=6)
      corrplot::corrplot(cors, tl.col = 'black', number.cex = .5,
                         main = paste(study_id[[i]],"deconvolution. Tissue:", datExprTissue$FACTOR_tissue[which(datExprTissue$FACTOR_studyID == study_id[[i]])[1]]),
                         mar=c(0,0,3,0),
                         order='hclust',method='color', addCoef.col = "black",
                         addrect = 3, rect.lwd = 5,
                         rect.col = 'black', outline = T,
                         tl.srt = 45, tl.cex = 0.75)
      dev.off()
    }
    
    # final design matrix for differential expression analysis
    # predictors$FACTOR_dx = as.factor(predictors$FACTOR_dx)
    # predictors$FACTOR_dx = relevel(predictors$FACTOR_dx, ref = 'CTL')
    # design = model.matrix( ~ -1 + ., predictors)
    design = model.matrix( ~ ., predictors)
    design = design[,!grepl("CTL", colnames(design))]
    
    # fit the linear model (log2 expression as response variable)
    lmFit = lm(as.matrix(y) ~ -1 + design)
    
    # extract summary statistics
    summary = summary(lmFit)
    # stdSummary = summary(stdLmFit)
    # format summary statistics into a table
    tidy_table = lapply(summary, function(x) broom::tidy(x))
    names(tidy_table) = colnames(y)
    # tidy_stdTable = lapply(stdSummary, function(x) broom::tidy(x))
    # names(tidy_stdTable) = colnames(y)
    
    # squash into a big table
    big_table = ldply(tidy_table)
    big_table$N_cases = Nca
    big_table$N_controls = Nco
    big_table$N = nrow(x)
    # stdTable = ldply(tidy_stdTable)
    # big_table$BetaStd = stdTable$estimate
    # big_table$SEStd = stdTable$std.error
    
    nGenes = length(unique(big_table$.id))
    
    names(big_table)[names(big_table) %in% ".id"] = "GeneSymbol"
    save_results[[i]]  = big_table
    names(save_results)[[i]] = study_id[[i]]
    
  }
  
  
  # difference between case-controls for cell types
  # names(wb_coef_set) = study_id
  if(length(wb_coef_set) > 0){
    ##TODO a cat or message here
    wb_coef_all = ldply(wb_coef_set)
    
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
      if (length(temp$StudyID) < 5){
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
    cellMeta
    write.csv(cellMeta, file = paste(metafolder,"/",tissue,"_deconvolution_meta-sumstats.csv",sep=""),row.names = FALSE)
  }
    
  ## merge study-wise sum stats into table
  ## TODO a cat or message here
  mergestats = ldply(save_results)
  names(mergestats)[names(mergestats) %in% ".id"] = "studyID"
  mergestats$term = gsub("designFACTOR_|design", "", mergestats$term)
  
  unique(mergestats$studyID) ## TODO fix up this with years
  unique(mergestats$term)
  # mergestats$studyID[grepl("AddNeuroMed1", mergestats$studyID)] = "AddNeuroMed Batch 1 (year)"
  # mergestats$studyID[grepl("AddNeuroMed2", mergestats$studyID)] = "AddNeuroMed Batch 2 (year)"
  # mergestats$studyID[grepl("ADNI", mergestats$studyID)] = "ADNI (year)"
  # mergestats$studyID[grepl("Bai", mergestats$studyID)] = "Bai et al. (year)"
  # mergestats$studyID[grepl("Chen", mergestats$studyID)] = "Chen et al. (year)"
  # mergestats$studyID[grepl("Leandro", mergestats$studyID)] = "Leandro et al. (year)"
  # mergestats$studyID[grepl("Maes", mergestats$studyID)] = "Maes et al. (year)"
  # mergestats$studyID[grepl("Samsudin", mergestats$studyID)] = "Samsudin et al. (year)"
  # mergestats$studyID[grepl("Scherzer", mergestats$studyID)] = "Scherzer et al. (year)"
  
  
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  genes = data.frame(genes)
  genes$SYMBOL = select(org.Hs.eg.db,keys=as.character(genes$gene_id), keytype='ENTREZID', columns="SYMBOL")$SYMBOL
  genes$SYMBOL = gsub("[-]", ".", genes$SYMBOL)
  
  mergestats = mergestats[mergestats$GeneSymbol %in% genes$SYMBOL, ]
  
  fwrite(mergestats,
         file = paste(metafolder, "/",tissue,"_AD_brain_perStudyDGEsummary.csv",sep=""),
         quote = F, row.names= F,sep=",")
  
  
  ## TODO fix these graph thingies
  graph_df = mergestats
  # graph_df$term[graph_df$term %in% "psychosisYes"] = "Psychosis (Yes/No)"
  graph_df$term[graph_df$term %in% "gendermale"] = "Sex (Female/Male)"
  # graph_df$term[graph_df$term %in% "tobaccoYes"] = "Tobacco (Yes/No)"
  graph_df$term[graph_df$term %in% "dxAD"] = "Diagnosis (AD/CTL)"
  graph_df$term[graph_df$term %in% "ethnicitywhite"] = "White (Yes/No)"
  graph_df$term[graph_df$term %in% "age"] = "Age (years)"
  graph_df$term[graph_df$term %in% "tissue"] = "Sample Collected"
  # graph_df$term[graph_df$term %in% "Batch"] = "Array batch"
  non_sv = graph_df[!grepl("SV",graph_df$term), ]
  sv_df = graph_df[grepl("SV",graph_df$term), ]
  sv_df$term = factor(sv_df$term,levels=paste("SV",1:20,sep=""))
  
  combn = ldply(list(non_sv,sv_df))
  combn$term = factor(combn$term, levels=unique(combn$term))
  
  g = ggplot(combn[!grepl("Intercept", ignore.case=T, combn$term), ],
              aes(x = term, y = -log10(p.value), fill = factor(term))) +
              ylab(expression(paste("Differential expression, -log"[10],"(P-value)"))) + 
              facet_wrap(~studyID, ncol = 4) +
              xlab(NULL) + 
              geom_violin() + 
              guides(fill = FALSE) + 
              theme_bw() +
              geom_hline(yintercept = -log10(.05), col = "red", linetype ="dashed", lwd = 0.2) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
              geom_boxplot(width=0.05, outlier.shape = NA, fill = "lightgrey", col = "black")
  
  png(file = paste("./QCplots/",tissue,"_qcsva_ADbrMEGA_MERGESTUDY_violin.png", sep = ""),
      res = 300, units = "in", height = 8, width = 11.5)
  print(g)
  dev.off()
  
  
  agg_stats = data.table(combn)
  agg_stats = agg_stats[,list(LOGP = mean(-log10(p.value)), SD = sd(-log10(p.value)), K = length(p.value)),by=list(term, studyID)]
  agg_stats = agg_stats[!grepl("Intercept", agg_stats$term)]
  agg_stats$CI_LOW = agg_stats$LOGP - (1.96*agg_stats$SD/sqrt(agg_stats$K))
  agg_stats$CI_HIGH = (1.96*agg_stats$SD/sqrt(agg_stats$K)) + agg_stats$LOGP
  
  png(paste0("./QCplots/",tissue,"meanLOGP_term_DGE.png"),res=300,units="in",height=10,width=10)
  print(
    ggplot(agg_stats, aes(x = term, y = LOGP, fill = term)) +
    facet_wrap(~studyID, ncol = 3, scale = 'free_x') +
    geom_bar(stat = 'identity',col = 'black') +
    guides(fill = F) +
    theme_classic() +
    ylab(expression(paste("-log"[10],"(P-value)"))) + 
    xlab(NULL) +
    geom_errorbar(aes(ymin = CI_LOW, ymax = CI_HIGH),width=0.2) +
    theme(axis.text.x = element_text(hjust = 1, angle = 90), axis.title=element_text(size = 15)) + 
    geom_hline(aes(yintercept = -log10(.05)), col = 'red', linetype=2)
    )
  dev.off()
  
  
  # Random effect meta-analysis
  
  ## perform meta-analysis of genes (REML model)
  ## TODO cat or message here
  
  genes = unique(mergestats$GeneSymbol)
  genes = genes[!grepl("_", genes)]
  
  res_save = list()
  loo_save = list()
  for( i in 1:length(genes)){
    
    cat("\rMeta-analysis:",i,"/",length(genes))
    sub = mergestats[mergestats$GeneSymbol %in% genes[[i]], ]
    sub = sub[grepl("dx", sub$term), ]
    
    sub_combn = data.frame(estimate = c(sub$estimate),
                           # BetaStd = c(sub$BetaStd),
                           SE = c(sub$std.error),
                           # SEstd = c(sub$SEStd),
                           # N = c(sub$N),
                           N_case = c(sub$N_cases),
                           N_control = c(sub$N_controls),
                           studyID = sub$studyID)
    
    if(nrow(sub_combn) < 2) next
    
    
    if(nrow(sub) <= 1) next
    
    direction = sign(sub_combn$estimate)
    direction = ifelse(direction == 1, "+", "-")
    direction = paste(direction, collapse= "", sep = "")
    
    res = NA
    
    res = try(suppressWarnings(metafor::rma(yi = sub_combn$estimate,
                                            sei = sub_combn$SE, 
                                            verbose = FALSE,
                                            weighted = TRUE, method = "DL")))
    
    # res.std = try(suppressWarnings(metafor::rma(yi = sub_combn$BetaStd,
    #                                             sei = sub_combn$SEstd, 
    #                                             verbose = FALSE,
    #                                             weighted = TRUE, 
    #                                             method = "DL")))
    
    # try to obtain results with fewer input studies
    if( unique(grepl("Error",res)) == T ) {
      
      cat(genes[[i]],"\nsub_combn = sub_combn[sub_combn$N > 30, ];\n")
      sub_combn = sub_combn[sub_combn$N > 30, ];
      
      res = try(suppressWarnings(metafor::rma(yi = sub_combn$estimate,
                                              sei = sub_combn$SE, 
                                              verbose = FALSE,
                                              weighted = TRUE, 
                                              method = "DL")));
      
      # res.std = try(suppressWarnings(metafor::rma(yi = sub_combn$BetaStd,
      #                                             sei = sub_combn$SEstd, 
      #                                             verbose = FALSE,
      #                                             weighted = TRUE, 
      #                                             method = "DL")))
      
    } 
  
    if( unique(grepl("Error",res)) == T ){
      cat("Giving up on",i,"\n")
      next
    }
    
    if(nrow(sub_combn) > 1){
      loo_tmp = leave1out(res)
      loo_save[[i]] = data.frame(studyID = sub_combn$studyID, loo_tmp, GeneSymbol = genes[[i]])
    }
    
    res = data.frame(GeneSymbol = sub$GeneSymbol[[1]],
                     Direction = direction, 
                     Nstudy = nrow(sub_combn),
                     Nsample = sum(sub_combn$N),
                     Ncase = sum(sub_combn$N_case),
                     Ncontrol = sum(sub_combn$N_control),
                     Log2FC = res$beta, 
                     SE = res$se, 
                     # BetaStd = res.std$beta,
                     # SEstd = res.std$se,
                     P = res$pval,
                     HetQ = res$QE,
                     HetP = res$QEp,
                     study_kept = paste(sub_combn$studyID, collapse =  ", ", sep = ""))
    
    ##reduce Nstudy RE both ANM batches
    if(length(grep("AddNeuroMed",sub_combn$studyID)) > 1){
      res$Nstudy = res$Nstudy - 1
    }
    
    res_save[[i]] = res
    names(res_save)[i] = genes[i]
  }
  
  foo = llply(res_save,class)
  foo = unlist(foo)
  res_df = ldply(res_save[which(foo == "data.frame")])
  
  min(res_df$P)
  
  fwrite(data.table(res_df), 
         sep = "\t",
         file = paste(metafolder,"/",tissue,"_prefreeze_qcAD_brain_meta.txt", sep=""),
         quote  = F, row.names = F)
  write(names(res_save)[which(foo != "data.frame")],
        file = paste(metafolder,"/",tissue,"_prefreeze_qcAD_brain_meta_NULLmetaforRMA.txt", sep=""))
  
  
  # ab = res_df[res_df$GeneSymbol %in% gene_filter$GeneSymbol, ]
  
  ## at least 2 studies are contributing to the gene
  res_df = res_df[res_df$Nstudy >= 2, ]
  res_df = res_df[order(res_df$P,decreasing=F),]
  res_df$FDR = p.adjust(res_df$P, "fdr")
  res_df$BONF = p.adjust(res_df$P, "bonferroni")
  res_df$GeneSymbol = gsub("[.]", "-", res_df$GeneSymbol)
  
  dim(res_df)
  
  head(res_df)
  
  # Location of genes 
  require(Homo.sapiens)
  locs <- select(Homo.sapiens, keys=as.character(res_df$GeneSymbol),keytype="SYMBOL", columns=c("ENTREZID", "TXCHROM","TXSTART","TXEND","TXSTRAND"))
  locs$WIDTH = abs(locs$TXSTART - locs$TXEND)
  locs = locs[order(locs$WIDTH,decreasing = T), ]
  locs = locs[!duplicated(locs$SYMBOL), ]
  locs$LOC = paste(locs$TXCHROM, ":",locs$TXSTART,"-",locs$TXEND,"(",locs$TXSTRAND,")",sep="")
  locs = locs[!is.na(locs$TXCHROM),]
  loc_df = locs[,c("SYMBOL","LOC", "ENTREZID")]
  colnames(loc_df)[1] = "GeneSymbol"
  
  
  # Combine results with chromosome map
  res_df = merge(loc_df, res_df, by="GeneSymbol",all.y=T)
  res_df = res_df[order(res_df$P, decreasing = F), ]
  
  fwrite(data.table(res_df), 
         sep = "\t",
         file = paste(metafolder,"/",tissue,"_freeze_qcAD_brain_meta_Nmin4_LeaveOneOut.txt",sep=""), quote  = F, row.names = F)
  
  
  # qq-plot
  # association p-values
  assoc = data.frame(P = res_df$P, source='DGE Meta-analysis')
  # this is not the solution
  # assoc = data.frame(P = res_df$BONF, source='DGE Meta-analysis')
  observed = assoc$P
  chisq1 <- qchisq(1-observed, 1)
  observed <- sort(assoc$P)
  lobs <- -(log10(as.numeric(observed)))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  medianchi = median(chisq1,na.rm=T)
  lambdaout = medianchi/.454
  
  assoc = assoc[order(assoc$P, decreasing = F), ]
  assoc$lobs = lobs
  assoc$lexp = lexp
  
  ci = .95
  N = length(assoc$P)
  observed = -log10(sort(assoc$P))
  expected = -log10(1:N / N)
  clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
  cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  assoc$clower = clower
  assoc$cupper = cupper
  
  # heterogeneity p-values
  het = data.frame(P = res_df$HetP, source='Heterogeneity')
  
  observed <- sort(het$P)
  lobs <- -(log10(as.numeric(observed)))
  
  # qqplot statistics
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  
  het = het[order(het$P, decreasing = F), ]
  het$lobs = lobs
  het$lexp = lexp
  
  ci = .95
  N = length(het$P)
  observed = -log10(sort(het$P))
  expected = -log10(1:N / N)
  clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
  cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  het$clower = clower
  het$cupper = cupper
  
  set = ldply(list(assoc,het))
  
  mclab = substitute(paste("Median ", chi^2, "=", MC, ", ", lambda , "=", LD), list(MC=format(medianchi, digits = 3), LD=format(lambdaout,digits=3)))
  
  qplot = ggplot(assoc, aes(x = lexp, y = lobs)) +
    geom_point(colour="black", fill= 'dodgerblue', shape = 21, size = 2, stroke = 0.3) + 
    geom_abline(slope = 1, intercept = 0, col = 'black', lwd = 0.3) +
    theme_classic() + 
    ylab(expression(paste("Observed -log"[10]," P-value"))) +
    xlab(expression(paste("Expected -log"[10]," P-value"))) +
    scale_fill_discrete(NULL) +
    geom_line(aes(x = lexp, y = clower), colour ='grey' ,lwd = 0.75) +
    geom_line(aes(x = lexp, y = cupper), colour = 'grey', lwd = 0.75) +
    geom_ribbon(aes(x = lexp, ymin = clower, ymax = cupper), fill="grey", alpha="0.2") +
    ggtitle(mclab)  +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
  
  # qplot
  png(paste0("./QCplots/",tissue,"_qcADbr_meta_qqplot.png"),res=300,units="in",height=5,width=5)
  print(qplot)
  dev.off()
  
  
  ## Manhattan plot
  trim = res_df[!is.na(res_df$LOC),]
  trim = trim[grepl("[:]", trim$LOC), ]
  locs = strsplit(trim$LOC, "[:]")
  chr = lapply(locs, function(x) x[[1]])
  bp = lapply(locs, function(x) x[[2]])
  pos = data.frame(CHR=unlist(chr),BP=unlist(bp))
  pos$GeneSymbol = trim$GeneSymbol
  pos$bonf = trim$BONF
  pos$fdr = trim$FDR
  pos$P = trim$P
  pos$BP = gsub("[(-)]", "", pos$BP)
  pos$BP = gsub("[(+)]", "", pos$BP)
  pos$START = unlist(lapply(strsplit(pos$BP, "[-]"), function(x)x[[1]]))
  pos$END = unlist(lapply(strsplit(pos$BP, "[-]"), function(x)x[[2]]))
  pos$CHR = gsub("chr","",pos$CHR)
  pos$CHR[pos$CHR %in% "X"] = 23
  pos$CHR[pos$CHR %in% "Y"] = 24
  
  col = data.frame(CHR = 1:24, col = c("royalblue", "orange"))
  
  sub = merge(pos, col, by = "CHR")
  sub$CHR = as.integer(sub$CHR)
  sub$START = as.integer(sub$START)
  sub$pos = NA
  sub$pos = ifelse(sub$CHR  == 1, as.integer(sub$START), NA)
  
  median_pos = list()
  chr_grab = unique(col$CHR)
  for( i in 1:length(chr_grab)){
    
    if(i == 1){
      prior_max = min(sub$pos[sub$CHR == 1]) - 1
    }
    
    if(i > 1){
      k = i - 1
      prior_max = max(sub$pos[sub$CHR == k])
    }
    
    chr_seq = sub$START[sub$CHR == chr_grab[[i]]]
    true_min = min(chr_seq)
    diff = chr_seq - true_min
    new_chr_seq = diff + prior_max + 1
    
    sub$pos[sub$CHR == i] <- new_chr_seq
    
    median_pos[[i]] =  (min(sub$pos[sub$CHR %in% chr_grab[[i]]]) + max(sub$pos[sub$CHR %in% chr_grab[[i]]]))/2
  }
  
  names(median_pos) = 1:24
  
  median_pos = ldply(median_pos)
  median_pos$.id[median_pos$.id %in% c(23,24)] = c("X","Y")
  
  
  sub = sub[order(sub$CHR, sub$pos), ]
  sub$SYMBOL = ifelse(sub$bonf < .05, sub$GeneSymbol, NA)
  
  require(ggrepel)
  
  png(paste0("./QCplots/",tissue,"_qcADvCTL_brain_ggplot_manhattan.png"), res = 300, units = "in", height = 5, width = 9)
  # manhattan plot
  print(
  try(
    ggplot(sub, aes(x = sub$pos, y = -log10(sub$P))) + 
    geom_point(size = 0.7, col = as.character(sub$col)) + 
    ylim(min = 0, max = 1.2*max(-log10(sub$P))) + 
    xlab("Genomic coordinate") + 
    theme_classic() +
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    scale_x_continuous(name="Genomic coordinate", breaks=median_pos$V1, labels=median_pos$.id) +
    geom_hline(yintercept = -log10(.05/nrow(sub)), col = "black", lwd = 0.5, linetype = "dashed") +
    geom_text_repel(aes(label = sub$SYMBOL), fontface = 'italic', size = 3, col = "black") +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))
  # end plot
  ))
  dev.off()
  
  
  # volcano plot
  res_df$FDR = p.adjust(res_df$P, "fdr")
  res_df$BONF = p.adjust(res_df$P, "bonferroni")
  psize = -log10(res_df$P)
  psize = (psize - min(psize))/(max(psize) - min(psize)) + 0.5
  col = rep("lightgrey", nrow(res_df))
  col[which(abs(res_df$Log2FC) > 0.2)] = "dodgerblue3"
  col[which(res_df$FDR < .05)] = "orange"
  col[which(res_df$BONF < .05)] = "firebrick3"
  
  res_df$vLabel = ifelse(abs(res_df$Log2FC) > .4 | res_df$BONF < .05, res_df$GeneSymbol, NA)
  res_df$vLabel = gsub("[.]", "-", res_df$vLabel)
  
  signCounts = table(sign(res_df$Log2FC))
  paste("n = ", signCounts, sep = "")
  
  png(paste0("./QCplots/",tissue,"_qcADbr_meta-volcano.png"),res=300,units="in",height = 5, width = 6)
  print(
    ggplot(res_df, aes(x = Log2FC, y = -log10(P))) + 
    geom_point(size = psize, col = col) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), panel.border = element_rect(fill=NA,size = 1.5)) +
    geom_hline(yintercept = -log10(.05), col = 'orange', linetype = "dashed") + 
    xlab( expression(paste("log"[2]," fold-change in gene expression"))) +
    ylab( expression(paste("-log"[10],"(P-value)"))) + 
    geom_hline(aes(yintercept = -log10(.05/nrow(res_df))), col = 'red', lwd = 0.25,linetype=1) +
    geom_vline(aes(xintercept = 0), col = 'black',lwd=0.3, linetype=2) +
    ggrepel::geom_text_repel(aes(label = res_df$vLabel), size = 3, segment.size = 0.2, segment.colour = "grey", fontface = 3)
  )
  dev.off()
  
  top_df = res_df[!is.na(res_df$vLabel),]
  
  fwrite(top_df,paste(metafolder,"/",tissue,"_brain_meta_significant_log2_and_pvals.csv", sep=""))
  
  
  ## Forest plots 
  pdf(paste(forestfolder,"/",tissue,"_FORESTPLOT_all.pdf",sep=""))
  
  res_filter = res_df[order(res_df$P,decreasing=F), ]
  topgenes = unique(res_filter$GeneSymbol[res_filter$FDR < .1])
  
  mergestats$GeneSymbol = gsub("[.]", "-", mergestats$GeneSymbol)
  
  if(length(topgenes)>1) for( i in 1:length(topgenes)){
    
    # subset meta-analysis results
    meta = res_filter[res_filter$GeneSymbol %in% topgenes[[i]], ]
    keep_studies = unlist(strsplit(as.character(meta$study_kept), ", "))
    
    meta = data.frame(Log2FC = meta$Log2FC,
                      P = meta$P, SE = meta$SE, 
                      studyID = "Pooled result", 
                      Ncase = meta$Ncase, 
                      Ncontrol = meta$Ncontrol,
                      N = meta$Nsample,
                      col = "blue", 
                      Source = "Meta-analysis")
    
    # first grab the study-specific expression profiles
    sub = mergestats[mergestats$GeneSymbol %in% topgenes[[i]], ]
    sub = sub[grepl("dx", sub$term), ]
    sub = sub[sub$studyID %in% keep_studies, ] # retain correct studies for forest plot
    
    
    if(nrow(sub) > 0){
      
      indiv = data.frame(Log2FC = sub$estimate, 
                         P = sub$p.value,
                         SE = sub$std.error, 
                         studyID = sub$studyID,
                         Ncase = sub$N_cases, 
                         Ncontrol = sub$N_controls,
                         # N = sub$N,
                         col = "red", 
                         Source = "Studies");
      
      results = ldply(list(indiv,meta))} else{results = meta}
    
    results$P_label = paste("P = ", format(scientific=T,digits=3,results$P), sep = "")
    results$CI_LOW = results$Log2FC - (1.96 * results$SE)
    results$CI_HIGH = results$Log2FC  + (1.96 * results$SE)
    
    locus = res_df$LOC[res_df$GeneSymbol %in% topgenes[[i]]]
    
    g = ggplot(results, aes(x = Log2FC, y = (studyID))) + 
      geom_point(col = results$col) + 
      theme_bw() +
      xlim(min(floor(results$CI_LOW)), max(ceiling(results$CI_HIGH))) +
      xlab(expression(paste("Log"[2], " fold change (95% CI)"))) +
      ylab(NULL) +
      ggtitle( topgenes[[i]]) +
      geom_vline(aes(xintercept = 0),linetype="dashed", col = "grey") +
      facet_grid(Source~., scales = "free_y", space = "free_y") + 
      geom_errorbarh(xmin = results$CI_LOW, xmax = results$CI_HIGH, col = results$col, height = NA) + 
      theme( strip.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
      geom_text(aes(vjust = -1, label = results$P_label))
    
    # png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/brain/forestplots/RANK_",i,"_FORESTPLOT_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
    print(g)
    # dev.off()
    
  } else print("No top genes.")
  dev.off()
  
  
  
  ## Directionality gene score
  res_df = res_df[order(res_df$P, decreasing = F), ]
  
  majority = max(table(sign(res_df$Log2FC)))
  broom::tidy(binom.test(x = majority, n = nrow(res_df), p = 0.5, alternative = 'two.sided'))
  
  direction = sign(res_df$Log2FC)
  dirAll = sign(sum(sign(res_df$Log2FC)))
  
  binomStats = as.list(rep(NA, length(direction)))
  for( i in 1:length(direction)){
    cat("\rDirection test:",round(i/length(direction)*100,2),"%         ")
    if(i == 1) next
    tmp = direction[1:i]
    sign_count = sum(tmp == dirAll) # total number of up(or down)-regulated genes
    propneg = sign_count/length(tmp) # total number of genes in bin
    binomStats[[i]] = broom::tidy(binom.test(x = sign_count, n = length(tmp), p = 0.5, "greater")) # binomial sign test 
  }
  binomStats_df = ldply(binomStats[-1])
  binomStats_df$rank = 1:nrow(binomStats_df) ## rank of gene expression p value
  binomStats_df = binomStats_df[order(binomStats_df$p.value,decreasing=F),] ## now sort by binomail test p value
  binomStats_df$label = NA
  binomStats_df$label[[1]] = paste("p = ",format(binomStats_df$p.value[[1]],digits=2,scientific=T),sep="")
  binomStats_df = binomStats_df[order(binomStats_df$rank,decreasing=F),]
  binomStats_df$bin = dplyr::ntile(binomStats_df$parameter, 1)
  cat("\n")
  
  png(paste("./QCplots/",tissue,"_AD_brain_meta_signtest.png",sep=""), res = 300, units = "in", height = 5, width = 7.5)
  print(
    ggplot(binomStats_df, aes(x = parameter, y = estimate)) + 
    geom_line(col = 'royalblue') + 
    facet_wrap(~bin, scales = 'free_x', ncol = 5) + 
    ggtitle(paste("Sign test, ",ifelse(dirAll>0,"up","down"),"-regulated genes.",sep="")) +
    theme_classic() + 
    theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
    geom_hline(aes(yintercept = table(direction)[[(dirAll+3)/2]]/length(direction)), linetype=2,col='red')
  )
  dev.off()
  
  gc()
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  # Mega-analysis 
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  if(!dir.exists(megafolder)){
    dir.create(megafolder)
  }
  
  cat("\nMega-analysis using standardized expression values.")
  
  # standardized gene expression values per study (mean = 0, sd = 1)
  collapse = fread(paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_OutliersRemoved_allstudies.txt"))
  collapse$FACTOR_age = as.numeric(sub("\\+","",collapse$FACTOR_age))
  collapse = collapse[grep("CTL|AD",collapse$FACTOR_dx),]
  collapse$FACTOR_dx = factor(collapse$FACTOR_dx, levels = c("CTL","AD"))
  collapse = as.data.frame(collapse)
  ## eths to white/nonwhite
  collapse$FACTOR_ethnicity[collapse$FACTOR_ethnicity == ""] = NA
  foo = unique(collapse$FACTOR_ethnicity)
  foo = foo[!grepl("white",foo)]
  foo = foo[!is.na(foo)]
  collapse$FACTOR_ethnicity[collapse$FACTOR_ethnicity %in% foo] = "nonwhite"
  collapse$FACTOR_ethnicity = factor(collapse$FACTOR_ethnicity,levels=c("nonwhite","white"))
  
  study_id = as.character(unique(collapse$FACTOR_studyID))
  
  datScaled = collapse  ##it's already scaled, ain't it?
  # datScaled = list()
  # for(x in 1:length(study_id)){ ## TODO is this... redoing the scaling?
  #   
  #   datExpr = collapse[collapse$FACTOR_studyID %in% study_id[[x]], ]
  #   
  #   datGenes = datExpr[,!grepl("FACTOR_", colnames(datExpr))]
  #   datGenes = datGenes[,!colSums(is.na(datGenes)) > 0]
  #   
  #   # variance checker
  #   varcheck = apply(datGenes, 2, sd)
  #   lowvariance = varcheck[varcheck < 1e-02]
  #   
  #   datGenes = datGenes[,!colnames(datGenes) %in% names(lowvariance),]
  #   
  #   datGenes = scale(datGenes, scale = T, center = T)
  #   
  #   datScaled[[x]] = data.frame(datExpr[,grepl("FACTOR_", colnames(datExpr))], datGenes)
  #   
  # }
  # datScaled = ldply(datScaled)
  # ## TODO just ad/ctl plox
  # datScaled = datScaled[datScaled$FACTOR_dx %in% c("AD","CTL"),]
  # dim(datScaled)
  
  # SVA with z-standardized data set
  exprs = datScaled[,!colSums(is.na(datScaled)) > 0]
  exprs = exprs[,!grepl("FACTOR_", colnames(exprs))]
  dim(exprs)
  
  # Identify SVA
  require(sva)
  mod = model.matrix(~FACTOR_dx + FACTOR_studyID, data = datScaled) # model with known factors and covariates
  mod0 = model.matrix(~1,data=datScaled[as.numeric(rownames(mod)), ]) # intercept only model
  
  # # exprs = exprs[as.numeric(rownames(mod)), ]
  # exprs = exprs(ExpressionSet(t(exprs)))
  # 
  # n.sv = num.sv(exprs, mod ,method="leek") # Leek method for asymptotic modeling approach
  # message("\nNumber of Surrogate variables detected: ", n.sv)
  
  foo = exprs(ExpressionSet(t(exprs)))
  n.sv = max(num.sv(foo,mod,method='leek'),3) ## leek is sometimes producing no SVs.
  svobj = sva(foo,mod, n.sv = n.sv)
  svdf = data.frame(NULL)
  svdf = as.data.frame(svobj$sv)
  colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
  
  # ADD SURROGATE VARIABLES TO EXPRESSION SET
  datAll = data.frame(svdf, datScaled)
  
  datAll$FACTOR_dx = relevel(as.factor(datAll$FACTOR_dx),ref="CTL")
  
  ## linear regression models
  frm = " ~ FACTOR_dx + FACTOR_studyID + SV1"
  cat("\nBeginning linear regression models.\nFormula:",frm,"\n")
  
  require(lmerTest)
  
  genes = names(collapse)[names(collapse) %in% colnames(datAll)]
  genes = genes[-grep("FACTOR_",genes)]
  
  lme_results=list()
  lme_failure=character()
  for(x in 1:length(genes)){
    cat("\rLinear model:", x, "of", length(genes))
    # formula = formula(paste(genes[[x]], " ~ FACTOR_dx + SV1 + (1|FACTOR_studyID)"))
    formula = formula(paste(genes[[x]], frm))
    lmeFit = NULL
    lmeFit = try(lm(formula, data = datAll))
    
    # CHECK IF LM FIT, SKIP IF FAILED
    if(class(lmeFit) == 'try-error'){
      cat("  lm() failure for",genes[[x]],"\n")
      lme_failure = c(lme_failure,genes[[x]])
      next
    }
    
    Coefs = broom::tidy(summary(lmeFit)$coefficients)
    Coefs = Coefs[grepl("FACTOR_dxAD", Coefs$.rownames), ]
    Coefs = data.frame(GeneSymbol = genes[[x]], Coefs)
    lme_results[[x]] = Coefs
  }
  lme_results = ldply(lme_results)
  colnames(lme_results)[grepl("Pr..", colnames(lme_results))] = "P"
  lme_results$FDR = p.adjust(lme_results$P, 'fdr')
  lme_results$BONF = p.adjust(lme_results$P, 'bonferroni')
  lme_results = lme_results[order(lme_results$P, decreasing = F), ]
  
  write(lme_failure,paste(megafolder,"/",tissue,"_ADbr_Linear_Model_Failed_Genes.txt",sep=""))
  fwrite(lme_results,paste(megafolder,"/",tissue,"_ADbr_Linear_Model_Results.txt",sep=""))
  
  head(lme_results)
  
  psize = -log10(lme_results$P)
  psize = (psize - min(psize))/(max(psize) - min(psize)) + 0.5
  col = rep("lightgrey", nrow(lme_results))
  col[which(abs(lme_results$Estimate) > 0.2)] = "dodgerblue3"
  col[which(lme_results$FDR < .05)] = "orange"
  col[which(lme_results$BONF < .05)] = "firebrick3"
  
  lme_results$vLabel = NA
  for(foo in 1:length(lme_results$vLabel)){
    if(abs(lme_results$Estimate[foo]) > .7 | lme_results$BONF[foo] < 1e-7){
      lme_results$vLabel[foo] = as.character(lme_results$GeneSymbol[foo])
    }
  }
  lme_results$vLabel = gsub("[.]", "-", lme_results$vLabel)
  
  signCounts = table(sign(lme_results$Estimate))
  paste("n = ", signCounts, sep = "")
  
  png(paste0("./QCplots/",tissue,"_qcADbr_mega_volcano.png"),res=300,units="in",height = 5, width = 6)
  print(
    ggplot(lme_results, aes(x = Estimate, y = -log10(lme_results$P))) + 
    geom_point(size = psize, col = col) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), panel.border = element_rect(fill=NA,size = 1.5)) +
    geom_hline(yintercept = -log10(.05), col = 'orange', linetype = "dashed") + 
    xlab( expression(paste("log"[2]," fold-change in gene expression"))) +
    ylab( expression(paste("-log"[10],"(P-value)"))) + 
    geom_hline(aes(yintercept = -log10(.05/nrow(lme_results))), col = 'red', lwd = 0.25,linetype=1) +
    geom_vline(aes(xintercept = 0), col = 'black',lwd=0.3, linetype=2) +
    ggrepel::geom_text_repel(aes(label = lme_results$vLabel), size = 3, segment.size = 0.2, segment.colour = "grey", fontface = 3)
  )
  dev.off()
  
  
  # qq-plot
  
  # association p-values
  assoc = data.frame(P = lme_results$P, source='DGE Mega-analysis')
  # This is not the solution.
  # assoc = data.frame(P = lme_results$BONF, source='DGE Mega-analysis')
  observed = assoc$P
  chisq1 <- qchisq(1-observed, 1)
  observed <- sort(assoc$P)
  lobs <- -(log10(as.numeric(observed)))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  medianchi = median(chisq1,na.rm=T)
  lambdaout = medianchi/.454
  
  assoc = assoc[order(assoc$P, decreasing = F), ]
  assoc$lobs = lobs
  assoc$lexp = lexp
  
  ci = .95
  N = length(assoc$P)
  observed = -log10(sort(assoc$P))
  expected = -log10(1:N / N)
  clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
  cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  assoc$clower = clower
  assoc$cupper = cupper
  
  mclab = substitute(paste("Median ", chi^2, "=", MC, ", ", lambda , "=", LD), list(MC=format(medianchi, digits = 3), LD=format(lambdaout,digits=3)))
  
  qplot = ggplot(assoc, aes(x = lexp, y = lobs)) +
    geom_point(colour="black", fill= 'dodgerblue', shape = 21, size = 2, stroke = 0.3) + 
    geom_abline(slope = 1, intercept = 0, col = 'black', lwd = 0.3) +
    theme_classic() + 
    ylab(expression(paste("Observed -log"[10]," P-value"))) +
    xlab(expression(paste("Expected -log"[10]," P-value"))) +
    scale_fill_discrete(NULL) +
    geom_line(aes(x = lexp, y = clower), colour ='grey' ,lwd = 0.75) +
    geom_line(aes(x = lexp, y = cupper), colour = 'grey', lwd = 0.75) +
    geom_ribbon(aes(x = lexp, ymin = clower, ymax = cupper), fill="grey", alpha="0.2") +
    ggtitle(mclab)  +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
  
  # qplot
  
  png(paste0("./QCplots/",tissue,"_qcADbr_mega_qqplot.png"),res=300,units="in",height=5,width=5)
  print(qplot)
  dev.off()
  
  gc()
  
  ## ******************************************************
  ##   now compare!
  ## ******************************************************
  
  cat("\nPerforming correlation test between meta- and mega-analysis.")
  
  meta_results = fread(paste(metafolder,"/",tissue,"_freeze_qcAD_brain_meta_Nmin4_LeaveOneOut.txt",sep=""))
  mega_results = fread(paste(megafolder,"/",tissue,"_ADbr_Linear_Model_Results.txt",sep=""))
  
  meta_results$t.value = meta_results$Log2FC / meta_results$SE
  names(meta_results)[which(names(meta_results) == ".id")] = "GeneSymbol"
  
  rownames(meta_results) = meta_results$GeneSymbol
  rownames(mega_results) = mega_results$GeneSymbol
  
  meta_results = meta_results[order(meta_results$GeneSymbol),]
  mega_results = mega_results[order(mega_results$GeneSymbol),]
  
  sink(paste0("./",tissue,"_ADbr_metamega_trimmed_genes.txt"))
    cat("\nMeta genes not in mega:",meta_results$GeneSymbol[!(meta_results$GeneSymbol %in% mega_results$GeneSymbol)])
    cat("\nMega genes not in meta:",mega_results$GeneSymbol[!(mega_results$GeneSymbol %in% meta_results$GeneSymbol)])
  sink()
  
  meta_trim = meta_results[(meta_results$GeneSymbol %in% mega_results$GeneSymbol),]
  mega_trim = mega_results[(mega_results$GeneSymbol %in% meta_results$GeneSymbol),]
  
  cors = cor.test(meta_trim$t.value,mega_trim$t.value)
  print(cors)
  
  core = cor.test(meta_trim$Log2FC,mega_trim$Estimate)
  print(core)
  
  corp = cor.test(-log10(meta_trim$P),-log10(mega_trim$P))
  print(corp)
  
  sink(paste0(file="./",tissue,"_ADbr_metamega_correlations.txt"))
    print(cors)
    print(core)
    print(corp)
  sink()
}








