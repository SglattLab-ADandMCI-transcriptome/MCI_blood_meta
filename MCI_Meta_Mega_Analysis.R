setwd("~/PsychGENe/MCI_blood_meta/")

## these can be set to do AD instead.  tissues can be c(...)
analysislabel = "MCI"
caselabel = "MCI"
controllabel = "CTL"
covariateslist = c("FACTOR_dx", "FACTOR_sex", "FACTOR_age","FACTOR_race")
## this MCI study is only going to be on whole blood.
tissues = c("whole_blood")

## meta analysis for AD/MCI or whatever
## GCH w/JH and WB

##### TODO binomial test thing.  get a number. understand it.
## TODO deconvolution categories or pcs?

# load these packages (install if needed)
require(plyr)
require(ggplot2)
require(ggrepel)
require(data.table)
library(limma)
require(BRETIGEA)
require(AnnotationDbi)
require(sva)
require(lmerTest)
require(hgu133a.db)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(Homo.sapiens)
require(metafor)

datafolder = "./data_for_analysis/"
metafolder = "./meta_analysis"
megafolder = "./mega_analysis"
forestfolder = "./forestplots"

if(!dir.exists(metafolder)){
  dir.create(metafolder)
}

if(!dir.exists(forestfolder)){
  dir.create(forestfolder)
}

datfiles = list.files(datafolder,"_GeneExpression_allstudies.txt")
covfiles = list.files(datafolder,"_SampleFactors_allstudies.txt")
scalefiles = list.files(datafolder,"_ScaledWithFactors_OutliersRemoved_allstudies.txt")

dattis = sub("_GeneExpression_allstudies.txt","",datfiles)
covtis = sub("_SampleFactors_allstudies.txt","",covfiles)
scaletis = sub("_ScaledWithFactors_OutliersRemoved_allstudies.txt","",scalefiles)

if(all(dattis == covtis, covtis == scaletis)){
  cat("Tissue files found: ", dattis,"\nTissues to be analyzed: ",tissues,"\n")
} else {
  stop("Tissue files are... not right.")
}

tissue = tissues[1]
for(tissue in tissues){
  message("Beginning analysis: ",tissue)

  ## read in data files
  datExprTissue = fread(paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"),
                        data.table=F, stringsAsFactors = F)
  datExprCovs = fread(paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"),
                      data.table=F, stringsAsFactors = F)
  datExprCovs$FACTOR_age = as.numeric(sub("\\+","",datExprCovs$FACTOR_age))
  
    
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
  all(datExprCovs$FACTOR_sampleID == row.names(abundances))  ##TRUE
  
  col = factor(datExprCovs$FACTOR_studyID)
  g = ggplot(abpc3, aes(x = `cells PC1`, y = `cells PC2`,
                            col = col)) +
      scale_color_brewer('Study ID', palette = 'Spectral')+
      geom_point(size = 3) +
      theme_minimal() +
      theme(panel.border = element_rect(size=1,fill=NA),
            axis.title = element_text(size =12),
            axis.text=element_text(color='black',size=12))
  png(paste("./QCPlots/",tissue,"_PCA_deconvolution.png", sep =""), res=300,units="in",height = 6, width = 10)
  print(g)
  dev.off()
  ## After this QC we are going to actually use collapsed categories
  groups = c("B.cells|Plasma", "T.cells", "NK.cells", "Monocytes|Macro", "Dendritic", "Mast", "Eosino|Neutrophils")
  group_name = c("B.cells", "T.cells", "NK.cells", "Monocytes", "Dendritic.cells", "Mast.cells", "Granulocytes")
  ab_collapse = list()
  for(i in 1:length(groups)){
    include = grep(groups[i], names(abundances))
    foo = rowSums(abundances[,include])
    foo = data.frame(foo)
    names(foo) = group_name[i]
    ab_collapse[[i]] = foo
  }
  ab_collapse = data.frame(ab_collapse)
  names(ab_collapse) = group_name
  
  ## DGE per tissue per study, non-scaled data
  # genes = datExprTissue$SYMBOL
  # datExprTissue = data.frame(datExprCovs,t(datExprTissue[,-1]))
  genes = names(datExprTissue)[-1]
  datExprTissue = data.frame(datExprCovs,datExprTissue[,-1])
  names(datExprTissue) = c(names(datExprCovs),genes)
  datExprTissue = datExprTissue[grep(paste0(controllabel,"|",caselabel,"$"),datExprTissue$FACTOR_dx),]
  if(sum(datExprTissue$FACTOR_dx == caselabel)==0){
    message("No ",caselabel," cases found for tissue: ", tissue)
    next()
  }
  
  study_id = unique(datExprTissue$FACTOR_studyID)
  for(study in study_id){
    foo = sum(datExprTissue$FACTOR_dx[datExprTissue$FACTOR_studyID == study]==caselabel)
    if(foo==0){
      datExprTissue = datExprTissue[-which(datExprTissue$FACTOR_studyID == study),]
    }
  }
  
  print(table(datExprTissue$FACTOR_dx))
  print(table(datExprTissue$FACTOR_studyID, datExprTissue$FACTOR_dx))
  
  study_id = unique(datExprTissue$FACTOR_studyID)
  save_results = list()
  for( i in 1:length(study_id)){
    cat("\nStudy-wise differential expression analysis:",tissue,i,"~",study_id[[i]])
    expr.tmp = datExprTissue[datExprTissue$FACTOR_studyID %in% study_id[[i]],]
    
    N = table(expr.tmp$FACTOR_dx)
    Nca = N[[2]]
    Nco = N[[1]]
    cat("\n   Detected", Nca,"cases and",Nco,"controls")
    
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
      x = data.frame(x, y[,colnames(y) %in% names(wrong_class)]);
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
    PredListNames = paste(covariateslist, collapse= "|", sep="")
    predictors = x[,grep(PredListNames, colnames(x))]
    predictors = data.frame(predictors)
    
    # set level of disease groups 
    predictors$FACTOR_dx = as.factor(predictors$FACTOR_dx)
    predictors$FACTOR_dx = relevel(predictors$FACTOR_dx, ref = controllabel)
    
    counts = lapply(predictors,unique)
    count_class = unlist(lapply(counts,length))
    count_class = count_class[count_class <= 1,drop=F]
    
    if(length(count_class) > 0){predictors = predictors[,!(colnames(predictors) %in% names(count_class)),drop=F]}
    
    N = table(predictors$FACTOR_dx)
    Nca = N[[2]]
    Nco = N[[1]]
    cat("\n   Kept", Nca,"cases and",Nco,"controls\n")
    
    # remove genes with low variance
    var_filter = lapply(y, sd)
    low_var_filter = var_filter[is.na(var_filter) | var_filter < .001]
    
    if(length(low_var_filter) > 0){y = y[,!colnames(y) %in% names(low_var_filter)]}

    ##include blood cell abundances, without mast cells to avoid math singularity
    abindex = numeric()
    subjects = x$FACTOR_sampleID
    for(j in 1:length(subjects)){
      abindex[j] = which(abnames == subjects[j])
    }
    print("all(subjects == abnames[abindex])")
    print(all(subjects == abnames[abindex]))
    abtemp = ab_collapse[abindex,-grep("Mast",names(ab_collapse))]
    predictors = data.frame(predictors,abtemp)
    
    ## Surrogate variable analysis - default method 
    cat("Analysing for surrogate variables.\n")
    exprs = ExpressionSet(as.matrix(t(y)))
    mod = model.matrix(~ ., data= predictors) # model with known factors and covariates
    mod0 = model.matrix(~1,data=predictors) # intercept only model

    svobj = NULL
    foo = exprs(exprs)

    # bar = min(c(num.sv(foo,mod, method = 'be'), 3))
    # bar = num.sv(foo,mod,method='be')
    bar = num.sv(foo,mod,method='leek')  ##conservative approach
    cat("#SVs = ",bar,"\n")
    svobj = sva(foo,mod, n.sv = bar)
    svdf = as.data.frame(svobj$sv)

    ## include only the first 3
    # baz = min(bar,3)
    ## or include all
    baz = bar
    if(ncol(svdf) > 0){
      svdf = svdf[,1:baz,drop=F]
      colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
      predictors = data.frame(predictors, svdf)
    }

    # final design matrix for differential expression analysis
    # predictors$FACTOR_dx = as.factor(predictors$FACTOR_dx)
    # predictors$FACTOR_dx = relevel(predictors$FACTOR_dx, ref = controllabel)
    # design = model.matrix( ~ -1 + ., predictors)
    design = model.matrix( ~ ., predictors)
    design = design[,!grepl(controllabel, colnames(design))]
    
    # fit the linear model (arcsinh expression as response variable)
    lmFit = lm(as.matrix(y) ~ -1 + design)
    
    tidy_table = broom::tidy(lmFit)
    names(tidy_table)[1] = "GeneSymbol"
    
    big_table = data.frame(studyID=study_id[i],
                           tidy_table,
                           N_cases = Nca,
                           N_controls = Nco,
                           N = nrow(x)
    )
    save_results[[i]]  = big_table
    names(save_results)[[i]] = study_id[[i]]
  }
    
  ## merge study-wise sum stats into table
  cat("\nMerging sum stats.\n")
  mergestats = ldply(save_results)
  names(mergestats)[names(mergestats) %in% ".id"] = "studyID"
  mergestats$term = gsub("designFACTOR_|design", "", mergestats$term)
  
  cat("Studies and covariates:\n")
  print(unique(mergestats$studyID))
  print(unique(mergestats$term))
  
  genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  genes = data.frame(genes)
  genes$SYMBOL = select(org.Hs.eg.db,keys=as.character(genes$gene_id), keytype='ENTREZID', columns="SYMBOL")$SYMBOL
  genes$SYMBOL = gsub("[-]", ".", genes$SYMBOL)
  
  mergestats = mergestats[mergestats$GeneSymbol %in% genes$SYMBOL, ]
  
  fwrite(mergestats,
         file = paste(metafolder, "/",tissue,"_",analysislabel,"_perStudyDGEsummary.csv",sep=""),
         quote = F, row.names= F,sep=",")
  
  
  graph_df = mergestats
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
  
  png(file = paste("./QCplots/",tissue,"_qcsva_",analysislabel,"_violin.png", sep = ""),
      res = 300, units = "in", height = 8, width = 11.5)
  print(g)
  dev.off()
  
  
  agg_stats = data.table(combn)
  agg_stats = agg_stats[,list(LOGP = mean(-log10(p.value)), SD = sd(-log10(p.value)), K = length(p.value)),by=list(term, studyID)]
  agg_stats = agg_stats[!grepl("Intercept", agg_stats$term)]
  agg_stats$CI_LOW = agg_stats$LOGP - (1.96*agg_stats$SD/sqrt(agg_stats$K))
  agg_stats$CI_HIGH = (1.96*agg_stats$SD/sqrt(agg_stats$K)) + agg_stats$LOGP
  
  png(paste0("./QCplots/",tissue,"_",analysislabel,"_meanLOGP_term_DGE.png"),res=300,units="in",height=10,width=10)
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
  cat("\nGenewise random effect meta-analysis.\n")
  ## perform meta-analysis of genes (REML model)
  
  genes = unique(mergestats$GeneSymbol)
  genes = genes[!grepl("FACTOR_", genes)]
  
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
    
    ## must have at least 3 studies contributing
    if(nrow(sub_combn) < 3){
      cat("\rNot enough studies for",i,genes[[i]],"\n")
      next
    }
    if(nrow(sub) < 3){
      cat("\rNot enough studies for",i,genes[[i]],"\n")
      next
    }
    
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
      cat("\rGiving up on",i,genes[[i]],"\n")
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
                     arcsinh = res$beta, 
                     SE = res$se, 
                     # BetaStd = res.std$beta,
                     # SEstd = res.std$se,
                     P = res$pval,
                     HetQ = res$QE,
                     HetP = res$QEp,
                     study_kept = paste(sub_combn$studyID, collapse =  ", ", sep = ""))
    
    # ##reduce Nstudy RE both ANM batches
    # if(length(grep("AddNeuroMed",sub_combn$studyID)) > 1){
    #   res$Nstudy = res$Nstudy - 1
    # }
    
    res_save[[i]] = res
    names(res_save)[i] = genes[i]
  }
  
  foo = llply(res_save,class)
  foo = unlist(foo)
  res_df = ldply(res_save[which(foo == "data.frame")])

  min(res_df$P)
  
  fwrite(data.table(res_df), 
         sep = "\t",
         file = paste(metafolder,"/",tissue,"_",analysislabel,"_meta.txt", sep=""),
         quote  = F, row.names = F)
  write(names(res_save)[which(foo != "data.frame")],
        file = paste(metafolder,"/",tissue,"_",analysislabel,"_meta_NULLmetaforRMA.txt", sep=""))
  
  ## save the loo
  loo_df = ldply(loo_save)
  fwrite(loo_df,file = paste0(metafolder,"/",tissue,"_loosave.txt"))
  
  # res_df = fread("./meta_analysis/whole_blood_MCI_meta.txt",data.table=F)
  # ab = res_df[res_df$GeneSymbol %in% gene_filter$GeneSymbol, ]
  
  ## at least 3 studies are contributing to the gene
  res_df = res_df[res_df$Nstudy >= 3, ]
  res_df = res_df[order(res_df$P,decreasing=F),]
  res_df$FDR = p.adjust(res_df$P, "fdr")
  res_df$BONF = p.adjust(res_df$P, "bonferroni")
  res_df$GeneSymbol = gsub("[.]", "-", res_df$GeneSymbol)
  
  dim(res_df)
  
  head(res_df)
  
  # Location of genes 
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
         file = paste(metafolder,"/",tissue,"_",analysislabel,"_meta_mapped.txt",sep=""), quote  = F, row.names = F)
  
  
  # qq-plot
  cat("\nConstructing qq plot.")
  # association p-values
  assoc = data.frame(P = res_df$P, source='DGE Meta-analysis')
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
    geom_ribbon(aes(x = lexp, ymin = clower, ymax = cupper), fill="grey", alpha=0.2) +
    ggtitle(mclab)  +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
  
  # qplot
  png(paste0("./QCplots/",tissue,"_qc",analysislabel,"_meta_qqplot.png"),res=300,units="in",height=5,width=5)
  print(qplot)
  dev.off()
  
  
  ## Manhattan plot
  cat("\nConstructing manhattan plot.")
  trim = res_df[!is.na(res_df$LOC),]
  trim = trim[grepl("[:]", trim$LOC), ]
  locs = strsplit(trim$LOC, "[:]")
  chr = lapply(locs, function(x) x[[1]])
  bp = lapply(locs, function(x) x[[2]])
  pos = data.frame(CHR=unlist(chr),BP=unlist(bp))
  pos$GeneSymbol = trim$GeneSymbol
  pos$BONF = trim$BONF
  pos$FDR = trim$FDR
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
  sub$SYMBOL = ifelse(sub$FDR < .05, sub$GeneSymbol, NA)
  
  png(paste0("./QCplots/",tissue,"_qc",analysislabel,"_ggplot_manhattan.png"), res = 300, units = "in", height = 5, width = 9)
  # manhattan plot
  # the line is BONF and labels are FDR q<.05
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
  
  png(paste0("./QCplots/",tissue,"_qc",analysislabel,"_ggplot_manhattan_q.png"), res = 300, units = "in", height = 5, width = 9)
  # manhattan plot with q values
  # line and labels are FDR q<.05
  print(
    try(
      ggplot(sub, aes(x = sub$pos, y = -log10(sub$FDR))) + 
        geom_point(size = 0.7, col = as.character(sub$col)) + 
        ylim(min = 0, max = 1.2*max(-log10(sub$FDR))) + 
        xlab("Genomic coordinate") + 
        theme_classic() +
        ylab(expression(paste("-log"[10],"(Q-value)"))) +
        scale_x_continuous(name="Genomic coordinate", breaks=median_pos$V1, labels=median_pos$.id) +
        geom_hline(yintercept = -log10(.05), col = "black", lwd = 0.5, linetype = "dashed") +
        geom_text_repel(aes(label = sub$SYMBOL), fontface = 'italic', size = 3, col = "black") +
        theme(axis.title = element_text(size = 12), axis.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))
      # end plot
    ))
  dev.off()
  

  # volcano plot
  cat("\nConstructing volcano plot.")
  res_df$FDR = p.adjust(res_df$P, "fdr")
  res_df$BONF = p.adjust(res_df$P, "bonferroni")
  psize = -log10(res_df$P)
  psize = (psize - min(psize))/(max(psize) - min(psize)) + 0.5
  col = rep("lightgrey", nrow(res_df))
  # col[which(abs(res_df$arcsinh) > 0.1)] = "dodgerblue3"
  col[which(res_df$FDR < .05)] = "orange"
  col[which(res_df$BONF < .05)] = "firebrick3"
  
  ## decide what to label and save as significant
  res_df$vLabel = ifelse(res_df$P < .9, res_df$GeneSymbol, NA)
  res_df$vLabel = gsub("[.]", "-", res_df$vLabel)
  top_df = res_df[!is.na(res_df$vLabel),]
  fwrite(top_df,paste(metafolder,"/",tissue,"_",analysislabel,"_meta_significant_arcsinh_and_pvals.csv", sep=""))
  
  ## top pval only, to reduce clutter if needed
  # topD = order(abs(res_df$arcsinh),decreasing=T)
  # topP = order(res_df$FDR)
  # keepLabels = unique(c(topD[1:20],topP[1:30]))
  keepLabels = which(res_df$FDR < .05)
  dropLabels = c(1:nrow(res_df))
  dropLabels = dropLabels[-which(dropLabels %in% keepLabels)]
  res_df$vLabel[dropLabels] = NA
  
  signCounts = table(sign(res_df$arcsinh))
  paste("n = ", signCounts, sep = "")
  
  png(paste0("./QCplots/",tissue,"_qc",analysislabel,"_meta-volcano.png"),res=300,units="in",height = 5, width = 6)
  print(
    ggplot(res_df, aes(x = arcsinh, y = -log10(P))) + 
    geom_point(size = psize, col = col) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), panel.border = element_rect(fill=NA,size = 1.5)) +
    geom_hline(yintercept = -log10(.05), col = 'orange', linetype = "dashed") + 
    xlab( expression(paste("arcsinh score change in gene expression"))) +
    ylab( expression(paste("-log"[10],"(P-value)"))) + 
    geom_hline(aes(yintercept = -log10(.05/nrow(res_df))), col = 'red', lwd = 0.25,linetype=1) +
    geom_vline(aes(xintercept = 0), col = 'black',lwd=0.3, linetype=2) +
    ggrepel::geom_text_repel(aes(label = res_df$vLabel), size = 3, segment.size = 0.2, segment.colour = "grey", fontface = 3)
  )
  dev.off()
  

  
  ## Forest plots
  cat("\nConstructing meta forest plots.")
  pdf(paste(forestfolder,"/",tissue,"_FORESTPLOT_all.pdf",sep=""))
  
  res_filter = res_df[order(res_df$P,decreasing=F), ]
  topgenes = unique(res_filter$GeneSymbol[res_filter$FDR < .1])
  
  mergestats$GeneSymbol = gsub("[.]", "-", mergestats$GeneSymbol)
  
  if(length(topgenes)>1) for( i in 1:length(topgenes)){
    # subset meta-analysis results
    meta = res_filter[res_filter$GeneSymbol %in% topgenes[[i]], ]
    keep_studies = unlist(strsplit(as.character(meta$study_kept), ", "))
    
    meta = data.frame(arcsinh = meta$arcsinh,
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
      
      indiv = data.frame(arcsinh = sub$estimate, 
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
    results$CI_LOW = results$arcsinh - (1.96 * results$SE)
    results$CI_HIGH = results$arcsinh  + (1.96 * results$SE)
    
    locus = res_df$LOC[res_df$GeneSymbol %in% topgenes[[i]]]
    
    g = ggplot(results, aes(x = arcsinh, y = (studyID))) + 
      geom_point(col = results$col) + 
      theme_bw() +
      xlim(min(floor(results$CI_LOW)), max(ceiling(results$CI_HIGH))) +
      xlab(expression(paste("arcsinh score change (95% CI)"))) +
      ylab(NULL) +
      ggtitle(paste0(topgenes[[i]]," - Contribution")) +
      geom_vline(aes(xintercept = 0),linetype="dashed", col = "grey") +
      facet_grid(Source~., scales = "free_y", space = "free_y") + 
      geom_errorbarh(xmin = results$CI_LOW, xmax = results$CI_HIGH, col = results$col, height = NA) + 
      theme( strip.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
      geom_text(aes(vjust = -1, label = results$P_label))
    
    png(paste("./forestplots/RANK_",i,"_FORESTPLOT_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
    print(g)
    dev.off()
    print(g)
    
  } else print("No top genes.")
  dev.off()
  
  cat("\nConstructing LOO forest plots.")
  pdf(paste(forestfolder,"/",tissue,"_LOO_FORESTPLOT_all.pdf",sep=""))
  
  loo_df$GeneSymbol = gsub("[.]", "-", loo_df$GeneSymbol)
  
  foo = which(loo_df$GeneSymbol %in% topgenes)
  loo_filter = loo_df[foo, ]

  if(length(topgenes)>1) for( i in 1:length(topgenes)){
    # subset meta-analysis results
    meta = res_filter[res_filter$GeneSymbol %in% topgenes[[i]], ]
    keep_studies = unlist(strsplit(as.character(meta$study_kept), ", "))
    
    meta = data.frame(arcsinh = meta$arcsinh,
                      P = meta$P, SE = meta$SE, 
                      studyID = "Pooled result", 
                      Ncase = meta$Ncase, 
                      Ncontrol = meta$Ncontrol,
                      N = meta$Nsample,
                      col = "blue", 
                      Source = "Meta-analysis")
    
    
    sub = loo_filter[loo_filter$GeneSymbol %in% topgenes[[i]], ]
    sub = sub[sub$studyID %in% keep_studies, ] # retain correct studies for forest plot

    if(nrow(sub) > 0){
      indiv = data.frame(arcsinh = sub$estimate, 
                         P = sub$pval,
                         SE = sub$se, 
                         studyID = sub$studyID,
                         # Ncase = sub$N_cases, 
                         # Ncontrol = sub$N_controls,
                         col = "red", 
                         Source = "Studies");
      
      results = ldply(list(indiv,meta))
    }else{
      results = meta
    }
    
    results$P_label = paste("P = ", format(scientific=T,digits=3,results$P), sep = "")
    results$CI_LOW = results$arcsinh - (1.96 * results$SE)
    results$CI_HIGH = results$arcsinh  + (1.96 * results$SE)
    
    locus = res_df$LOC[res_df$GeneSymbol %in% topgenes[[i]]]
    
    g = ggplot(results, aes(x = arcsinh, y = (studyID))) + 
      geom_point(col = results$col) + 
      theme_bw() +
      xlim(min(floor(c(results$CI_LOW,-1))), max(ceiling(c(results$CI_HIGH,1)))) +
      xlab(expression(paste("arcsinh score change (95% CI)"))) +
      ylab(NULL) +
      ggtitle(paste0(topgenes[[i]]," - Leave One Out")) +
      geom_vline(aes(xintercept = 0),linetype="dashed", col = "grey") +
      facet_grid(Source~., scales = "free_y", space = "free_y") + 
      geom_errorbarh(xmin = results$CI_LOW, xmax = results$CI_HIGH, col = results$col, height = NA) + 
      theme( strip.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
      geom_text(aes(vjust = -1, label = results$P_label))
    
    png(paste("./forestplots/RANK_",i,"_LOO_FORESTPLOT_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
    print(g)
    dev.off()
    print(g)
    
  } else print("No top genes.")
  dev.off()
  
  
  
  ## Directionality gene score
  cat("\nConstructing directionality gene score and sign test plot.")
  res_df = res_df[order(res_df$P, decreasing = F), ]  ##order by diffex Pval
  
  majority = max(table(sign(res_df$arcsinh)))
  broom::tidy(binom.test(x = majority, n = nrow(res_df), p = 0.5, alternative = 'two.sided'))
  
  direction = sign(res_df$arcsinh)
  dirAll = sign(sum(sign(res_df$arcsinh)))
  
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
  binomStats_df = binomStats_df[order(binomStats_df$p.value,decreasing=F),] ## now sort by binomial test p value
  binomStats_df$label = NA
  binomStats_df$label[[1]] = paste("p = ",format(binomStats_df$p.value[[1]],digits=2,scientific=T),sep="")
  binomStats_df = binomStats_df[order(binomStats_df$rank,decreasing=F),]
  binomStats_df$bin = dplyr::ntile(binomStats_df$parameter, 1)
  cat("\n")
  
  png(paste("./QCplots/",tissue,"_qc",analysislabel,"_meta_signtest.png",sep=""), res = 300, units = "in", height = 5, width = 7.5)
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
  
  fwrite(binomStats_df,paste0(metafolder,"/",tissue,"_",analysislabel,"_meta_binomStats.csv"),row.names=T)
  
  gc()
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  # Mega-analysis 
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  
  ## no mega-analysis after consulting with SG and JH
}








