setwd("~/PsychGENe/brain/")

# Merge expression data (with some other housekeeping)
# GCH w/JH and WB

# Update gene symbols for AD gene expression data sets
qcFolder = "./QCplots/"

# tissues = c("hippocampus","frontal_cortex","temporal_cortex","cerebellum","whole_blood")
tissues = "whole_blood"

require(data.table)
require(plyr)
require(org.Hs.eg.db)
require(AnnotationDbi)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(ggplot2)

cat("Compiling gene symbols.\n")
# https://www.gencodegenes.org/releases/22.html
# get list of gene symbols
genes = fread("./references/gencode.v22.annotation.gtf.gz",header=F,stringsAsFactors=FALSE)
genes = data.frame(genes)
genes = genes[genes$V3 %in% "gene", ]
split = strsplit(genes$V9, "[; ]")
split = lapply(split, function(x) x[[11]])
symbols = unlist(split)

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

genes$SYMBOL <- unlist(lapply(genes$V9, extract_attributes, "gene_name"))
# names(genes)[names(genes) %in% "symbol"] = "SYMBOL"
conv = select(org.Hs.eg.db, keys=as.character(genes$SYMBOL), keytype="SYMBOL", columns="ENTREZID")
genes = merge(genes, conv, by='SYMBOL')
genes = genes[!is.na(genes$ENTREZID), ]
genes$width = abs(genes$V4 - genes$V5)
genes = genes[order(genes$width, decreasing = T), ]
genes = genes[!duplicated(genes$SYMBOL), ]

### Begin updating tables
cat("\nReading file, scaling, separating tissues, and updating gene symbols.\n")
rawall = fread("normalized_data/ADMCI_merged.txt", header = T, stringsAsFactors = F, data.table = F)

## race cleanup
foo = rawall$FACTOR_race
foo[foo == "Caucasian"] = "white"
foo[foo == "caucasian"] = "white"
foo[foo == "White"] = "white"
foo[foo == "Japanese"] = "asian"
foo[foo == "other"] = "unknown"
foo = gsub(" not","",foo)
rawall$FACTOR_race = factor(foo)
print(table(rawall$FACTOR_race))

## sex cleanup
foo = rawall$FACTOR_sex
foo[grep("NA",foo)] = NA
f = grep("female",foo)
foo[grep("male",foo)] = "male"
foo[grep("M",foo)] = "male"
foo[grep("F",foo)] = "female"
foo[f] = "female"
rawall$FACTOR_sex = factor(foo)
print(table(rawall$FACTOR_sex))

## tissue cleanup
foo = rawall$FACTOR_tissue
foo[foo == "hippocampus - CA1"] = "hippocampus"
foo[foo == "HIP"] = "hippocampus"
foo[foo == "PFC"] = "frontal_cortex"
foo[foo == "FCX"] = "frontal_cortex"
foo[foo == "CR"] = "cerebellum"
foo[foo == "CBE"] = "cerebellum"
foo[foo == "CER"] = "cerebellum"
foo[foo == "TCX"] = "temporal_cortex"
foo[foo == "Whole blood"] = "whole_blood"
foo[foo == "whole blood"] = "whole_blood"
rawall$FACTOR_tissue = factor(foo)
print(table(rawall$FACTOR_tissue))

names(rawall)[which(names(rawall) == "Sample_ID")] = "FACTOR_sampleID"
for(tissue in tissues){
  message("TISSUE: ",tissue)
  rawdat = rawall[which(rawall$FACTOR_tissue==tissue),]
  study_id = unique(rawdat$FACTOR_studyID)
  df_list = list()
  cat("Updating gene symbols.\n")
  for(study in study_id){
    # read in the file and match genes
    cat(study,"\n")
    readIn = rawdat[which(rawdat$FACTOR_studyID == study),]
    datExpr = readIn[,!grepl("FACTOR_", colnames(readIn))]
    rownames(datExpr) = readIn$FACTOR_sampleID
    
    genes_in_data = colnames(datExpr)
    genes_in_data = gsub("[.]", "-", genes_in_data)
   
    nomatch = genes_in_data[!genes_in_data %in% genes$SYMBOL]
    nomatch_conv = select(org.Hs.eg.db, keys=as.character(nomatch), keytype='ALIAS', columns='ENTREZID')
    nomatch_conv$updated_hgnc = select(org.Hs.eg.db, keys=as.character(nomatch_conv$ENTREZID), keytype='ENTREZID', columns='SYMBOL')$SYMBOL
    nomatch_conv$mismatch = nomatch_conv$ALIAS == nomatch_conv$updated_hgnc
    noentrez = nomatch_conv[is.na(nomatch_conv$ENTREZID), ]
    
    datExpr = datExpr[,!colnames(datExpr) %in% noentrez$SYMBOL]
    datExpr = as.data.frame(t(datExpr))
    datExpr = data.frame(SYMBOL = rownames(datExpr), datExpr,check.names=F)
    names(datExpr) = c("SYMBOL",readIn$FACTOR_sampleID)
    
    df_list = c(df_list,list(datExpr))
  }
  cat("\nMerging.\n")
  datAll = Reduce(function(x,y) merge(x,y,by='SYMBOL',all = TRUE), df_list)
  dim(datAll)
  rm(df_list)
  
  # cat("\nWriting merged numeric table to file.")
  # if(!dir.exists("./data_for_analysis/")){ dir.create("./data_for_analysis/")}
  # fwrite(datAll, file=paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"),quote=F,row.names=F,sep="\t")
  # # datAll = fread("./data_for_analysis/GeneExpression_allstudies.txt")
  # 
  ### Bring in sample factors.
  cat("\nBringing in sample factors.")
  sf_list = list()
  for (study in study_id){
    cat("\nWorking on study:", study,"\n")
    readIn = rawdat[which(rawdat$FACTOR_studyID == study),]
    datExpr = readIn[,grepl("FACTOR_", colnames(readIn))]
    rownames(datExpr) = readIn$FACTOR_sampleID
    
    if(nrow(datExpr) < 1 ) next
    
    cat(colnames(datExpr),"\n")
    
    factor_keep = paste("FACTOR_", c("dx","age","race","sex","tissue","sampleID","studyID"), sep = "")
    
    # datExpr$FACTOR_studyID = strsplit(basename(files[[x]]),"_")[[1]][[2]]
    
    sf_list = c(sf_list,list(datExpr[,colnames(datExpr) %in% factor_keep]))
  
  }
  cat("\nMerging.\n")
  sfall = ldply(sf_list)
  
  head(sfall)
  cat("nrow(sfall) == (ncol(datAll)-1)     ",nrow(sfall) == (ncol(datAll)-1))
  # 
  # cat("\nWriting merged factor table to file.\n")
  # fwrite(sfall, file=paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"),quote=F,row.names=F,sep="\t")
  # 
  
  ## rebind sample factors with datAll (phenos + gene expression)
  cat("\nBinding factors with expression data.\n")
  # sfall$FACTOR_sampleID = gsub(" 1", ".1", sfall$FACTOR_sampleID)
  
  rownames(datAll)=datAll$SYMBOL
  datAll = datAll[,!colnames(datAll) %in% "SYMBOL"]
  datAll = as.data.frame(t(datAll))
  datAll = data.frame(FACTOR_sampleID = as.character(rownames(datAll)), datAll)
  
  datExpr0 = merge(sfall, datAll, by='FACTOR_sampleID')
  
  cat("dim(datExpr0)     ",dim(datExpr0))
  
  datExpr0 = datExpr0[order(datExpr0$FACTOR_studyID), ]
  
  geneExpr = datExpr0[,!grepl("FACTOR_", colnames(datExpr0))]
  geneExpr = geneExpr[,colSums(is.na(geneExpr))==0]
  
  # PCA analysis
  cat("\nPerforming PCA analysis.\n")
  pca = prcomp(geneExpr)

  # calculate variance explained by PCs
  sdev=(pca$sdev^2)/sum(pca$sdev^2)*100
  sdev=paste("PC", 1:length(sdev), " (", round(sdev,2),"%)", sep="")
  
  pca = data.frame(pca$x)
  pca$studyID = datExpr0$FACTOR_studyID
  
  png(paste(qcFolder,tissue,"_PCA_original_values_merged_genes.png", sep =""), res=300,units="in",height = 6, width = 10)
  print(
    ggplot(pca, aes(x = PC1, y = PC2, col = factor(studyID))) +
    scale_color_brewer('Study ID', palette = 'Spectral')+
    geom_point(size = 3) +
    theme_minimal() +
    theme(panel.border = element_rect(size=1,fill=NA), axis.title = element_text(size =12), axis.text=element_text(color='black',size=12)) +
    xlab(sdev[1]) +
    ylab(sdev[2])
  )
  dev.off()
  
  
  # boxplot with original values
  cat("\nOutputting boxplot with original values.")
  col = data.frame(FACTOR_studyID = unique(as.factor(datExpr0$FACTOR_studyID)))
  col$color = rainbow(nrow(col))
  
  factors = datExpr0[,grepl("FACTOR_", colnames(datExpr0))]
  factors$rank = 1:nrow(factors)
  factors = merge(col, factors,by='FACTOR_studyID')
  factors = factors[order(factors$rank, decreasing = F), ]
  
  legend_col = factors[!duplicated(factors$FACTOR_studyID), ]
  
  png(file = paste(qcFolder,tissue,"_BOXPLOT_original_values_merged_genes.png",sep=""), res =300 , units = "in", width = 8, height = 5)
  boxplot(t(geneExpr), outline=FALSE, 
          xlab="Samples",
          ylab=expression(paste("Original values")),
          col = factors$color, las = 1, xaxt='n',
          border = factors$color, las = 1, xaxt='n',
          main = paste(tissue,"Original Values, Merged"))
  legend("bottomleft", pch = 15, ncol = 2, cex=0.75, legend = legend_col$FACTOR_studyID, col = legend_col$color, bty="n")
  dev.off()
  
  # Outliers and Standardize gene expression across samples per study (expr mean = 0, var = 1)
  cat("\nRemoving outliers and scaling across studies to mean=0, var=1.\n")
  studyid = unique(datExpr0$FACTOR_studyID)
  
  keep = list()
  keepo = list()
  keepf = list()
  pca_outlier_list = list()
  for(x in 1:length(studyid)){
    
    cat("\n\nWorking on study:", studyid[[x]])
    dat = datExpr0[datExpr0$FACTOR_studyID %in% studyid[[x]], ]
    dat_sf = dat[,grepl("FACTOR_",colnames(dat))]
    dat = dat[,!grepl("FACTOR_", colnames(dat))]
    rownames(dat) = dat_sf$FACTOR_sampleID
    names(dat) = toupper(names(dat))
    
    # PCA outlier detection in each study
    
    cat("\nPerforming PCA outlier detection.")
    variance = apply(dat, 2, function(x) var(x))
    pca = prcomp(dat[,which(variance > 0.02)], center=T,scale=T)
    
    pca_var = (pca$sdev^2)/sum(pca$sdev^2)*100
    pca_labels = paste("PC", 1:length(pca_var), " (", round(pca_var, 2), "%)", sep= "")
    
    pca = data.frame(pca$x)
    
    pca$dx = dat_sf$FACTOR_dx
    
    png(file = paste(qcFolder,tissue,"_PCA_outlier_detection_",studyid[[x]],".png",sep=""), res =300 , units = "in", width = 8, height = 5)
    print(
      ggplot(pca, aes(x = PC1, y = PC2, fill = factor(dx))) +
      geom_point(size = 3, alpha = 0.5, pch = 21, col = 'black') +
      scale_fill_brewer('Diagnosis', palette = 'Set1') +
      theme_minimal() +
      geom_text(aes(label = rownames(pca)),nudge_y=-2, nudge_x=2, size = 2) +
      xlab(pca_labels[1]) +
      ylab(pca_labels[2]) +
      ggtitle(paste("PCA for:", studyid[[x]])) +
      theme(panel.border = element_rect(fill = NA, size = 1))
    )
    dev.off()
    
    pca_outlier = list()
    for(y in 1:3){  ## take out past 4 sd for first 3 PCs.  Per Glatt JAACAP 2012
      mean_x = mean(pca[,y])
      sd_filter = sd(pca[,y]) * 4  ## number of SDs
      pca_outlier[[y]] = abs(pca[,y]) > (mean_x + sd_filter)
    }
    pca_outlier = data.frame(do.call(cbind, pca_outlier))
    pca_outlier = which(rowSums(pca_outlier) > 0)
    
    if(length(pca_outlier) > 0) {
      cat("\nOutliers:", paste(rownames(dat)[pca_outlier], sep=","))
    } else {
      cat("\nOutliers: ~none~")
    }
    pca_outlier_list[[x]] = rownames(dat)[pca_outlier]
    names(pca_outlier_list)[x] = studyid[[x]]
    
    # remove outliers from dat and dat_sf
    if(length(pca_outlier)>0){
      dat = dat[-pca_outlier,]
      dat_sf = dat_sf[-pca_outlier,]
    }
    
    ##save nonscaled data
    keepo[[x]] = data.frame(FACTOR_sampleID = dat_sf$FACTOR_sampleID, dat)
    keepf[[x]] = data.frame(dat_sf)
    
    # scale expression
    cat("\nScaling.")
    dat = dat[,colSums(is.na(dat))/nrow(datAll) < 0.1]
    dat = scale(dat, center = T, scale = T)
    
    # PCA on scaled data.
    {
      cat("\nPerforming PCA on scaled study data.")
      variance = apply(dat, 2, function(x) var(x))
      pca = prcomp(dat[,which(variance > 0.02)], center=T,scale=T)
      
      pca_var = (pca$sdev^2)/sum(pca$sdev^2)*100
      pca_labels = paste("PC", 1:length(pca_var), " (", round(pca_var, 2), "%)", sep= "")
      
      pca = data.frame(pca$x)
      
      pca$dx = dat_sf$FACTOR_dx
      
      png(file = paste(qcFolder,tissue,"_PCA_after_scaling_",studyid[[x]],".png",sep=""), res =300 , units = "in", width = 8, height = 5)
      print(
        ggplot(pca, aes(x = PC1, y = PC2, fill = factor(dx))) +
          geom_point(size = 3, alpha = 0.5, pch = 21, col = 'black') +
          scale_fill_brewer('Diagnosis', palette = 'Set1') +
          theme_minimal() +
          geom_text(aes(label = rownames(pca)),nudge_y=-2, nudge_x=2, size = 2) +
          xlab(pca_labels[1]) +
          ylab(pca_labels[2]) +
          ggtitle(paste("PCA for:", studyid[[x]],tissue)) +
          theme(panel.border = element_rect(fill = NA, size = 1))
      )
      dev.off()
    }
    
    # Combine sample factors with gene expression
    dat = data.frame(dat_sf, dat)
    
    keep[[x]] = dat
  }
  
  # bind all scaled expression data together
  cat("\nBinding data.")
  scaledDf = ldply(keep)
  orgD = ldply(keepo)
  facD = ldply(keepf)
  rm(keep)
  rm(keepo)
  rm(keepf)
  
  cat("\nWriting merged factor table to file.\n")
  fwrite(facD, file=paste0("./data_for_analysis/",tissue,"_SampleFactors_allstudies.txt"),quote=F,row.names=F,sep="\t")
  
  cat("\nWriting merged numeric table to file.")
  if(!dir.exists("./data_for_analysis/")){ dir.create("./data_for_analysis/")}
  fwrite(orgD, file=paste0("./data_for_analysis/",tissue,"_GeneExpression_allstudies.txt"),quote=F,row.names=F,sep="\t")
  
  # notate outliers
  sink(paste0("./data_for_analysis/",tissue,"_PCA_removed_outliers.txt"))
  print(pca_outlier_list)
  sink()
  
  ## R is leaking memory and can't handle doing these plots.
  ## not doing mega-analysis anyway so we don't need the scaled/merged table
  ## TODO someday split up the studies or something, we don't need all these 0s and NAs
  # rm(datAll)
  # rm(datExpr0)
  # rm(geneExpr)
  # rm(orgD)
  # rm(rawall)
  # rm(rawdat)
  # rm(readIn)
  # rm(sf_list)
  # rm(split)
  # gc()
  # 
  # 
  # 
  # 
  # 
  # # boxplot of scaled expression
  # cat("\nGenerating final boxplots and writing data.")
  # scaledExpr = scaledDf[,!grepl("FACTOR_", colnames(scaledDf))]
  # # boxplot(t(scaledExpr), outline=F, col = as.factor(scaledDf$FACTOR_studyID))
  # col = rainbow(length(unique(scaledDf$FACTOR_studyID)))
  # col = col[as.factor(scaledDf$FACTOR_studyID)]
  # png(file = paste(qcFolder,tissue,"_BOXPLOT_scaled_values_merged_genes.png",sep=""), res =300 , units = "in", width = 8, height = 5)
  # boxplot(t(scaledExpr), outline=FALSE, 
  #         xlab="Samples",
  #         ylab=expression(paste("Scaled values")),
  #         col = col,
  #         border = col,
  #         main = paste(tissue,"Scaled Gene Expression Values"))
  # legend("bottomleft", pch = 15, ncol = 2, cex=0.75, legend = unique(scaledDf$FACTOR_studyID), col = unique(col), bty="n")
  # dev.off()
  # 
  # fwrite(scaledDf,paste0("./data_for_analysis/",tissue,"_ScaledWithFactors_OutliersRemoved_allstudies.txt"))
  # 
  # # PCA analysis
  # cat("\nPerforming PCA analysis on scaled and merged data.\n")
  # pca = prcomp(scaledExpr[,colSums(is.na(scaledExpr))==0])
  # 
  # # calculate variance explained by PCs
  # sdev=(pca$sdev^2)/sum(pca$sdev^2)*100
  # sdev=paste("PC", 1:length(sdev), " (", round(sdev,2),"%)", sep="")
  # 
  # pca = data.frame(pca$x)
  # pca$studyID = scaledDf$FACTOR_studyID
  # 
  # png(paste(qcFolder,tissue,"_PCA_scaled_values_merged_genes.png", sep =""), res=300,units="in",height = 6, width = 10)
  # print(
  #   ggplot(pca, aes(x = PC1, y = PC2, col = factor(studyID))) +
  #     scale_color_brewer('Study ID', palette = 'Spectral')+
  #     geom_point(size = 3) +
  #     theme_minimal() +
  #     theme(panel.border = element_rect(size=1,fill=NA), axis.title = element_text(size =12), axis.text=element_text(color='black',size=12)) +
  #     xlab(sdev[1]) +
  #     ylab(sdev[2])
  # )
  # dev.off()
}



