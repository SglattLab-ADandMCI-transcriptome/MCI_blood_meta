## i'm just moving the stuff from there to here now.  this doesn't actually work. yet.

# require(dtangle)
# require(dtangle.data)
## https://umich.app.box.com/v/dtangledatapkg


# deconvolution SECTION
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





svobj = NULL
foo = exprs(exprs)
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


