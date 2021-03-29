setwd("~/PsychGENe/MCI_blood_meta/")

## New forest plots and other information from meta analysis
## GCH

require(data.table)
require(biomaRt)
require(plyr)
require(ggplot2)


loo = fread("./meta_analysis/whole_blood_loosave.txt", data.table=F)
metar = fread("./meta_analysis/whole_blood_MCI_meta_significant_arcsinh_and_pvals.csv", data.table=F)


topgenes = unique(metar$GeneSymbol[metar$FDR < .05])
toploo = loo[which(loo$GeneSymbol %in% topgenes),]

## order toploo
foo = list()
for(i in 1:length(topgenes)){
  foo[[i]] = which(toploo$GeneSymbol == topgenes[i])
}
foo = unlist(foo)
all(1:length(foo) %in% foo) #TRUE
toploo = toploo[foo,]

## "pval" no longer sig
abolished = toploo[which(toploo$pval > .05),]
fwrite(abolished, file="./meta_analysis/whole_blood_loo_abolished.csv")

## direction change check
direction = sign(metar$arcsinh[metar$FDR < .05])
direction = data.frame(GeneSymbol = topgenes,direction)
dirchange = logical()
for(i in 1:nrow(toploo)){
  foo = sign(toploo$estimate[i])
  dirchange[i] = !(foo == direction$direction[direction$GeneSymbol == toploo$GeneSymbol[i]])
}
any(dirchange) #FALSE


## annotation for genes
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
atts = listAttributes(mart)
descriptions = getBM(attributes = c("hgnc_symbol","description"),
              filters = "hgnc_symbol",
              values = topgenes,
              mart = mart)
## order descriptions
foo = list()
for(i in 1:length(topgenes)){
  foo[[i]] = which(descriptions$hgnc_symbol == topgenes[i])
}
foo = unlist(foo)
all(1:length(foo) %in% foo) #TRUE
descriptions = descriptions[foo,]
fwrite(descriptions, file="./meta_analysis/whole_blood_top_descriptions.csv")




## Forest plots
## TODO reorder so largest study is close to result, sort by study size
## TODO dot size proportional to weight in meta analysis, with weight written on right hand side
cat("\nConstructing meta forest plots.")

res_filter = metar[order(metar$P,decreasing=F), ]
topgenes = unique(res_filter$GeneSymbol[res_filter$FDR < .05])

mergestats = fread("./meta_analysis/whole_blood_MCI_perStudyDGEsummary.csv",data.table=F)
mergestats$GeneSymbol = gsub("[.]", "-", mergestats$GeneSymbol)
mergestats = mergestats[order(mergestats$N),]

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
                       N = sub$N,
                       col = "red", 
                       Source = "Studies");
    
    meta$N = sum(indiv$N)
    results = ldply(list(indiv,meta))
  } else{results = meta}
  
  results$P_label = paste("P = ", format(scientific=T,digits=3,results$P), sep = "")
  results$CI_LOW = results$arcsinh - (1.96 * results$SE)
  results$CI_HIGH = results$arcsinh  + (1.96 * results$SE)
  results$studyID = factor(results$studyID, levels = results$studyID)
  
  locus = metar$LOC[metar$GeneSymbol %in% topgenes[[i]]]
  
  g = ggplot(results, aes(x = arcsinh, y = studyID)) + 
    geom_point(col = results$col, shape=18, size = 1.2*(.7+results$N/max(results$N))) + 
    theme_bw() +
    xlim(min(floor(results$CI_LOW)), max(ceiling(results$CI_HIGH))) +
    xlab(expression(paste("arcsinh score change (95% CI)"))) +
    ylab(results$studyID) +
    ggtitle(paste0(topgenes[[i]]," - Contribution")) +
    geom_vline(aes(xintercept = 0),linetype="dashed", col = "grey") +
    facet_grid(Source~., scales = "free_y", space = "free_y") + 
    geom_errorbarh(xmin = results$CI_LOW, xmax = results$CI_HIGH, col = results$col, height = NA) + 
    theme( strip.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
    geom_text(aes(x=1, label = results$N),size=2) +
    geom_text(aes(vjust = -1, label = results$P_label))
  
  png(paste("./forestplots/",i,"_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
  print(g)
  dev.off()
  
} else print("No top genes.")


cat("\nConstructing LOO forest plots.")
toploo$GeneSymbol = gsub("[.]", "-", toploo$GeneSymbol)

foo = which(toploo$GeneSymbol %in% topgenes)
loo_filter = toploo[foo, ]

Ns = mergestats[,c(1,11)]
Ns = unique(Ns)

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
    for(q in 1:nrow(sub)){
      sub$N[q] = Ns$N[which(Ns$studyID == sub$studyID[q])]
    }
    indiv = data.frame(arcsinh = sub$estimate, 
                       P = sub$pval,
                       SE = sub$se, 
                       studyID = sub$studyID,
                       N = sub$N,
                       col = "red", 
                       Source = "Studies");
    
    meta$N = sum(indiv$N)
    results = ldply(list(indiv,meta))
  }else{
    results = meta
  }
  
  results$P_label = paste("P = ", format(scientific=T,digits=3,results$P), sep = "")
  results$CI_LOW = results$arcsinh - (1.96 * results$SE)
  results$CI_HIGH = results$arcsinh  + (1.96 * results$SE)
  results$studyID = factor(results$studyID, levels = results$studyID)
  
  locus = metar$LOC[metar$GeneSymbol %in% topgenes[[i]]]
  
  g = ggplot(results, aes(x = arcsinh, y = (studyID))) + 
    geom_point(col = results$col, shape=18, size = 1.2*(.7+results$N/max(results$N))) + 
    theme_bw() +
    xlim(min(floor(c(results$CI_LOW,-1))), max(ceiling(c(results$CI_HIGH,1)))) +
    xlab(expression(paste("arcsinh score change (95% CI)"))) +
    ylab(NULL) +
    ggtitle(paste0(topgenes[[i]]," - Leave One Out")) +
    geom_vline(aes(xintercept = 0),linetype="dashed", col = "grey") +
    facet_grid(Source~., scales = "free_y", space = "free_y") + 
    geom_errorbarh(xmin = results$CI_LOW, xmax = results$CI_HIGH, col = results$col, height = NA) + 
    theme( strip.text.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) + 
    geom_text(aes(x=1, label = results$N),size=2) +
    geom_text(aes(vjust = -1, label = results$P_label))
  
  png(paste("./forestplots/",i,"_LOO_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
  print(g)
  dev.off()
} else print("No top genes.")
