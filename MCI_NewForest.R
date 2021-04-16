setwd("~/PsychGENe/MCI_blood_meta/")

## New forest plots and other information from meta analysis
## GCH

## TODO GSER for top genes

require(data.table)
require(biomaRt)
require(plyr)
require(ggplot2)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(GO.db)
require(reactome.db)


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
bar = strsplit(descriptions[,2],"\\[")
bar = data.frame(bar)
bar = t(bar)
descriptions$description = bar[,1]
descriptions$sourcenotes = bar[,2]
fwrite(descriptions, file="./meta_analysis/whole_blood_top_descriptions.csv")




## Forest plots
cat("\nConstructing meta forest plots.")

res_filter = metar[order(metar$P,decreasing=F), ]
topgenes = unique(res_filter$GeneSymbol[res_filter$FDR < .05])

mergestats = fread("./meta_analysis/whole_blood_MCI_perStudyDGEsummary.csv",data.table=F)
mergestats$GeneSymbol = gsub("[.]", "-", mergestats$GeneSymbol)
mergestats = mergestats[order(mergestats$N),]

i=1
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
  results$percentlab = paste0(trunc(results$N / max(results$N) * 100),"%")
  results$percentlab[grep("100",results$percentlab)] = paste0("n=",max(results$N))
  
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
    geom_text(aes(x=1, label = results$percentlab),size=2) +
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

i=1
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
  results$percentlab = paste0(trunc(results$N / max(results$N) * 100),"%")
  results$percentlab[grep("100",results$percentlab)] = paste0("n=",max(results$N))
  
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
    geom_text(aes(x=1, label = results$percentlab),size=2) +
    geom_text(aes(vjust = -1, label = results$P_label))
  
  png(paste("./forestplots/",i,"_LOO_",topgenes[[i]],".png", sep = ""),res=300,units="in",height=5,width =5.5)
  print(g)
  dev.off()
} else print("No top genes.")




## Gene set enrichment analysis for interesting module eigengenes
## whole-cohort WGCNA
cat("\nPerforming gene set enrichment analysis.")

net = readRDS("./WGCNA_results/WGCNA.Rdata")
MEs = net$MEs
# svobj = readRDS("./WGCNA_results/ME_SVA_OBJ.Rdata")
datExpr = fread("./data_for_analysis/whole_blood_ScaledWithFactors_OutliersRemoved_allstudies.txt",data.table=F)
phenos = fread("./data_for_analysis/whole_blood_SampleFactors_allstudies.txt",data.table=F)
module = fread("./WGCNA_results/wgcna_module-membership.txt")

##Generate abundances categories
cat("Generating abundances categories\n")
abfiles = list.files("./deconvolution","abundances", full.names = T)
ablist = list()
for(file in 1:length(abfiles)){
  ablist[[file]] = fread(abfiles[file], data.table=F)
}
abundances = ldply(ablist)
abnames = abundances$V1
abundances = abundances[,-1]
row.names(abundances) = abnames
## use collapsed categories
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
ab_collapse = data.frame(row.names(ab_collapse),ab_collapse)
names(ab_collapse)[1] = "FACTOR_sampleID"

tempab = list()
for(i in 1:nrow(phenos)){
  tempab[[i]] = ab_collapse[which(ab_collapse$FACTOR_sampleID == phenos$FACTOR_sampleID[i]),]
}
ab_collapse = ldply(tempab)
phenos = data.frame(phenos,ab_collapse)

## Remove ROSMAP batch 3 due to low gene overlap with other studies
r3 = which(phenos$FACTOR_studyID=="ROSMAP3")
phenos = phenos[-r3,]
datExpr = datExpr[-r3,]

## Limit to Dx we care about
gooddx = grep("(MCI|CTL)$",phenos$FACTOR_dx)
phenos = phenos[gooddx,]
datExpr = datExpr[gooddx,]

datExpr = datExpr[,-grep("FACTOR_",names(datExpr))]
##remove any missingness
badgenes = which(colSums(is.na(datExpr)) > 1)
if(length(badgenes) > 0){
  datExpr = datExpr[,-badgenes]
}


# load gene names for hg19
genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes = as.data.frame(genes)
genename = AnnotationDbi::select(org.Hs.eg.db, keys = genes$gene_id, keytype="ENTREZID", columns=c("GENENAME","SYMBOL", "GO"))

# pull out gene annotations from reactome and gene ontology databases
reactome = AnnotationDbi::select(reactome.db, keys=as.character(unique(genename$ENTREZID)), keytype='ENTREZID', columns=c('PATHNAME', 'REACTOMEID'))
reactome$PATHNAME = gsub("Homo sapiens:", "REACTOME:", reactome$PATHNAME)
conv =  AnnotationDbi::select(org.Hs.eg.db, keys=as.character(unique(genename$ENTREZID)), keytype='ENTREZID', columns='SYMBOL')
reactome = merge(conv, reactome, by='ENTREZID')
colnames(reactome)[3] = 'PATHWAY'
reactome = reactome[,colnames(reactome) %in% c("SYMBOL", "PATHWAY")]

cat("GO heirarchy: Biological Processes, Cellular Component, Molecular Function")
go_terms = AnnotationDbi::select(GO.db, keys=as.character(genename$GO), keytype="GOID", columns=c("TERM", "ONTOLOGY"))
go_terms = go_terms[!duplicated(go_terms), ]
table(go_terms$ONTOLOGY)

names(genename)[names(genename) %in% "GO"] = "GOID"
go_terms = merge(go_terms, genename[,colnames(genename) %in% c("GOID", "SYMBOL")], by="GOID")
go_terms = go_terms[!duplicated(go_terms), ]
go_terms$PATHWAY = paste(go_terms$GOID,":", gsub("[,| ]", "_", go_terms$TERM), sep = "")
go_terms = go_terms[!is.na(go_terms$GOID), ]

go_terms = go_terms[,colnames(go_terms) %in% c("SYMBOL", "PATHWAY")]
go_terms = ldply(list(go_terms, reactome))
go_terms$in_data = ifelse(go_terms$SYMBOL %in% colnames(datExpr), "in", "out")

count_int = table(go_terms$PATHWAY, go_terms$in_data)
perc_overlap = count_int[,1]/rowSums(count_int)
perc_overlap = perc_overlap[perc_overlap > 0.5]

# Trim out some large umbrella sets
countz = table(factor(go_terms$PATHWAY))
lattice::histogram(as.numeric(countz))
plot(density(as.numeric(countz)))
paths = names(countz)
toobigcutoff = 100
toobig = paths[countz > toobigcutoff]
foo = go_terms[!(go_terms$PATHWAY %in% toobig),]
lattice::histogram(as.numeric(table(factor(foo$PATHWAY))))
sink("./WGCNA_results/toobig.txt")
  print(paste0("cutoff:", toobigcutoff))
  print(toobig)
sink()

newcountz = table(factor(foo$PATHWAY))
lattice::histogram(as.numeric(newcountz))
plot(density(as.numeric(newcountz)))
go_terms = foo


cat("\nPathways in REACTOME")
print(table(grepl("REACTOME", unique(go_terms$PATHWAY))))


# fit = lm(as.matrix(MEs) ~ svobj$sv + phenos$FACTOR_dx)
fit = lm(as.matrix(MEs) ~ phenos$FACTOR_dx)
thing = summary(fit)
foo = laply(thing, '[[', "coefficients")
sig_mods = list()
for(i in 1:ncol(MEs)){
  if(foo[i,ncol(foo),4] < .05) {sig_mods[[i]] = WGCNA::labels2colors(i-1)}
}
sig_mods = compact(sig_mods)

geneUniverse = colnames(datExpr)

# go_terms = go_terms[go_terms$SYMBOL %in% geneUniverse, ]
# msig=go_terms

msigg = go_terms[go_terms$SYMBOL %in% geneUniverse, ]
msigg = msigg[msigg$PATHWAY %in% names(perc_overlap), ]

# remove number of characters from long pathway names
ncharcount = nchar(msigg$PATHWAY)
toomanychar = which(ncharcount > 55)

substr = substr(x = msigg$PATHWAY[toomanychar], 1, 55)
substr = paste(substr, "...",sep="")
msigg$PATHWAY[toomanychar] = substr


# number of top ranked gene sets to report per SLIM set
topValue = 5e3


hypStats_save = list()
for( i in 1:length(sig_mods)){
  
  cat("\nRunning hypergeometric test for module:",sig_mods[[i]])
  
  gene_grab = module[module$color %in% sig_mods[[i]], ]
  gene_grab$symbol = as.character(gene_grab$symbol)
  gene_grab = gene_grab[gene_grab$symbol %in% msigg$SYMBOL, ]
  
  
  symbolSelect = gene_grab$symbol
  
  msig = msigg[!is.na(msig$SYMBOL), ]
  
  
  # filter sets
  set_count = table(as.character(msig$PATHWAY))
  
  list = unique(names(set_count))
  
  msig.keep = msig[msig$PATHWAY %in% list, ]
  
  msig.split = split(msig, msig$PATHWAY)
  
  # calculate hypergeometric input values
  universe = length(geneUniverse) # total # of genes in transcriptome 
  overlaps = lapply(msig.split, function(x) length(intersect(x$SYMBOL, symbolSelect)))
  geneSetSize = lapply(msig.split, nrow)
  geneSetOut = universe - unlist(geneSetSize)
  
  geneSetGenes = character()
  for(q in 1:length(msig.split)){
    doodad = unlist(msig.split[[q]][1])
    doodad = doodad[doodad %in% gene_grab$symbol]
    doodad = paste0(doodad,collapse = " ")
    geneSetGenes[q] = doodad
  }
  
  # hypergeomtric p-values
  hypStats = data.frame(
    Module = sig_mods[[i]],
    Genes = geneSetGenes,
    PATHWAY = names(msig.split),
    Population = universe, 
    Sample_success = unlist(overlaps),
    Population_success = unlist(geneSetSize),
    Population_fails = unlist(geneSetOut),
    Sample_size = length(unique(symbolSelect))
  )
  
  hypStats = hypStats[hypStats$Sample_success >= 1, ]
  
  # enrichment p-value test
  pvals = phyper(hypStats$Sample_success - 1, 
                 hypStats$Population_success, 
                 hypStats$Population_fails, 
                 hypStats$Sample_size, 
                 lower.tail = FALSE, log.p = FALSE)
  
  hypStats$P = pvals
  
  hypStats = hypStats[order(hypStats$P, decreasing = F), ]
  hypStats$FDR = p.adjust(hypStats$P, "fdr")
  hypStats$BONF = p.adjust(hypStats$P, "bonferroni")
  
  names(hypStats)[names(hypStats) %in% "PATHWAY"] = "SET_NAME"
  hypStats_save[[i]] = hypStats # save sum stats to list obj 
  
}

# combine hypergeometric test results across module eigengenes
hypStats_save = ldply(hypStats_save)

# calculate fold enrichment score
hypStats_save$FES = (hypStats_save$Sample_success/hypStats_save$Sample_size)/(hypStats_save$Population_success/hypStats_save$Population)

fwrite(hypStats_save,file=paste("./WGCNA_results/hypStats_new.csv",sep=""))

hypStats_trunc = hypStats_save[hypStats_save$FDR < .05,]
fwrite(hypStats_trunc,file=paste("./WGCNA_results/GSER_filtered_new.csv",sep=""))

for(i in 1:ncol(MEs)){
  png(paste0("./WGCNA_results/",names(MEs)[i],"_density.png"),res=300,units="in",height=5,width =5)
  plot(density(MEs[,i]))
  dev.off()
}
# g= ggplot(aes=MEs[1])  ##TODO work on this tutorial
