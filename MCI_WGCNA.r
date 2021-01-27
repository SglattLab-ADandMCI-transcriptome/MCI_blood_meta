setwd("~/PsychGENe/MCI_blood_meta/")

## WGCNA for MCI blood studies
## GCH

require(data.table)
require(WGCNA)
require(plyr)
require(ggplot2)
require(sva)
require(limma)
require(sas7bdat)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(GO.db)
require(reactome.db)

Wfolder = "./WGCNA_results"
Pfolder = "./WGCNA_results"

if(!dir.exists(Wfolder)){
  dir.create(Wfolder)
}

if(!dir.exists(Pfolder)){
  dir.create(Pfolder)
}

##  WGCNA COMMANDS
disableWGCNAThreads() # number of processors reported by system will be used

## Load raw data
cat("\nLoading raw data.")
datRaw = fread("./data_for_analysis/whole_blood_ScaledWithFactors_OutliersRemoved_allstudies.txt",data.table=F)

datExpr = datRaw[,-which(grepl("FACTOR_",names(datRaw)))]
genes = names(datExpr)
samples = datRaw$FACTOR_sampleID
phenos = datRaw[which(grepl("FACTOR_",names(datRaw)))]


# ####################
# #####test area
# datExpr = datExpr[-which(phenos$FACTOR_studyID=="Shigemizu20"),]
# phenos = phenos[-which(phenos$FACTOR_studyID=="Shigemizu20"),]
# ####################


badgenes = which(colSums(is.na(datExpr)) > 1)
if(length(badgenes) > 0){
  datExpr = datExpr[,-badgenes]
}
gooddx = grep("(MCI|CTL)$",phenos$FACTOR_dx)
datExpr = datExpr[gooddx,]
samples = samples[gooddx]
phenos = phenos[gooddx,]

## run goodsamplesgenes to ensure suitability of samples and genes
cat("\nRunning GoodSamplesGenes.")
gsg = goodSamplesGenes(datExpr, verbose = 3);
bad = c("Bad Genes: ",names(datExpr)[!gsg$goodGenes],"\n",
        "Bad Samples: ",rownames(datExpr)[!gsg$goodSamples])
# print(bad)
cat("\n",length(names(datExpr)[!gsg$goodGenes]),"bad genes and",
    length(rownames(datExpr)[!gsg$goodSamples]),"bad samples.\n")
cat("\nWriting bad samples and genes to file.\n")
write(bad, file = paste(Wfolder,"/WGCNA_SamplesGenes_Removed.txt",sep=""))

datExpr = datExpr[gsg$goodSamples,gsg$goodGenes]
genes = genes[gsg$goodGenes]
samples = samples[gsg$goodSamples]
phenos = phenos[gsg$goodSamples,]
rownames(datExpr) = samples

# ## ComBat with studies as batches, so must delete all genes with ANY NA
# foo = colSums(is.na(datExpr))==0
# write(names(datExpr)[!foo], file = paste(Wfolder,"/Genes_anymissing_Removed.txt",sep=""))
# datExpr = datExpr[,foo]
# met = t(datExpr)
# mod = model.matrix(~ FACTOR_dx,phenos)
# foo = ComBat(met,batch=phenos$FACTOR_studyID,mod=mod)
# genes = names(datExpr)
# samples = rownames(datExpr)
# datExpr = data.frame(t(foo))
# names(datExpr) = genes
# rownames(datExpr) = samples
#
# cat("\nWGCNA genes:",length(genes),"\n")

## one-step automated gene network analysis
# 1. find optimal soft-threshold power for network construction
cat("\nFinding soft threshold...\n")
# colnames(datExpr) = convertKeys(colnames(datExpr))

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr,
                        powerVector = powers,
                        corFnc="bicor",
                        networkType="signed",
                        verbose = 1)

sft0 = sft

sft0

png(paste(Pfolder,"/softthreshold.png",sep=""),res=300,units="in",height=6,width=6)

  par(mfrow=c(1,1))
  par(mar = c(5.1, 5.1,5.1,2.1))
  plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
       las = 1, cex.lab = 1.1, cex.axis = 1.1,
       xlab="Soft Threshold (power)",ylab=expression(paste("SFT, signed R"^2)),
       type="n",main=paste("Scale independence"))
  text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
       labels=powers,col="red")
  abline(h=0.80,col="dodgerblue3", lty=2)    #CHOOSE A  R^2 CUT-OFF OF H
  # plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
  #      xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
  # text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

cat("\nConstructing adjacency and assigning modules.\n")
adjacencyPre = adjacency((datExpr),
                         power=14,  #A Power used to build a consensus network
                         type="signed")
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

# MODULE ASSIGNMENTS
mColorh=NULL
for (ds in 0:4){
  tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
                      minClusterSize = (25), cutHeight = 0.9999,
                      deepSplit = ds, distM = dissTOMPre)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}

pdf(paste(Wfolder,"/WGCNA_Consensus_DeepSplit.pdf",sep=""), height=10,width=25);
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

# 3. set parameters for network algorithm
sft.power = 12; ##
deepSplit = 2;
minModuleSize = 30;

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,deepSplit]
table(modulesPRE)

#Check to see if network modules can be cut and merged...
cat("Calculating potential eigengenes.\n")
MEList = moduleEigengenes(datExpr, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs, use='p')
METree = hclust(as.dist(MEDiss), method ="average")

# 4.one-step automated network construction
cat("Constructing network.\n")
net = blockwiseModules(datExpr,
                       power = sft.power,
                       networkType = "signed",
                       deepSplit= deepSplit,
                       TOMType = "signed",
                       minModuleSize = minModuleSize,
                       minCoreKME = 0.5,
                       minCoreKMESize = minModuleSize/3,
                       minKMEtoStay=0,
                       reassignThreshold = 1e-6,
                       mergeCutHeight = 0.25,
                       detectCutHeight = 0.995,
                       numericLabels = TRUE,
                       corType = 'bicor',
                       pamRespectsDendro = FALSE,
                       pamStage = TRUE,
                       saveTOMs = TRUE,
                       verbose = 3,
                       maxBlockSize = 5000)

# 5. save the network to a .Rdata file for future use
saveRDS(net, file = paste(Wfolder,"/WGCNA.Rdata",sep=""))

net = readRDS(paste(Wfolder,"/WGCNA.Rdata",sep=""))

# 6. extract network meta-data and eigengenes
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
foo = order(as.numeric(substr(names(net$MEs),3,stop=4)))
MEs = net$MEs[,foo]

metadata = data.frame(colors = moduleColors, labels = paste("ME",net$colors, sep=""))
metadata = metadata[!duplicated(metadata$colors), ]

mergedColors = labels2colors(net$colors)

# PLOT DENDROGRAM
png(paste(Pfolder,"/wgcna_dendrogram.png",sep=""), res=300, units="in",height = 5.5, width = 8)
plotDendroAndColors(net$dendrograms[[1]], las = 1, cex.axis = 1.2, cex.lab = 1.1,
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors", cex.colorLabels = 1,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# SAVE MODULE MEMBERSHIP
module = list()
colors = unique(moduleColors)
for(y in 1:length(colors)){
  genesInModule = colnames(datExpr)[which(moduleColors %in% colors[[y]])]
  module[[y]] = data.frame(color = colors[[y]], label = unique(net$colors)[[y]], symbol = genesInModule)
}
module = ldply(module)

mcounts = data.frame(table(module$label))
mcounts = data.frame(c(0:(nrow(mcounts)-1)),labels2colors(c(0:(nrow(mcounts)-1))),mcounts$Freq)
names(mcounts) = c("ME","color","frequency")
fwrite(mcounts,file = paste0(Wfolder,"/wgcna_module-membership_counts.txt"))

fwrite(data.table(module),
       file = paste(Wfolder,"/wgcna_module-membership.txt", sep=""),
       quote = F, row.names = F, sep = "\t")

module = fread(paste(Wfolder,"/wgcna_module-membership.txt", sep=""),data.table=F)

## Generate plots for differential expresison among groups
cat("Plotting associations.\n")
phenos = phenos[phenos$FACTOR_sampleID %in% rownames(datExpr),]

## The boxplots are the interquartile range, whiskers are extended to 1.5 IQR
## The lines are medians, and the notches (from ggplot2 docs):
## In a notched box plot, the notches extend 1.58 * IQR / sqrt(n).
## This gives a roughly 95% confidence interval for comparing medians.
## See McGill et al. (1978) for more details.

toplot = list()
col = numeric()
for(i in 1:length(colors)){
  toplot[[i*2-1]] = MEs[phenos$FACTOR_dx == "CTL",i]
  toplot[[i*2]] = MEs[phenos$FACTOR_dx == "AD",i]
  names(toplot)[[i*2-1]] = paste(i-1,"CTL")
  names(toplot)[[i*2]] = paste(labels2colors(i-1),"AD")
  col[c(i*2-1,i*2)] = labels2colors(i-1)
  # png(paste(Pfolder,"/ME",i,".png",sep=""), res=300, units="in",height = 5.5, width = 8)
  # # plot(1:nrow(MEs),MEs[,i],##residuals(fit)[,i],
  # #      col=as.factor(phenos$FACTOR_dx),
  # #      ylab="Eigengene Value",
  # #      xlab="Sample")
  # # # abline(fit$coefficients[1,i],fit$coefficients[2,i])
  # # legend("topleft",legend = levels(phenos$FACTOR_dx),fill=c(1,2))
  # boxplot(MEs[,i] ~ phenos$FACTOR_dx, outline=F)
  # title(main = paste("Module eigengene",i,"(",colors[i],")"))
  # dev.off()
}
col[col=="black"] = "dimgray"

png(paste(Pfolder,"/ME_diffex.png",sep=""), res=300, units="in",height = 5.5, width = 8)
# par(mar = c(7, 4, 4, 2) + 0.1)
# boxplot(toplot,names=names(toplot),pars=list(las = 2),col=col,ylab = "ME Value")
# title(main = "Module Eigengene Differential Expression")

## ggplot boxplot module eigengene expression across affection groups
df_plot = data.frame(Dx = phenos$FACTOR_dx, MEs)
df_plot = data.frame(melt(df_plot))
colnames(df_plot)[2] = 'label'
module_df = module
module_df$label = paste("ME", module_df$label, sep = "")
module_df = module_df[!duplicated(module_df$label), ]
module_df = module_df[,!colnames(module_df) %in% c('symbol','sig')]
df_plot = merge(df_plot, module_df, by='label')
g = ggplot(df_plot, aes(x = label, y = value, fill = label, group = interaction(color, Dx))) +
  geom_boxplot(notch = TRUE,
               outlier.shape = NA) +
  theme_bw() +
  theme(legend.position="none",axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12,angle=90,hjust=1,vjust=.3),
        axis.title=element_text(size =12),
        panel.border = element_rect(size = 2, fill = NA),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  ylab("Module eigengene expression") +
  xlab("Module eigengenes") +
  scale_x_discrete(breaks=df_plot$label, labels=df_plot$color) +
  scale_fill_manual(values = labels2colors(0:(length(unique(df_plot$color))-1))) +
  labs(title = "Eigengene Differential Expression",
       subtitle = paste0(unique(df_plot$Dx)[1]," (left) vs ",unique(df_plot$Dx)[2]," (right)")) +
  coord_cartesian(ylim = c(-.06,.06))
print(g)
dev.off()

png(paste(Pfolder,"/ME_networks.png",sep=""), res=300, units="in",height = 5.5, width = 8)
plotEigengeneNetworks(MEs,setLabels=rownames(datExpr))
dev.off()


conn = intramodularConnectivity.fromExpr(datExpr,module$color)

##intramodularConnectivity and hub genes, with significance from the meta-analysis
cat("Comparing connectivity.\n")
fwrite(conn,paste(Wfolder,"/intramodularConnectivity.txt",sep=""))
meta = fread("./meta_analysis/whole_blood_MCI_meta.txt",data.table=F)
meta$t.value = meta$arcsinh / meta$SE
module$sig = NA
for(i in 1:length(module$symbol)){
  foo = meta$t.value[which(meta$GeneSymbol == module$symbol[i])]
  if(length(foo) == 1) module$sig[i] = foo
}

for (i in 1:(length(colors)-1)){ ## because grey isn't connected lol
  foo = module$label==i
  png(paste(Pfolder,"/module_connectivity_Tscore_",i,".png",sep=""), res=300, units="in",height = 5.5, width = 8)
    verboseScatterplot(
      conn$kWithin[foo],
      module$sig[foo],
      col=module$color[foo],
      main=paste("Module ",i,", ",labels2colors(i),":",sep=""),
      xlab = "Internal Connectivity", ylab = "T Score", abline = TRUE)
  dev.off()
}

png(paste(Pfolder,"/module_Tscore.png",sep=""), res=300, units="in",height = 5.5, width = 8)
par(las = 2)
plotModuleSignificance(
  module$sig,
  module$color,
  boxplot = TRUE,
  main = "Average T-Score by Module,",
  ylab = "T-Score",
  xlas = 1
  )
dev.off()


cat("Generating top hub genes per module.\n")
topHubs = list()
for (i in 1:length(colors)) {
  foo = module$label==i-1
  foo = data.frame(module$symbol[foo],conn$kWithin[foo],stringsAsFactors = F)
  foo = foo[order(foo$conn.kWithin.foo.,decreasing=T),]
  names(foo) = c("symbol","kWithin")
  topHubs[[i]] = head(foo,10)
  names(topHubs)[i]=labels2colors(i-1)
}
sink(paste(Wfolder,"/topHubs.txt",sep=""))
  print(topHubs)
sink()
cat("\nPrinted to topHubs.txt.")


## Create linear model with eigengenes
cat("\nCreating linear model.")
mod = model.matrix(~ FACTOR_dx,phenos)
fit = lm(as.matrix(MEs) ~ phenos$FACTOR_dx)

coefs = coef(summary(fit))
res = residuals(fit)
coefs = ldply(coefs)
coefs = coefs[seq(2,nrow(coefs),2),]
colnames(coefs)[grepl("Pr..", colnames(coefs))] = "P"
coefs$FDR = p.adjust(coefs$P, 'fdr')
coefs$BONF = p.adjust(coefs$P, 'bonferroni')
coefs$label = 0:(nrow(coefs)-1)
coefs$color = labels2colors(coefs$label)
coefs = coefs[order(coefs$P, decreasing = F), ]

png(paste(Wfolder,"/WGCNA_volcano.png",sep=""),res=300,units="in",height = 5, width = 6)
print(
  ggplot(coefs, aes(x = coefs$Estimate, y = -log10(coefs$P))) +
    geom_point(size = 1, col = coefs$color) +
    theme_classic() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), panel.border = element_rect(fill=NA,size = 1.5)) +
    geom_hline(yintercept = -log10(.05), col = 'orange', linetype = "dashed") +
    xlab( expression(paste("Model Eigengene Difference in Expression"))) +
    ylab( expression(paste("-log"[10],"(P-value)"))) +
    geom_hline(aes(yintercept = -log10(.05/nrow(coefs))), col = 'red', lwd = 0.25,linetype=1) +
    geom_vline(aes(xintercept = 0), col = 'black',lwd=0.3, linetype=2) +
    ggrepel::geom_text_repel(aes(label = coefs$color), size = 3, segment.size = 0.2, segment.colour = "grey", fontface = 3)
)
dev.off()

fwrite(coefs,file=paste(Wfolder,"/ME_LM.txt",sep=""))


## sva a couple SVs and account for them in the lm fit with eigengenes
# mod = model.matrix(~ FACTOR_dx,phenos)
mod = model.matrix(~ FACTOR_dx,phenos)
svobj = sva(t(MEs), mod)
fit = lm(as.matrix(MEs) ~ svobj$sv)

# boxplot(residuals(fit))
png(paste(Pfolder,"/ME_SVA_residuals.png",sep=""), res=300, units="in",height = 5.5, width = 8)
## ggplot boxplot module eigengene residuals across affection groups
df_plot = data.frame(Dx = phenos$FACTOR_dx, residuals(fit))
df_plot = data.frame(melt(df_plot))
colnames(df_plot)[2] = 'label'
module_df = module
module_df$label = paste("ME", module_df$label, sep = "")
module_df = module_df[!duplicated(module_df$label), ]
module_df = module_df[,!colnames(module_df) %in% c('symbol','sig')]
df_plot = merge(df_plot, module_df, by='label')

g = ggplot(df_plot, aes(x = label, y = value, fill = label, group = interaction(color, Dx))) +
  geom_boxplot(notch = TRUE,
               outlier.shape = NA) +
  theme_bw() +
  theme(legend.position="none",axis.text.y=element_text(size = 12),
        axis.text.x=element_text(size = 12,angle=90,hjust=1,vjust=.3),
        axis.title=element_text(size =12),
        panel.border = element_rect(size = 2, fill = NA),
        plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  ylab("Module eigengene residuals") +
  xlab("Module eigengenes") +
  scale_x_discrete(breaks=df_plot$label, labels=df_plot$color) +
  scale_fill_manual(values = labels2colors(0:(length(unique(df_plot$color))-1))) +
  labs(title = paste("Eigengene Residuals after",svobj$n.sv,"SV(s) removed"),
       subtitle = paste0(unique(df_plot$Dx)[1]," (left) vs ",unique(df_plot$Dx)[2]," (right)")) +
  coord_cartesian(ylim = c(-.06,.06))
print(g)
dev.off()

save(svobj,file=paste(Wfolder,"/ME_SVA_OBJ.Rdata",sep=""))
fwrite(coef(summary(fit)),file=paste(Wfolder,"/ME_SVA_LM.txt",sep=""))


## Gene set enrichment analysis for interesting module eigengenes
cat("\nPerforming gene set enrichment analysis.")

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

go_terms = go_terms[go_terms$SYMBOL %in% colnames(datExpr), ] # only look at genes in the datExpr object
go_terms = go_terms[go_terms$PATHWAY %in% names(perc_overlap), ]

cat("\nPathways in REACTOME")
print(table(grepl("REACTOME", unique(go_terms$PATHWAY))))


## If there are no sig mods, make sure that there's only case and control in dx
## decide some significant modules
# fit = lm(as.matrix(MEs) ~ svobj$sv)
# thing = residuals(fit)
# sig_mods = list()
# for(i in 1:ncol(MEs)){
#   # foo = t.test(datExpr[phenos$FACTOR_dx == "AD",module$color %in% color],datExpr[phenos$FACTOR_dx == "CTL",module$color %in% color])
#   foo = t.test(thing[phenos$FACTOR_dx == "AD",i],thing[phenos$FACTOR_dx == "CTL",i])
#   if(foo$p.value < .05) {sig_mods[[i]] = labels2colors(i-1)}
# }
# sig_mods = compact(sig_mods)
#---------
# sig_mods = unique(as.character(graphDF$Module[graphDF$FDR < .05]))
# sig_mods = as.character(metadata$colors[metadata$labels %in% sig_mods])
#-----------
fit = lm(as.matrix(MEs) ~ svobj$sv + phenos$FACTOR_dx)
thing = summary(fit)
foo = laply(thing, '[[', "coefficients")
sig_mods = list()
for(i in 1:ncol(MEs)){
  if(foo[i,ncol(foo),4] < .05) {sig_mods[[i]] = labels2colors(i-1)}
}
sig_mods = compact(sig_mods)

geneUniverse = colnames(datExpr)

go_terms = go_terms[go_terms$SYMBOL %in% geneUniverse, ]
msig=go_terms

# remove number of characters from long pathway names
ncharcount = nchar(msig$PATHWAY)
toomanychar = which(ncharcount > 55)

substr = substr(x = msig$PATHWAY[toomanychar], 1, 55)
substr = paste(substr, "...",sep="")
msig$PATHWAY[toomanychar] = substr


# number of top ranked gene sets to report per SLIM set
topValue = 5e3


hypStats_save = list()
for( i in 1:length(sig_mods)){
  
  cat("\nRunning hypergeometric test for module:",sig_mods[[i]])
  
  gene_grab = module[module$color %in% sig_mods[[i]], ]
  gene_grab$symbol = as.character(gene_grab$symbol)
  gene_grab = gene_grab[gene_grab$symbol %in% msig$SYMBOL, ]
  
  
  symbolSelect = gene_grab$symbol
  
  msig=go_terms
  
  msig = msig[!is.na(msig$SYMBOL), ]
  
  
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

fwrite(hypStats_save,file=paste(Wfolder,"/hypStats.csv",sep=""))

hypStats_trunc = hypStats_save[hypStats_save$FDR < .05,]
fwrite(hypStats_trunc,file=paste(Wfolder,"/GSER_filtered.csv",sep=""))
