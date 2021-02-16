setwd("~/PsychGENe/MCI_blood_meta/")

## TODO gsg??
## TODO soft threshold
#### maybe just 12 according to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
##### with a picksoft just to verify reasonableness
## TODO deconvolution???


## WGCNA network preservation analysis for MCI blood studies
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

Wfolder = "./WGCNA_npa_results"
Pfolder = "./WGCNA_npa_results"

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


## Remove ROSMAP batch 3 due to low gene overlap with other studies
datExpr = datExpr[-which(phenos$FACTOR_studyID=="ROSMAP3"),]
phenos = phenos[-which(phenos$FACTOR_studyID=="ROSMAP3"),]

# badgenes = which(colSums(is.na(datExpr)) > 1)
# if(length(badgenes) > 0){
#   datExpr = datExpr[,-badgenes]
# }


gooddx = grep("(MCI|CTL)$",phenos$FACTOR_dx)
datExpr = datExpr[gooddx,]
samples = samples[gooddx]
phenos = phenos[gooddx,]

datCTL = datExpr[which(phenos$FACTOR_dx == "CTL"),]
datMCI = datExpr[which(phenos$FACTOR_dx == "MCI"),]

phenosCTL = phenos[which(phenos$FACTOR_dx == "CTL"),]
phenosMCI = phenos[which(phenos$FACTOR_dx == "MCI"),]


## run goodsamplesgenes to ensure suitability of samples and genes
cat("\nRunning GoodSamplesGenes for CTL.")
gsg = goodSamplesGenes(datCTL, verbose = 3);
bad = c("Bad Genes: ",names(datCTL)[!gsg$goodGenes],"\n",
        "Bad Samples: ",rownames(datCTL)[!gsg$goodSamples])
# print(bad)
cat("\n",length(names(datCTL)[!gsg$goodGenes]),"bad genes and",
    length(rownames(datCTL)[!gsg$goodSamples]),"bad samples.\n")
cat("\nWriting bad samples and genes to file.\n")
write(bad, file = paste(Wfolder,"/WGCNA_SamplesGenes_Removed_CTL.txt",sep=""))

datCTL = datCTL[gsg$goodSamples,gsg$goodGenes]
phenosCTL = phenosCTL[gsg$goodSamples,]
rownames(datCTL) = phenosCTL$FACTOR_sampleID

cat("\nRunning GoodSamplesGenes for MCI.")
gsg = goodSamplesGenes(datMCI, verbose = 3);
bad = c("Bad Genes: ",names(datMCI)[!gsg$goodGenes],"\n",
        "Bad Samples: ",rownames(datMCI)[!gsg$goodSamples])
# print(bad)
cat("\n",length(names(datMCI)[!gsg$goodGenes]),"bad genes and",
    length(rownames(datMCI)[!gsg$goodSamples]),"bad samples.\n")
cat("\nWriting bad samples and genes to file.\n")
write(bad, file = paste(Wfolder,"/WGCNA_SamplesGenes_Removed_MCI.txt",sep=""))

datMCI = datMCI[gsg$goodSamples,gsg$goodGenes]
phenosMCI = phenosMCI[gsg$goodSamples,]
rownames(datMCI) = phenosMCI$FACTOR_sampleID



multiExpr = list(CTL = list(data = datCTL), MCI = list(data = datMCI))

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



## get ctl modules
## one-step automated gene network analysis
# 1. find optimal soft-threshold power for network construction
cat("\nFinding soft threshold...\n")
# colnames(datExpr) = convertKeys(colnames(datExpr))

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datCTL,
                        powerVector = powers,
                        RsquaredCut = .8,
                        corFnc="bicor",
                        networkType="signed",
                        verbose = 1)

sft0 = sft

sft0

png(paste(Pfolder,"/softthreshold_CTL.png",sep=""),res=300,units="in",height=6,width=6)

  par(mfrow=c(1,1))
  par(mar = c(5.1, 5.1,5.1,2.1))
  plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
       las = 1, cex.lab = 1.1, cex.axis = 1.1,
       xlab="Soft Threshold (power)",ylab=expression(paste("SFT, signed R"^2)),
       type="n",main=paste("Scale independence"))
  text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
       labels=powers,col="red")
  text(15,.2,paste("Power used:",sft0$powerEstimate))
  abline(h=0.80,col="dodgerblue3", lty=2)    #CHOOSE A  R^2 CUT-OFF OF H
  # plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
  #      xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
  # text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

cat("\nConstructing adjacency and assigning modules.\n")
adjacencyPre = adjacency(datCTL,
                         # power=14,  #A Power used to build a consensus network
                         power = sft0$powerEstimate,
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
sft.power = sft0$powerEstimate
sft.power = 12; ##
deepSplit = 2;
minModuleSize = 30;

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,deepSplit]
table(modulesPRE)

#Check to see if network modules can be cut and merged...
cat("Calculating potential eigengenes.\n")
MEList = moduleEigengenes(datCTL, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs, use='p')
METree = hclust(as.dist(MEDiss), method ="average")

# 4.one-step automated network construction
cat("Constructing network.\n")
net = blockwiseModules(datCTL,
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
                       saveTOMFileBase = "CTL",
                       verbose = 3,
                       maxBlockSize = 5000)

# 5. save the network to a .Rdata file for future use
saveRDS(net, file = paste(Wfolder,"/WGCNA_CTL.Rdata",sep=""))

net = readRDS(paste(Wfolder,"/WGCNA_CTL.Rdata",sep=""))

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
  genesInModule = colnames(datCTL)[which(moduleColors %in% colors[[y]])]
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



## MCI module network
# 1. find optimal soft-threshold power for network construction
cat("\nFinding soft threshold...\n")
# colnames(datExpr) = convertKeys(colnames(datExpr))

powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datMCI,
                        powerVector = powers,
                        RsquaredCut = .8,
                        corFnc="bicor",
                        networkType="signed",
                        verbose = 1)

sft0 = sft

sft0

png(paste(Pfolder,"/softthreshold_MCI.png",sep=""),res=300,units="in",height=6,width=6)

par(mfrow=c(1,1))
par(mar = c(5.1, 5.1,5.1,2.1))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
     las = 1, cex.lab = 1.1, cex.axis = 1.1,
     xlab="Soft Threshold (power)",ylab=expression(paste("SFT, signed R"^2)),
     type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],
     labels=powers,col="red")
text(15,.2,paste("Power used:",sft0$powerEstimate))
abline(h=0.80,col="dodgerblue3", lty=2)    #CHOOSE A  R^2 CUT-OFF OF H
# plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
# text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

cat("\nConstructing adjacency and assigning modules.\n")
adjacencyPre = adjacency(datMCI,
                         # power=14,  #A Power used to build a consensus network
                         power = sft0$powerEstimate,
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

pdf(paste(Wfolder,"/WGCNA_Consensus_DeepSplit_MCI.pdf",sep=""), height=10,width=25);
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

# 3. set parameters for network algorithm
sft.power = sft0$powerEstimate
sft.power = 12; ##
deepSplit = 2;
minModuleSize = 30;

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,deepSplit]
table(modulesPRE)

#Check to see if network modules can be cut and merged...
cat("Calculating potential eigengenes.\n")
MEList = moduleEigengenes(datMCI, colors=modulesPRE)
MEsMCI=MEList$eigengenes
MEDiss = 1-cor(MEsMCI, use='p')
METree = hclust(as.dist(MEDiss), method ="average")

# 4.one-step automated network construction
cat("Constructing MCI network.\n")
netMCI = blockwiseModules(datMCI,
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
                       saveTOMFileBase = "MCI",
                       verbose = 3,
                       maxBlockSize = 5000)

# 5. save the network to a .Rdata file for future use
saveRDS(netMCI, file = paste(Wfolder,"/WGCNA_MCI.Rdata",sep=""))
netMCI = readRDS(paste(Wfolder,"/WGCNA_MCI.Rdata",sep=""))

# 6. extract network meta-data and eigengenes
moduleLabelsMCI = netMCI$colors
moduleColorsMCI = labels2colors(netMCI$colors)
foo = order(as.numeric(substr(names(netMCI$MEs),3,stop=4)))
MEsMCI = netMCI$MEs[,foo]

metadataMCI = data.frame(colors = moduleColorsMCI, labels = paste("ME",netMCI$colors, sep=""))
metadataMCI = metadataMCI[!duplicated(metadataMCI$colors), ]

mergedColorsMCI = labels2colors(netMCI$colors)

# PLOT DENDROGRAM
png(paste(Pfolder,"/wgcna_dendrogram_MCI.png",sep=""), res=300, units="in",height = 5.5, width = 8)
plotDendroAndColors(netMCI$dendrograms[[1]], las = 1, cex.axis = 1.2, cex.lab = 1.1,
                    mergedColorsMCI[netMCI$blockGenes[[1]]],
                    "Module colors", cex.colorLabels = 1,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# SAVE MODULE MEMBERSHIP
moduleMCI = list()
colors = unique(moduleColorsMCI)
for(y in 1:length(colors)){
  genesInModule = colnames(datMCI)[which(moduleColorsMCI %in% colors[[y]])]
  moduleMCI[[y]] = data.frame(color = colors[[y]], label = unique(netMCI$colors)[[y]], symbol = genesInModule)
}
moduleMCI = ldply(moduleMCI)

mcounts = data.frame(table(moduleMCI$label))
mcounts = data.frame(c(0:(nrow(mcounts)-1)),labels2colors(c(0:(nrow(mcounts)-1))),mcounts$Freq)
names(mcounts) = c("ME","color","frequency")
fwrite(mcounts,file = paste0(Wfolder,"/wgcna_module-membership_counts_MCI.txt"))

fwrite(data.table(moduleMCI),
       file = paste(Wfolder,"/wgcna_module-membership_MCI.txt", sep=""),
       quote = F, row.names = F, sep = "\t")

moduleMCI = fread(paste(Wfolder,"/wgcna_module-membership_MCI.txt", sep=""),data.table=F)












## network preservation analysis part
message("Beginning network preservation analysis.")
multiColor = list(CTL = module$color, MCI = moduleMCI$color)
mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1,
                        networkType = "signed",
                        verbose = 3)

saveRDS(mp, file = paste0(Wfolder,"/modulePreservation.RData"))

mp = readRDS(paste0(Wfolder,"/modulePreservation.RData"))

## TODO this is just copypasted from a tutorial, i need to figure out what i want
## TODO something with the labels, i should recreate this in ggplot
sink(paste0(Pfolder,"/Preservation.txt"))
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1],
                 mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
sink()

##TODO dying here right now, Inf zsummary
png(paste(Pfolder,"/Preservation.png",sep=""), res=300, units="in",
    height = 5.5, width = 11)
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
# sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2){
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.1);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();






## save preservation statistics
for(i in 1:4){
  thing = mp$preservation[i][[1]][[1]][[2]]
  thingname = names(mp$preservation[i])
  fwrite(thing,paste0(Pfolder,"/Preservation_Statistics_",thingname,".csv"),row.names=T)
}






## TODO this is mad too -- do we need to gsg the main group, then split it?
## create and save overlap table
## TODO kME style
## TODO can use adjacency matrices from above?

# ## TODO i guess this doesn't work if it's not the same labels?  where do the labels go?
# overlap = overlapTable(mergedColors,mergedColorsMCI) # rows label1, columns label2
# odf = data.frame(overlap$countTable)
# fwrite(odf,paste0(Pfolder,"/overlaps.csv"),row.names=T)
# opdf = data.frame(overlap$pTable)
# fwrite(opdf,paste0(Pfolder,"/overlaps_p.csv"),row.names=T)

connCTL = intramodularConnectivity.fromExpr(datCTL,module$color)
connMCI = intramodularConnectivity.fromExpr(datMCI,moduleMCI$label)

overlapKME = overlapTableUsingKME(datCTL, datMCI,
                                  module$color, moduleMCI$label)



# stop("stuff after here has more than one group")

### saving a copy of the environment
# save.image("~/psychgene/wgcnaworking.rdata.RData")




png(paste(Pfolder,"/ME_networks_CTL.png",sep=""), res=300, units="in",height = 5.5, width = 8)
plotEigengeneNetworks(MEs,setLabels=rownames(datCTL))
dev.off()

png(paste(Pfolder,"/ME_networks_MCI.png",sep=""), res=300, units="in",height = 5.5, width = 8)
plotEigengeneNetworks(MEsMCI,setLabels=rownames(datMCI))
dev.off()


# conn = intramodularConnectivity.fromExpr(datExpr,module$color)
# 
# ##intramodularConnectivity and hub genes, with significance from the meta-analysis
# cat("Comparing connectivity.\n")
# fwrite(conn,paste(Wfolder,"/intramodularConnectivity.txt",sep=""))
# meta = fread("./meta_analysis/whole_blood_MCI_meta.txt",data.table=F)
# meta$t.value = meta$arcsinh / meta$SE
# module$sig = NA
# for(i in 1:length(module$symbol)){
#   foo = meta$t.value[which(meta$GeneSymbol == module$symbol[i])]
#   if(length(foo) == 1) module$sig[i] = foo
# }
# 
# for (i in 1:(length(colors)-1)){ ## because grey isn't connected lol
#   foo = module$label==i
#   png(paste(Pfolder,"/module_connectivity_Tscore_",i,".png",sep=""), res=300, units="in",height = 5.5, width = 8)
#     verboseScatterplot(
#       conn$kWithin[foo],
#       module$sig[foo],
#       col=module$color[foo],
#       main=paste("Module ",i,", ",labels2colors(i),":",sep=""),
#       xlab = "Internal Connectivity", ylab = "T Score", abline = TRUE)
#   dev.off()
# }
# 
# png(paste(Pfolder,"/module_Tscore.png",sep=""), res=300, units="in",height = 5.5, width = 8)
# par(las = 2)
# plotModuleSignificance(
#   module$sig,
#   module$color,
#   boxplot = TRUE,
#   main = "Average T-Score by Module,",
#   ylab = "T-Score",
#   xlas = 1
#   )
# dev.off()


cat("Generating top hub genes per module.\n")
topHubs = list()
for (i in 1:length(colors)) {
  foo = module$label==i-1
  foo = data.frame(module$symbol[foo],connCTL$kWithin[foo],stringsAsFactors = F)
  foo = foo[order(foo$connCTL.kWithin.foo.,decreasing=T),]
  names(foo) = c("symbol","kWithin")
  topHubs[[i]] = head(foo,10)
  names(topHubs)[i]=labels2colors(i-1)
}
sink(paste(Wfolder,"/CTLtopHubs.txt",sep=""))
  print(topHubs)
sink()
cat("\nPrinted to CTLtopHubs.txt.")

topHubs = list()
for (i in 1:length(colors)) {
  foo = moduleMCI$label==i-1
  foo = data.frame(moduleMCI$symbol[foo],connMCI$kWithin[foo],stringsAsFactors = F)
  foo = foo[order(foo$connMCI.kWithin.foo.,decreasing=T),]
  names(foo) = c("symbol","kWithin")
  topHubs[[i]] = head(foo,10)
  names(topHubs)[i]=labels2colors(i-1)
}
sink(paste(Wfolder,"/MCItopHubs.txt",sep=""))
print(topHubs)
sink()
cat("\nPrinted to MCItopHubs.txt.")


### saving a copy of the environment
save.image("~/psychgene/wgcnaworking.rdata.RData")


stop("blarg")


## Gene set enrichment analysis for interesting module eigengenes
cat("\nPerforming gene set enrichment analysis for CTL modules.")

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
go_terms$in_data = ifelse(go_terms$SYMBOL %in% colnames(datCTL), "in", "out")

count_int = table(go_terms$PATHWAY, go_terms$in_data)
perc_overlap = count_int[,1]/rowSums(count_int)
perc_overlap = perc_overlap[perc_overlap > 0.5]

go_terms = go_terms[go_terms$SYMBOL %in% colnames(datCTL), ] # only look at genes in the datCTL object
go_terms = go_terms[go_terms$PATHWAY %in% names(perc_overlap), ]

cat("\nPathways in REACTOME")
print(table(grepl("REACTOME", unique(go_terms$PATHWAY))))


## TODO instead of 'significant' we need to go with 'disrupted'
## pValuesHypergeo is a feature in overlapsusingKME, I think by module
## though that may be more for figuring out which modules correspond to which
## we can also sort by/note the zsummary and median rank

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
fit = lm(as.matrix(MEs) ~ svobj$sv + phenos$FACTOR_dx)   ## TODO svobj is what
### ahhh TODO this whole section is focused on dx comparison, gotta change it
## maybe just output by kwithin, but include a bit on the mci/severity?
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
  
  # hypergeomtric p-values
  hypStats = data.frame(
    Module = sig_mods[[i]],
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
