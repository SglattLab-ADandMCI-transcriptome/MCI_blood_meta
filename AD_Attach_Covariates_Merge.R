setwd("~/PsychGENe/MCI_blood_meta/")

## to get covariates, name genes, and median-average duplicate genes
## GCH w/JH

# control probes will just drop later as part of the later symbol check

## for convenience
cmString = "_merged.txt"
covsfolder = "./normalized_data/"
datafolder = "./normalized_data/"
lastname = "Hearn"  #placeholder

# make the new folder
if(!dir.exists(covsfolder)){ dir.create(covsfolder)}

require(data.table)
require(plyr)
require(org.Hs.eg.db)
require(data.table)

#list of platforms
platform = fread("./references/platform.txt")


## remove duplicates by averaging, do this after getting names
dupav <- function(x) {
  todrop = vector()
  toadd = as.data.frame(matrix(nrow=0,ncol=ncol(x)))
  names(toadd)=names(x)
  foo = 0
  names = x$PROBEID
  ln = length(unique(names))
  cat("Removing duplicates.  Total gene names:", length(unique(names)),
      " Total lines:",length(names),"\n")
  for(name in unique(names)) {
    l = which(names==name)
    if (length(l) > 1) {
      y = apply(as.matrix(x[l,-1,drop=F]),2,median)
      y = data.frame(matrix(c(0,y),nrow=1))
      y[1] = name
      names(y)=names(x)
      toadd = rbind(toadd,y,deparse.level = 0)
      todrop = c(todrop,l)
      foo = foo+1
      cat(foo,"of a possible",ln,"\r")
    }
  }
  print(paste("number of names with duplicates:",foo))
  if(length(todrop) > 0){
    x = x[-todrop,,drop=FALSE] #drop=FALSE just in case we're left with only one nonduplicate column
  }
  x = rbind(x,toadd)
  gc()
  return(x)
}

dupadd <- function(x) {
  # x = fread("~/psychgene/dupaddtest.csv")
  todrop = vector()
  toadd = as.data.frame(matrix(nrow=0,ncol=ncol(x)))
  names(toadd)=names(x)
  foo = 0
  names = x$PROBEID
  ln = length(unique(names))
  cat("Removing duplicates.  Total gene names:", length(unique(names)),
      " Total lines:",length(names),"\n")
  for(name in unique(names)) {
    l = which(names==name)
    if (length(l) > 1) {
      y = x[l,-1]
      y = sinh(y)
      y = apply(y,2,sum)
      y = asinh(y)
      y = data.frame(matrix(c(0,y),nrow=1))
      y[1] = name
      names(y)=names(x)
      toadd = rbind(toadd,y,deparse.level = 0)
      todrop = c(todrop,l)
      foo = foo+1
      cat(foo,"of a possible",ln,"\r")
    }
  }
  print(paste("number of names with duplicates:",foo))
  if(length(todrop) > 0){
    x = x[-todrop,,drop=FALSE] #drop=FALSE just in case we're left with only one nonduplicate column
  }
  x = rbind(x,toadd)
  gc()
  return(x)
}

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

covfiles = list.files(datafolder, "covariates.txt", full.names = T)
datafiles = list.files(datafolder, "normalizedProbes.txt", full.names = T)

cat("Now averaging duplicates and attaching covariates.\n")

alldata = list()
for(i in 1:length(covfiles)){
  message(covfiles[i])
  studyID = sub(".*data/","",covfiles[i])
  studyID = sub("_.*","",studyID)
  covs = fread(covfiles[i],data.table=F)
  data = fread(datafiles[i],data.table=F)
  datanames = names(data)[-1]

  ##update gene symbols
  cat("Updating gene symbols.\n")
  # read in the file and match genes
  cat(studyID,"\n")
  genes_in_data = data$PROBEID
  genes_in_data = gsub("[.]", "-", genes_in_data)
  
  nomatch = genes_in_data[!genes_in_data %in% genes$SYMBOL]
  nomatch_conv = select(org.Hs.eg.db, keys=as.character(nomatch), keytype='ALIAS', columns='ENTREZID')
  nomatch_conv$updated_hgnc = select(org.Hs.eg.db, keys=as.character(nomatch_conv$ENTREZID), keytype='ENTREZID', columns='SYMBOL')$SYMBOL
  nomatch_conv$mismatch = nomatch_conv$ALIAS == nomatch_conv$updated_hgnc
  noentrez = nomatch_conv[is.na(nomatch_conv$ENTREZID), ]

  newnames = genes_in_data
  for(q in 1:length(newnames)){
    thing = nomatch_conv$updated_hgnc[which(nomatch_conv$ALIAS == newnames[q])]
    if(length(thing) == 1 && !is.na(thing[1])){
      newnames[q] = thing
    }
  }
  cat("Updated",sum(!genes_in_data %in% newnames),"symbols.\n")
  sink(paste0("./QCplots/",studyID,"_symbol_update.txt"))
  print(data.frame(genes_in_data,newnames))
  sink()
  
  data$PROBEID = newnames
  
  if(platform$platform[which(platform$study == studyID)] == "array"){
    data = dupav(data)
  } else {
    data = dupadd(data)
  }
  
  
  probenames = data$PROBEID
  
  data = data[,-1]
  data = t(data)
  data = data.frame(data,row.names = datanames)
  names(data) = probenames
  
  newcovs = list()
  for(n in 1:length(datanames)){
    newcovs[[n]] = covs[which(covs$Sample_ID == datanames[n])[1],]
  }
  covs = ldply(newcovs)
  covs$FACTOR_age = gsub("+","",covs$FACTOR_age)
  
  thing = data.frame(covs,FACTOR_studyID = studyID, data)
  
  alldata[[i]] = thing
}
writethis = ldply(alldata)

print(table(writethis$FACTOR_studyID,writethis$FACTOR_tissue))
print(table(writethis$FACTOR_studyID,writethis$FACTOR_dx))

fwrite(writethis,paste0(covsfolder,"ADMCI",cmString))

studies = unique(writethis$FACTOR_studyID)
foo = colSums(is.na(writethis))
cat("Overlap has",sum(foo==0),"genes\n")
for(study in studies){
  bar = writethis[-which(writethis$FACTOR_studyID==study),]
  foo = colSums(is.na(bar))
  bar = writethis[which(writethis$FACTOR_studyID==study),]
  baz = colSums(is.na(bar))
  cat(study,"has",sum(baz==0),"genes, and without it the overlap has",sum(foo == 0),"\n")
}
sink("./QCplots/geneoverlapstable.txt")
foo = colSums(is.na(writethis))
cat("Overlap has",sum(foo==0),"genes\n")
for(study in studies){
  bar = writethis[-which(writethis$FACTOR_studyID==study),]
  foo = colSums(is.na(bar))
  bar = writethis[which(writethis$FACTOR_studyID==study),]
  baz = colSums(is.na(bar))
  cat(study,"has",sum(baz==0),"genes, and without it the overlap has",sum(foo == 0),"\n")
}
sink()
