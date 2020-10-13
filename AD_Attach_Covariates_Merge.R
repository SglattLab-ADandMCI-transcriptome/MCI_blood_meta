setwd("~/PsychGENe/MCI_blood_meta/")

## to get covariates, name genes, and median-average duplicate genes
## GCH w/JH

# control probes will just drop later as part of the later symbol check

## for convenience
cmString = "_merged.txt"
covsfolder = "./normalized_data/"
datafolder = "./normalized_data/"
lastname = "Hearn"  #placeholder

#list of platforms

# make the new folder
if(!dir.exists(covsfolder)){ dir.create(covsfolder)}

require(data.table)
require(plyr)

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
  x = x[-todrop,,drop=FALSE] #drop=FALSE just in case we're left with only one nonduplicate column
  x = rbind(x,toadd)
  gc()
  return(x)
}

cat("Now averaging duplicates and attaching covariates.\n")

covfiles = list.files(datafolder, "covariates.txt", full.names = T)
datafiles = list.files(datafolder, "normalizedProbes.txt", full.names = T)

alldata = list()
for(i in 1:length(covfiles)){
  message(covfiles[i])
  studyID = sub(".*data/","",covfiles[i])
  studyID = sub("_.*","",studyID)
  covs = fread(covfiles[i],data.table=F)
  data = fread(datafiles[i],data.table=F)
  datanames = names(data)[-1]
  
  data = dupav(data)
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

fwrite(writethis,paste0(covsfolder,"ADMCI",cmString))
