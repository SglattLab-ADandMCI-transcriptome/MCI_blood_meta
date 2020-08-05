adni_import_start_time = Sys.time()
## imports ADNI data.  it has already been RMA normalized.

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "ADNI"

cat("Reading data files from ADNI\n")

## extract and isolate expression data
adni.expr = read.csv(paste0(rawlocation,"blood/adni/data/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),header=FALSE,stringsAsFactors = F)

agestuff = rbind(adni.expr[8,],adni.expr[3,]) #need this later for age calculation

cat("~~Extracting data and naming probes~~\n")

colnames(adni.expr) = adni.expr[3,]
adni.expr = adni.expr[-c(1:9),-c(1:2,748)] # the final column is a descriptor of the genes
adni.norm = data.frame(t(adni.expr),stringsAsFactors = F)
foo = adni.expr[,1]
bar = sub("\\s\\|\\|.*","",foo)
colnames(adni.norm) = bar ##changes multiply-gened probe names to just one gene
adni.norm = adni.norm[-1,]

rm(foo)
rm(bar)

adni.median=adni.norm


## extract covariates
cat("~~Extracting covariates~~\n")
adni.covs = read.csv(paste0(rawlocation,"blood/adni/data/ADNIMERGE.csv"),stringsAsFactors = F)
IDs = rownames(adni.norm)

genders = matrix(nrow = 0, ncol = 2)
for(foo in IDs){
  quux = adni.covs$PTGENDER[adni.covs$PTID %in% foo]
  gender = c(foo,unique(quux))
  genders = rbind(genders,gender)
}
genders = tolower(genders)
rownames(genders) = genders[,1]
genders = as.data.frame(genders[-745,2])
colnames(genders) = "PTGENDER"
table(genders)
# genders
# female   male 
# 336    408  
336/(336+408)
# [1] 0.4516129

ages = matrix(nrow = 0, ncol = 2)
for(foo in IDs){
  age = adni.covs$AGE[adni.covs$PTID %in% foo]
  exam = adni.covs$EXAMDATE[adni.covs$PTID %in% foo]
  bloodyear = agestuff[1,agestuff[2,]==foo]
  extra = as.numeric(bloodyear) - as.numeric(substr(exam[1],1,4))
  realage = c(foo,age[1]+extra)
  ages = rbind(ages,realage)
}
rownames(ages) = ages[,1]
ages = ages[,2,drop=F]
colnames(ages) = "BLOODAGE"
ages=as.numeric(ages)
summary(ages)
# BLOODAGE    
# Min.   :55.00  
# 1st Qu.:69.20  
# Median :75.10  
# Mean   :74.59  
# 3rd Qu.:80.00  
# Max.   :93.60  

Dx = matrix(nrow = 0, ncol = 2)
for(foo in IDs){
  diag = adni.covs$DX_bl[adni.covs$PTID %in% foo]
  dg = c(foo,unique(diag))
  Dx = rbind(Dx,dg)
}
rownames(Dx) = Dx[,1]
Dx = as.data.frame(Dx[-745,2])
colnames(Dx) = "DX_bl"
table(Dx)
# AD   CN EMCI LMCI 
# 43  260  215  226 
Dx[Dx=="CN"] = "CTL"
Dx[Dx=="EMCI"] = "MCI"
Dx[Dx=="LMCI"] = "MCI"
table(Dx)

eths = matrix(nrow = 0, ncol = 3)
for(foo in IDs){
  eth = adni.covs$PTRACCAT[adni.covs$PTID %in% foo]
  hl = adni.covs$PTETHCAT[adni.covs$PTID %in% foo]
  eh = c(foo,unique(eth),unique(hl))
  eths = rbind(eths,eh)
}
rownames(eths) = eths[,1]
eths = as.data.frame(eths[-745,-1])
table(eths)
# V1                    V2     
# Am Indian/Alaskan:  2   Hisp/Latino    : 17  
# Asian            : 11   Not Hisp/Latino:724  
# Black            : 30   Unknown        :  3  
# Hawaiian/Other PI:  2                        
# More than one    :  7                        
# Unknown          :  2                        
# White            :690                        
foo = paste(eths[,1],eths[,2])
foo = matrix(foo,ncol=1)
rownames(foo) = rownames(eths)
eths = foo
table(eths)
# V1     
# White Not Hisp/Latino            :675  
# Black Not Hisp/Latino            : 29  
# White Hisp/Latino                : 13  
# Asian Not Hisp/Latino            : 11  
# More than one Not Hisp/Latino    :  5  
# Am Indian/Alaskan Not Hisp/Latino:  2  
# (Other)                          :  9

foo = gsub(" Not Hisp/Latino","",eths)
foo = gsub(" Unknown","",foo)
foo[grep("Latino",foo)] = "latino"
foo[grep("Hawaiian",foo)] = "islander"
foo[grep("than one",foo)] = "unknown"
foo[grep("Alaskan",foo)] = "nativeamerican"
eths = tolower(foo)

covs = data.frame(IDs, Dx,genders,ages,eths,tissue = rep("whole blood",times=nrow(genders)))
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")


foo = as.data.frame(t(sapply(adni.norm, as.numeric)))
foo = 2^foo      ## RMA log2s it
foo = asinh(foo) ## convert to arcsinh space
Exprs = data.frame(PROBEID = colnames(adni.norm),foo)
colnames(Exprs) = c("PROBEID",rownames(adni.norm))


rm(foo)
rm(adni.expr)
rm(adni.covs)
rm(adni.median)
rm(adni.norm)

rm(age)
rm(bloodyear)
rm(dg)
rm(diag)
rm(eh)
rm(eth)
rm(exam)
rm(extra)
rm(gender)
rm(hl)
rm(IDs)
rm(quux)
rm(realage)

rm(ages)
rm(Dx)
rm(agestuff)
rm(eths)
rm(genders)


adni_import_processing_time = Sys.time() - adni_import_start_time
rm(adni_import_start_time)
cat("ADNI import processing time: ",format(adni_import_processing_time), "\n")
