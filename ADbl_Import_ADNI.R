adni_import_start_time = Sys.time()
## imports ADNI data.  it has already been RMA normalized.

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "ADNI"

cat("Reading data files from ADNI\n")

## extract and isolate expression data
adni.expr = read.csv(paste0(rawlocation,"blood/adni/data/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),header=FALSE,stringsAsFactors = F)

agestuff = rbind(adni.expr[8,],adni.expr[3,]) #need this later for age calculation
viscode = adni.expr[c(2,3),]
phases = adni.expr[1,]

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

## this is incorrect
# Dx = matrix(nrow = 0, ncol = 2)
# for(foo in IDs){
#   diag = adni.covs$DX_bl[adni.covs$PTID %in% foo]
#   dg = c(foo,unique(diag))
#   Dx = rbind(Dx,dg)
# }
# rownames(Dx) = Dx[,1]
# Dx = as.data.frame(Dx[-745,2])
# colnames(Dx) = "DX_bl"
# table(Dx)
# # AD   CN EMCI LMCI 
# # 43  260  215  226 
# Dx[Dx=="CN"] = "CTL"
# Dx[Dx=="EMCI"] = "MCI"
# Dx[Dx=="LMCI"] = "MCI"
# table(Dx)

## Get diagnosis
## special thanks to JH and JH  ... for helping me get it right
#samples are from ADNI2 and ADNIGO phases
#1: fix visit code for ADNI2
viscode = t(viscode)
viscode = data.frame(viscode)
names(viscode) = viscode[1,]
viscode = viscode[-c(1:3,nrow(viscode)),]
table(viscode$Visit)
# bl m03 m12 m36 m48 m60 v02 v03 v04 v05 v06 v11 
# 1 111   1   1   6 123  53   1 357   8   1  72  10 
viscode$Visit[viscode$Visit == "v03"] = "bl"

visit = fread(paste0(rawlocation,"blood/adni/Assessments/ADNI2_VISITID.csv"), stringsAsFactors = F)
Dx = data.frame(SubjectID = rownames(adni.norm))
for(i in 1:nrow(Dx)){
  sIDi = grep(Dx$SubjectID[i],adni.covs$PTID)
  sID = adni.covs$RID[sIDi]
  sID = unique(sID)
  if(length(sID)>1) stop("bad ptID/RID match")
  Dx$RID[i] = sID
}
table(Dx$RID %in% visit$RID)
# FALSE  TRUE 
# 28   716 
all(viscode$SubjectID == Dx$SubjectID) # TRUE
Dx = data.frame(Dx,Visit = viscode$Visit)

## which are adni2 and not specifically for blood, translate viscode
notgo = grepl("v",Dx$Visit)
adni2 = Dx[notgo,]
adni2$mer = paste(adni2$RID, adni2$Visit, sep = "_")
nrow(adni2) #92
visit$mer = paste(visit$RID,visit$VISCODE, sep = "_")
mervisit = merge(adni2,visit, by = "mer")
nrow(mervisit) #92
for(i in 1:nrow(Dx)){
  if(Dx$SubjectID[i] %in% mervisit$SubjectID){
    Dx$Visit[i] = mervisit$VISCODE2[which(mervisit$SubjectID == Dx$SubjectID[i])]
  }
}
table(Dx$Visit)
# bl   m03   m06   m12   m36   m48   m60   m72   m84 scmri 
# 468     9     1     8     6   134   102    13     2     1 

#2: create coID as phase_visit_PTID
phases = data.frame(Phase = t(phases))
phases = phases[-c(1:3,nrow(phases)),]
Dx$Phase = phases
for(i in 1:nrow(Dx)){
  Dx$coID[i] = paste0(Dx$Phase[i],"_",Dx$Visit[i],"_",Dx$SubjectID[i])
}
for(i in 1:nrow(adni.covs)){
  adni.covs$coID[i] = paste0(adni.covs$COLPROT[i],"_",adni.covs$VISCODE[i],"_",adni.covs$PTID[i])
}
table(Dx$coID %in% adni.covs$coID)
# FALSE  TRUE 
# 1   743 
Dx$coID[which(Dx$Visit == "scmri")] = "ADNI2_bl_009_S_4564" ##see below
for(i in 1:nrow(Dx)){
  DX_bl = NA
  DX_bl = adni.covs$DX_bl[which(adni.covs$coID == Dx$coID[i])]
  DX = adni.covs$DX[which(adni.covs$coID == Dx$coID[i])]
  # if(length(DX_bl)==0) DX_bl = "Unknown"
  # if(length(DX)==0) DX = "Unknown"
  Dx$Dx_bl[i] = DX_bl
  Dx$Dx_new[i] = DX
  Dx$AGE_bl[i] = adni.covs$AGE[which(adni.covs$coID == Dx$coID[i])]
  Dx$Years_bl[i] = adni.covs$Years_bl[which(adni.covs$coID == Dx$coID[i])]
  #TODO updates etc
}
table(Dx$Dx_bl,Dx$Dx_new)
#              CN Dementia MCI Unknown
# AD       43   0        0   0       0
# CN      175  70        1  14       0
# EMCI    103   1        0 111       0
# LMCI    128   4       49  44       0
# Unknown   0   0        0   0       1
Dx$Dx_bl[which(Dx$SubjectID == "009_S_4564")] = "MCI"
for(i in 1:nrow(Dx)){
  if(Dx$Dx_bl[i] == "AD") Dx$Dx[i] = "AD"
  if(Dx$Dx_bl[i] == "CN") Dx$Dx[i] = "CTL"
  if(Dx$Dx_bl[i] == "EMCI") Dx$Dx[i] = "MCI"
  if(Dx$Dx_bl[i] == "LMCI") Dx$Dx[i] = "MCI"
  if(Dx$Dx_new[i] == "Dementia") Dx$Dx[i] = "AD"
  if(Dx$Dx_new[i] == "MCI") Dx$Dx[i] = "MCI"
  if(Dx$Dx_new[i] == "CN") Dx$Dx[i] = "CTL"
}
table(Dx$Dx)
# AD  CN MCI 
# 116 247 381 

## Get age at blood draw
Dx$AGE = Dx$AGE_bl + Dx$Years_bl

## TODO age and dx down below

#### 009_S_4564 has a visit code of "scmri", but only one entry in adnimerge:
# > adni.covs[4630,]
# RID       PTID VISCODE SITE COLPROT ORIGPROT   EXAMDATE DX_bl  AGE PTGENDER PTEDUCAT        PTETHCAT PTRACCAT PTMARRY
# 4630 4564 009_S_4564      bl    9   ADNI2    ADNI2 2012-03-20  LMCI 64.9   Female       12 Not Hisp/Latino    White Married
# APOE4 FDG PIB AV45 CDRSB ADAS11 ADAS13 MMSE RAVLT_immediate RAVLT_learning RAVLT_forgetting RAVLT_perc_forgetting FAQ MOCA
# 4630     1  NA  NA   NA     3      7     12   27              41              6               10                   100   8   24
# EcogPtMem EcogPtLang EcogPtVisspat EcogPtPlan EcogPtOrgan EcogPtDivatt EcogPtTotal EcogSPMem EcogSPLang EcogSPVisspat
# 4630     1.375          1             1          1     1.16667            1     1.10256      3.25       1.75          1.25
# EcogSPPlan EcogSPOrgan EcogSPDivatt EcogSPTotal FLDSTRENG FSVERSION Ventricles Hippocampus WholeBrain Entorhinal Fusiform
# 4630        1.6     3.16667            2     2.28571                             NA          NA         NA         NA       NA
# MidTemp ICV DX EXAMDATE_bl CDRSB_bl ADAS11_bl ADAS13_bl MMSE_bl RAVLT_immediate_bl RAVLT_learning_bl RAVLT_forgetting_bl
# 4630      NA  NA     2012-03-20        3         7        12      27                 41                 6                  10
# RAVLT_perc_forgetting_bl FAQ_bl FLDSTRENG_bl FSVERSION_bl Ventricles_bl Hippocampus_bl WholeBrain_bl Entorhinal_bl
# 4630                      100      8                                      NA             NA            NA            NA
# Fusiform_bl MidTemp_bl ICV_bl MOCA_bl EcogPtMem_bl EcogPtLang_bl EcogPtVisspat_bl EcogPtPlan_bl EcogPtOrgan_bl
# 4630          NA         NA     NA      24        1.375             1                1             1        1.16667
# EcogPtDivatt_bl EcogPtTotal_bl EcogSPMem_bl EcogSPLang_bl EcogSPVisspat_bl EcogSPPlan_bl EcogSPOrgan_bl EcogSPDivatt_bl
# 4630               1        1.10256         3.25          1.75             1.25           1.6        3.16667               2
# EcogSPTotal_bl FDG_bl PIB_bl AV45_bl Years_bl Month_bl Month M          update_stamp                coID
# 4630        2.28571     NA     NA      NA        0        0     0 0 2017-08-04 22:50:48.0 ADNI2_bl_009_S_4564




#### TODO above like age

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
foo = eths[,1]

foo[grep("Hawaiian",foo)] = "pacificislander"
foo[grep("than one",foo)] = "other"
foo[grep("Alaskan",foo)] = "nativeamerican"
eths = tolower(foo)

covs = data.frame(IDs, Dx$Dx,genders,Dx$AGE,eths,tissue = rep("whole_blood",times=nrow(genders)))
names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")
covs[1:10,]


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
# rm(dg)
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
