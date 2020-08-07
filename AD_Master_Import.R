## AD brain import master script
## Subscripts to import particular ones
## Unification of naming and QA is the goal here
## GCH

setwd("~/PsychGENe/brain/")

# install if needed
# BiocManager::install(c("affy","oligo","plyr","limma","pd.hg.u133.plus.2", "R.utils",
#                      "pd.hg.u133a","pd.hg.u95av2","illuminaio","illuminaHumanv3.db",
#                      "illuminaHumanv4.db","hgu95av2.db","hgu133plus2.db", "GEOquery",
#                      "hgu133a.db","hgu133b.db","gcrma","BiocParallel","readxl","biomaRt",
#                      "edgeR"))

require(affy)
require(oligo)
require(plyr)
require(dplyr)
require(ggplot2)
require(data.table)
require(gcrma)
library(limma)
require(illuminaio)
require(illuminaHumanv3.db)
require(illuminaHumanv4.db)
require(GEOquery)
require(biomaRt)
require(edgeR)

## Will have raw data on a flash drive eventually
rawlocation = "D:/brain/data/"

# make QC plot folder
QCfolder = paste(getwd(),"/QCplots",sep="")
if(!dir.exists(QCfolder)){ dir.create(QCfolder)}

# make normalized expression folder
NormOut = paste(getwd(),"/normalized_data",sep="")
if(!dir.exists(NormOut)){ dir.create(NormOut)}

scripts = list.files("./", "AD(br|bl)_Import.*(r|R)")
cat("Scripts:\n")
print(scripts)

for(script in scripts){
  message(script)
  studyname = "STUDY_NAME"
  Exprs = data.frame()
  covs = data.frame()
  source(script)  ## output in "Exprs" study name in "studyname" covariates in "covs"
  ## quantile normalize then arcsinh
  ## PROBID as column 1, other column names as subjects
  ## names(covs) = c("Sample_ID","FACTOR_dx","FACTOR_sex","FACTOR_age","FACTOR_race","FACTOR_tissue")
  
  # Boxplot of expression
  if(ncol(Exprs) > 1){
    cat("\nGenerating plot for normalized data")
    png(paste(QCfolder,"/BOXPLOTexprnormed_",studyname,".png", sep =""), res=300,units="in",height = 6, width = 10)
    boxplot(Exprs[,-1], 
            col = 'tan', 
            las = 2, 
            xaxt = "n",
            outline = F, 
            ylab = expression("arcsinh expression (quantile normalized)"),
            main = paste0("Expression data from:",studyname))
    axis(side = 1, 
         labels = paste("S",1:(ncol(Exprs)-1),sep=""), 
         at = 1:(ncol(Exprs)-1), las = 2,cex.axis = 0.6)
    dev.off()
  }
  
  if(ncol(Exprs) > 0){
    # Write normalized expression (probe level) to data file
    cat("\nWriting normalized expression data\n")
    fwrite(Exprs,
           file = paste(NormOut,"/",studyname,"_normalizedProbes.txt", sep =""),
           quote = F, row.names = F, sep= "\t")}
  
  if(nrow(covs) > 0){
    # Write covariates to data file
    cat("\nWriting covariates\n")
    fwrite(covs,
           file = paste(NormOut,"/",studyname,"_covariates.txt", sep =""),
           quote = F, row.names = F, sep= "\t")}
  
  gc()
}
