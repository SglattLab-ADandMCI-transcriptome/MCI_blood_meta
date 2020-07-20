## Import data from GSE29378 aka Miller13
## GCH

# some subjects get removed for having no pheno data

require(data.table)
require(GEOquery)
require(limma)
require(forcats)

if(!exists("rawlocation")) stop("No rawlocation defined!  Run from the master import script!")

studyname = "Miller13"

cat("Reading text file from GSE29378 aka Miller13\n")
data = fread(paste0(rawlocation,"GSE29378/GSE29378_non-normalized.txt.gz"),data.table=F)

nsubj = (ncol(data) - 1) /2
Exprs = data[,(1:nsubj)*2]
Exprs = normalizeQuantiles(Exprs)
# Exprs = log2(Exprs)
Exprs = asinh(Exprs)

Exprs = data.frame(PROBEID = data$ID_REF, Exprs)

rawcovs = getGEO(filename = paste0(rawlocation,"GSE29378/GSE29378_series_matrix.txt.gz"))

covs = data.frame(Sample_ID = rawcovs@phenoData@data$geo_accession,
                  FACTOR_tissue = rawcovs@phenoData@data$source_name_ch1,
                  FACTOR_dx = rawcovs@phenoData@data$title,
                  FACTOR_age = rawcovs@phenoData@data$`age(years):ch1`,
                  FACTOR_sex = rawcovs@phenoData@data$`gender:ch1`,
                  FACTOR_ethnicity = "unknown")

covs$FACTOR_dx = unlist(lapply(strsplit(as.character(covs$FACTOR_dx)," "), `[[`,1))
covs$FACTOR_dx[covs$FACTOR_dx == "Control"] = "CTL"
print(covs$FACTOR_sex)
covs$FACTOR_sex = factor(covs$FACTOR_sex)
levels(covs$FACTOR_sex) = c("female","male")
print(covs$FACTOR_sex)

CA3 = grep("CA3",names(Exprs))
Exprs = Exprs[,-CA3]
CA3 = grep("CA3",covs$FACTOR_tissue)
covs = covs[-CA3,]
covs$FACTOR_tissue = "HIP"
print(head(covs))

# duplicates
datanames = names(Exprs)
datanames = sub("^.*?\\.","",datanames)
datanames = sub("\\..*","",datanames)
dup = which(!isUnique(datanames))
rep2 = grepl("rep2",names(Exprs)[dup])
out = dup[!rep2]
Exprs = Exprs[,-out]

## remove subjects without pheno data
datanames = names(Exprs)
datanames = sub("^.*?\\.","",datanames)
datanames = sub("\\..*","",datanames)
bads = logical()
bads[1] = F
for(i in 2:length(datanames)){
  gname = which(rawcovs@phenoData@data$`subjectnumber:ch1` == datanames[i])
  gnome = which(rawcovs@phenoData@data$geo_accession[gname] %in% covs$Sample_ID)
  arg = rawcovs@phenoData@data$geo_accession[gname[gnome]]
  if(length(arg) == 1){
    names(Exprs)[i] = arg
    bads[i] = F
  }else{
    names(Exprs)[i] = paste0("bad",i)
    bads[i] = T
  }
}
bad = which(bads)
Exprs = Exprs[,-bad]

Exprs$PROBEID = rawcovs@featureData@data$ILMN_Gene

rm(rawcovs)
rm(data)
