# https://www.jianshu.com/p/982db59ec2fe
# https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247490092&idx=1&sn=a9f9e5d2a997f9bcdacee6a46006cfec&chksm=fa46df11cd315607eb0a0152d5db7e475bfe1a126f571b0fc22b1eed51d3b477ad1dcf1c38c5#rd
# https://xenabrowser.net/datapages/

# example: COAD
# clinlic data
# https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222

# Methy data
# https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FHumanMethylation450&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222

# Curated survival data
# https://xenabrowser.net/datapages/?dataset=survival%2FCOAD_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

options(stringsAsFactors = F)

rm(list = ls())

cancer.top <- "BRCA" # Breast Invasive Carcinoma
cancername <- "TNBC"
cancer.top <- "COAD" # Colon Adenocarcinoma
cancer.top <- "LUAD" # Lung Adenocarcinoma
cancer.top <- "PRAD" # Prostate Adenocarcinoma

if(!exists('cancername')){
  cancername <- cancer.top
}

prj.name <- "Tumor.vs.Ctrl"
path <- "./NotOnlyPairSample"
workingpath <- paste(path,cancername,sep='/')

if(!dir.exists(workingpath))dir.create(workingpath)
setwd(workingpath)

if(!dir.exists("Rdata"))dir.create("Rdata")
if(!dir.exists("figure"))dir.create("figure")

library(data.table)
beta.raw = data.table::fread(paste("../../",cancer.top,"/raw_data/TCGA.",cancer.top,".sampleMap_HumanMethylation450.gz",sep = ""), data.table = F)

if(cancername=='TNBC'){
  cli <- fread("data/private/caseID_clinical.txt")
  # cli <- cli[cli$tnbc=="Yes" & (cli$race=="WHITE" | cli$race=="BLACK OR AFRICAN AMERICAN"),]
  cli <- cli[cli$tnbc=="Yes",]
  cli <- cli[,c("bcr_patient_barcode","race")]
  beta.raw = beta.raw[,str_sub(colnames(beta.raw),1,12) %in% c(cli$bcr_patient_barcode,'sample')]
}

library(tibble)
beta.raw = column_to_rownames(beta.raw,"sample")
library(stringr)
colnames(beta.raw)= str_sub(colnames(beta.raw),1,15)

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = FALSE)
cat('raw: probes samples','\n')
print(dim(beta.raw))
sink()

pd <- data.frame(sampleID=colnames(beta.raw))
pd$bcr_patient_barcode = str_sub(pd$sampleID,1,12)
pd$Sample_Group = ifelse(as.numeric(str_sub(pd$sampleID,14,15))<10,"Tumor","Ctrl")
rownames(pd) <- pd$sampleID

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n')
print(table(pd$Sample_Group))
sink()

beta=as.matrix(beta.raw)
# BiocManager::install("impute")
library(impute)
beta=impute.knn(beta) 
sum(is.na(beta))

beta=beta$data
beta=beta+0.00001

# BiocManager::install("ChAMP")
library(ChAMP)
myLoad=champ.filter(beta = beta ,pd = pd) 

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','filted: probes samples','\n')
print(dim(myLoad$beta))
sink()

save(myLoad,workingpath,cancer.top,cancername,prj.name,file = paste('./Rdata/',prj.name,'_step1_Load.Rdata',sep = ""))

# gz1 <- gzfile(paste('./Rdata/',prj.name,'_step1_beta.gz',sep = ""), "w")
# write.csv(myLoad[[1]], gz1)
# close(gz1)
