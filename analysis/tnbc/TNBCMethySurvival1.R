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
prj.name <- "TNBC"
path <- "."
cancername <- "BRCA"
setwd(paste(path,cancername,sep = "/"))
if(!dir.exists("raw_data"))dir.create("raw_data")
if(!dir.exists("Rdata"))dir.create("Rdata")
if(!dir.exists("figure"))dir.create("figure")

library(data.table)

# only for TNBC
cli <- fread("data/private/caseID_clinical.txt")
cli <- cli[cli$tnbc=="Yes" & (cli$race=="WHITE" | cli$race=="BLACK OR AFRICAN AMERICAN"),]
cli <- cli[,c("bcr_patient_barcode","race")]
# end

methy.value = data.table::fread(paste("./raw_data/TCGA.",cancername,".sampleMap_HumanMethylation450.gz",sep = ""), data.table = F)
methy.value[1:4,1:4]
library(tibble)
methy.value = column_to_rownames(methy.value,"sample")
colnames(methy.value)= str_sub(colnames(methy.value),1,15)
methy.value = methy.value[,as.numeric(str_sub(colnames(methy.value),14,15)) < 10]
colnames(methy.value)= str_sub(colnames(methy.value),1,12)
methy.value = methy.value[,colnames(methy.value) %in% cli$bcr_patient_barcode]
cli = cli[cli$bcr_patient_barcode %in% colnames(methy.value),]
cli$group_list = ifelse(cli$race=="BLACK OR AFRICAN AMERICAN","AA","EA")
cli = column_to_rownames(cli,"bcr_patient_barcode")

table(cli$group_list)

beta=as.matrix(methy.value)
# BiocManager::install("impute")
library(impute)
beta=impute.knn(beta) 
sum(is.na(beta))

beta=beta$data
beta=beta+0.00001

# BiocManager::install("ChAMP")
library(ChAMP)
myLoad=champ.filter(beta = beta ,pd = cli) 
dim(myLoad$beta)
save(myLoad,cancername,path,prj.name,file = paste('./Rdata/',prj.name,'_step1_myLoad.Rdata',sep = ""))
