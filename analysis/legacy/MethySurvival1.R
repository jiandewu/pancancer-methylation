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
path <- "."
cancername <- "BRCA"
prj.name="BRCA"
setwd(paste(path,cancername,sep = "/"))
if(!dir.exists("raw_data"))dir.create("raw_data")
if(!dir.exists("Rdata"))dir.create("Rdata")
if(!dir.exists("figure"))dir.create("figure")

library(data.table)
pd <- fread(paste("./raw_data/TCGA.",cancername,".sampleMap_",cancername,"_clinicalMatrix",sep = ""))
#colnames(pd)
table(pd$anatomic_neoplasm_subdivision)

pd <- pd[,c("sampleID","anatomic_neoplasm_subdivision")]
library(stringr)
pd$sampleID = str_sub(pd$sampleID,1,15)
pd$patient = str_sub(pd$sampleID,1,12)
pd$group_list = ifelse(as.numeric(str_sub(pd$sampleID,14,15))<10,"Tumor","Normal")
table(pd$group_list)

tp = pd$patient[pd$group_list =="Normal"]
np = pd$patient[pd$group_list =="Tumor"]
okp = intersect(tp,np)

pd <- pd[pd$patient %in% okp,]
rownames(pd) <- pd$sampleID; dim(pd)

library(textshape)
# only for TNBC
cli <- fread("data/private/caseID_clinical.txt")
cli <- cli[cli$tnbc=="Yes",]
cli <- cli[,c("bcr_patient_barcode","race")]

black = cli$bcr_patient_barcode[cli$race=="BLACK OR AFRICAN AMERICAN"]
white = cli$bcr_patient_barcode[cli$race=="WHITE"]
allrace = c(black, white)
cli = column_to_rownames(cli,"bcr_patient_barcode")
# end

a = data.table::fread(paste("./raw_data/TCGA.",cancername,".sampleMap_HumanMethylation450.gz",sep = ""), data.table = F)
a[1:4,1:4]
library(tibble)
a = column_to_rownames(a,"sample")
colnames(a)= str_sub(colnames(a),1,15)

k1 = str_sub(colnames(a),1,12) %in% pd$patient
table(k1)
a = a[,k1]

# 
tmp1 = str_sub(colnames(a),1,12) %in% allrace
a = a[,tmp1]
#

patient = str_sub(colnames(a),1,12)
group_list = ifelse(as.numeric(str_sub(colnames(a),14,15))<10,"Tumor","Normal")
table(group_list)

tp = patient[group_list =="Normal"]
np = patient[group_list =="Tumor"]
okp = intersect(tp,np);length(okp)

a = a[,patient %in% okp]; dim(a)

pd = pd[match(str_sub(colnames(a),1,15),pd$sampleID),]
identical(str_sub(colnames(a),1,15),pd$sampleID)


beta=as.matrix(a)
# BiocManager::install("impute")
library(impute)
beta=impute.knn(beta) 
sum(is.na(beta))

beta=beta$data
beta=beta+0.00001

# BiocManager::install("ChAMP")
library(ChAMP)
myLoad=champ.filter(beta = beta ,pd = pd) 
dim(myLoad$beta)
save(path,cancername,myLoad,file = './Rdata/step1_myLoad.Rdata')
