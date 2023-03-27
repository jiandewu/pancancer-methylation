# https://www.jianshu.com/p/1aa756ef73c7

rm(list=setdiff(ls(), c('prj.name','workingpath','cancer.top','cancername','myLoad')))
if (!exists('myLoad')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step1_Load.Rdata',sep = ''))
}

setwd(workingpath)

norm_file = paste("./Rdata/",prj.name,"_step2_Norm.Rdata",sep = "")
if(!file.exists(norm_file)){
  library(ChAMP)
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=8)
  save(myNorm,file = norm_file)
} else {
  load(norm_file)
}

num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
table(num.na)
names(num.na) = colnames(myNorm)
delsamples = names(num.na[num.na>0])
keep = setdiff(colnames(myNorm),c(delsamples))
myNorm = myNorm[,keep]

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','normalized: probes samples','\n')
print(dim(myNorm))
sink()

pd = myLoad$pd
pd = pd[pd$sampleID %in% keep,]

## clinical table
library(TCGAbiolinks)

clin.file <- paste("../../",cancer.top,"/raw_data/TCGA-",cancer.top,"_clinical.csv",sep = "")
if(!file.exists(clin.file)){
  clin <- GDCquery_clinic(paste("TCGA-",cancer.top,sep = ''), type = "clinical", save.csv = TRUE)
  file.copy(from = paste("TCGA-",cancer.top,"_clinical.csv",sep = ""),
            to = paste("../../",cancer.top,"/raw_data/TCGA-",cancer.top,"_clinical.csv",sep = ""))
  file.remove(from = paste0("TCGA-",cancer.top,"_clinical.csv"))
} else {
  clin = data.table::fread(clin.file, data.table = F)
}

clin.succ <- clin[c('bcr_patient_barcode','race','gender','age_at_index','vital_status','ajcc_pathologic_stage')]
clin.succ = clin.succ[clin.succ$bcr_patient_barcode %in% pd$bcr_patient_barcode,]

library(dplyr)
library(stringr)
clin.succ <- clin.succ %>% mutate(race = str_replace(race, 'american indian or alaska native|asian|not reported', 'other'))
clin.succ <- clin.succ %>% mutate(ajcc_pathologic_stage = str_replace(ajcc_pathologic_stage, 'A|B|C', ''))
clin.succ$age_range = ifelse(clin.succ$age_at_index<65,"<65","≥65")
## ------------------
pd = pd[pd$bcr_patient_barcode %in% clin.succ$bcr_patient_barcode,]

library(plyr)
pd <- join(pd,clin.succ,by = 'bcr_patient_barcode')
rownames(pd) <- pd$sampleID

myNorm = myNorm[,str_sub(colnames(myNorm),1,12) %in% pd$bcr_patient_barcode]

identical(row.names(pd),colnames(myNorm))

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','intergrated clinical: probes samples','\n')
print(dim(myNorm))

print(table(pd$race,pd$Sample_Group))
print(table(pd$age_range,pd$Sample_Group))
print(table(pd$gender,pd$Sample_Group))
print(table(pd$ajcc_pathologic_stage,pd$Sample_Group))

print(table(clin.succ$race))
print(table(clin.succ$age_range))
print(table(clin.succ$gender))
print(table(clin.succ$ajcc_pathologic_stage))
sink()

save(cancer.top,cancername,prj.name,workingpath,pd,myNorm,file = paste('./Rdata/',prj.name,'_step2_filterNAClinNorm_pd.Rdata',sep = ""))

beta.norm <- t(myNorm)

# BiocManager::install("FactoMineR")
library(FactoMineR)
beta.norm.pca <- PCA(beta.norm, graph = FALSE) 
# BiocManager::install("factoextra")
library(factoextra) 
png(file=paste('./figure/',prj.name,'_pca.png',sep = ""))
# dev.list()
dev.set(which = 2)
par(cex=1, cex.main=2, cex.lab = 1.5, cex.sub=0.8)
fviz_pca_ind(beta.norm.pca,
             geom.ind = "point", 
             col.ind = pd$Sample_Group, 
             addEllipses = TRUE, 
             legend.title = "Groups",
             title = paste(cancername,prj.name,'PCA',sep = ' - '))
dev.copy(which = 4)
# dev.cur()
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
ac=data.frame(group=pd$Sample_Group)
rownames(ac)=colnames(myNorm)  

library(pheatmap)
png(file=paste('./figure/',prj.name,'_orig_heatmap.png',sep = ""))
dev.set(which = 2)
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac)
dev.copy(which = 4)
dev.off()

png(file=paste('./figure/',prj.name,'_orig_heatmap_2.png',sep = ""))
dev.set(which = 2)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F)
dev.copy(which = 4)
dev.off()


