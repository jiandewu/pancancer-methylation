# https://www.jianshu.com/p/1aa756ef73c7

rm(list = ls())
# path <- "."
# cancername <- "PRAD"
load('./BRCA/Rdata/step1_myLoad.Rdata')
setwd(paste(path,cancername,sep = "/"))
myLoad$beta[1:4,1:4]
myLoad$pd[1:10,1:2]

# norm_file = paste("./Rdata/",prj.name,"_step2_champ_myNorm.Rdata",sep = "")
norm_file = paste("./Rdata/","step2_champ_myNorm.Rdata",sep = "")
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
dt = names(num.na[num.na>0])
library(stringr)
dn = str_replace(dt,"-01","-11")
keep = setdiff(colnames(myNorm),c(dt,dn))
myNorm = myNorm[,keep]
pd = myLoad$pd
# pd = pd[row.names(pd) %in% keep,]
pd = pd[pd$sampleID %in% keep,]
identical(row.names(pd),colnames(myNorm))

dat <- t(myNorm)

group_list=pd$group_list
table(group_list)

# BiocManager::install("FactoMineR")
library(FactoMineR)
dat.pca <- PCA(dat, graph = TRUE) 
# BiocManager::install("factoextra")
library(factoextra) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE, 
             legend.title = "Groups")

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)  

library(pheatmap)
# load("./Rdata/step2_filtered_pd_myNorm.Rdata")
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac)


png(file="./figure/Heatmap.png",width=900, height=900)
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac)
dev.off()
# dev.off()
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F)

save(cg,ac,cancername,pd,myNorm,cancername,path,prj.name,file = paste('./Rdata/',prj.name,'_step2_filtered_pd_myNorm.Rdata',sep = ""))

