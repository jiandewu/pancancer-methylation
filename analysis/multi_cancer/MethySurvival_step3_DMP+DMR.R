
rm(list=setdiff(ls(), c('cancer.top','prj.name','workingpath','cancername','pd','myNorm')))
if (!exists('pd') | !exists('myNorm')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step2_filterNAClinNorm_pd.Rdata',sep = ''))
}

setwd(workingpath)

library(ChAMP)

DMP <- champ.DMP(beta=myNorm,
                 pheno=pd$Sample_Group,
                 adjPVal=1,
                 adjust.method="BH",
                 compare.group=NULL,
                 arraytype="450K")
DMP <- DMP[[1]]
DMP = data.frame("CpG"=rownames(DMP),DMP)
DMP$DMPrefCtrl <- ifelse(DMP$adj.P.Val<0.05, 
                         ifelse(DMP$logFC<0,'Hypo','Hyper'),'NOT') 

table(DMP$DMPrefCtrl)

# save(DMP,df_DMP,file = paste('./Rdata/',prj.name,'_step3_DMP.Rdata',sep = ""))

# DMRs
DMRBumphunter <- champ.DMR(beta=myNorm,pheno=pd$Sample_Group,method="Bumphunter")
DMRBumphunter <- DMRBumphunter[[1]]
DMRBumphunter = data.frame("region"=rownames(DMRBumphunter),DMRBumphunter)
DMRBumphunter$DMRrefCtrl <- ifelse(DMRBumphunter$fwer<0.05, 
                      ifelse(DMRBumphunter$value<0,'Hypo','Hyper'),'NOT') 
table(DMRBumphunter$DMRrefCtrl)
# value:  Mean value of all smooth coef in one candidate bump. is the average of the estimated coefficient 
#         (usually representing the difference between two groups) in the region
# area - is the absolute value of the sum of the estimated coefficients in the region
# For 2 groups, value column is the average of the coefficients of the linear model 
# calculated for each site within the DMR while area column is the sum of all coefficients. 
# This can be easily checked by calculating the difference, between the 2 groups, 
# of the mean value (M or beta) for each methylation site of the DMR and then averaging 
# (value column) or summing (area column) these differences.
# Thus, a positive value indicates more methylation in the non-control group

# cluster - the genome is divided into smaller regions, called clusters. We search for regions in each cluster to avoid going over a long stretch of the genome with no data. The value provided here is the ID for the cluster.
# indexStart - is the index of the probe in the original methylation matrix where the region starts
# indexEnd- is the index of the probe in the original methylation matrix where the region ends
# L: Numbers of CpGs contained in candidate bump.
# L - is the number of probes in the region
# clusterL - number of probes in the cluster
# p.value:  Proportion of random bumps show most CpGs and higher mean value then this DMR. These do need multiple comparison corrections.
# fwer:  Proportion that a random run would generate one such DMR shows most CpGs and higher mean value. These do not need multiple comparison corrections
# p.valueArea:  Proportion of random bumps show higher abs sum value then this DMR.
# fwerArea: Proportion that a random run would generate one such DMR show higher abs sum value.

# A cluster is defined purely based on the array design (the spacing between probes).  
# It is a collection of probes where the between probe distance (on be a cluster.
# A region is (always) a subset of a cluster and is a set of probes which are
# differentially methylated. So that depends heavily on the actual observed data.


library(dplyr)
library(stringr)
DMRBumphunter <- DMRBumphunter %>% mutate(seqnames = str_replace(seqnames, 'chr', ''))
DMRBumphunter$seqnames = as.integer(DMRBumphunter$seqnames)

library(sqldf)
DMP = sqldf("select DMP.*, DMRBumphunter.region, DMRBumphunter.DMRrefCtrl, DMRBumphunter.L as noCpGs, DMRBumphunter.fwer as DMPfwer from DMP
                left join DMRBumphunter
                on DMP.CHR = DMRBumphunter.seqnames and DMP.MAPINFO between DMRBumphunter.start and DMRBumphunter.end")
DMP <- DMP %>% mutate(DMRDMPrefCtrl = case_when(is.na(DMRrefCtrl) & DMPrefCtrl == 'NOT' ~ 'NOT BOTH',
                                                is.na(DMRrefCtrl) ~ 'DMP ONLY',
                                                DMPrefCtrl == 'NOT' ~ 'DMR ONLY',
                                                DMPrefCtrl != DMRrefCtrl ~ 'INV MUTE',
                                                TRUE ~ paste(DMPrefCtrl,'BOTH')))

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','*********** DMP: adj p-value < 0.05, |deltaBeta| > 0\n')
print(table(DMP$DMPrefCtrl))
cat('\n')
print(table(cut(DMP$adj.P.Val,breaks=c(0,10^-45,10^-35,10^-25,10^-15,10^-5,0.005,0.01,0.05,1))))
cat('\n','*********** DMR: adj p-value < 0.05\n')
print(table(DMRBumphunter$DMRrefCtrl))
cat('\n')
print(table(cut(DMRBumphunter$fwer,breaks=c(0,10^-45,10^-35,10^-25,10^-15,10^-5,0.005,0.01,0.05,1))))
cat('\n','*********** DMP integrate DMR: adj p-value < 0.05, |deltaBeta| > 0; DMR.p-value < 0.05\n')
print(table(DMP$DMRDMPrefCtrl))
sink()

write.table(DMRBumphunter, 
            file = paste('./figure/',prj.name,'_DMR.Bumphunter.txt',sep = ""), 
            sep="\t", quote=F, row.names=F)
write.table(DMP, 
            file = paste('./figure/',prj.name,'_DMP.DMR.txt',sep = ""), 
            sep="\t", quote=F, row.names=F)

## DMRcate ------------------------
#myDMR <- champ.DMR(beta=myNorm,pheno=pd$Sample_Group,method="DMRcate")

# https://github.com/YuanTian1991/ChAMP/issues/18
DMRcate_manual_run = function(beta_matrix,
                              pheno,
                              arraytype,
                              dist = 2, 
                              mafcut = 0.05,
                              fdr = 0.05,
                              lambda = 1000,
                              C = 2)
{
  
  require(DMRcate)
  myMs <- logit2(beta_matrix)
  myMs <- rmSNPandCH(myMs, dist = dist, mafcut = mafcut)
  design <- model.matrix(~ pheno)
  
  if(arraytype == "450K"){
    myannotation <- cpg.annotate(datatype = "array",
                                 fdr = 0.05, 
                                 myMs, 
                                 design = design,
                                 coef = ncol(design), 
                                 analysis.type = "differential",
                                 annotation = c(array = "IlluminaHumanMethylation450k", 
                                                annotation = "ilmn12.hg19"), 
                                 what = "M")
  } else {
    myannotation <- cpg.annotate(datatype = "array",
                                 fdr = 0.05, 
                                 myMs, 
                                 design = design,
                                 coef = ncol(design), 
                                 analysis.type = "differential",
                                 annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                                annotation = "ilm10b4.hg19"), 
                                 what = "M")
  }
  
  #commenting lines with issues
  # M <- do.call("cbind", lapply(myannotation, as.data.frame)) #don't know why this code is here
  # colnames(M) <- names(myannotation) #don't know why this code is here
  
  #removed unused mc.cores argument
  dmrcoutput <- dmrcate(myannotation, lambda = lambda, C = C) #mc.cores = cores - this argument no longer applies 
  #data(dmrcatedata) - unsused line of code
  DMR <- as.data.frame(extractRanges(dmrcoutput, genome = "hg19"))
  rownames(DMR) <- paste("DMR", 1:nrow(DMR), sep="_")
  
  message("DMRcate detected ", nrow(DMR)," DMRs with mafcut as = ", mafcut, ".")
  if(nrow(DMR) == 0) 
    stop("No DMR detected.")
  
  return(DMR)
}

#to run
DMRcateDMRs <- DMRcate_manual_run(beta=myNorm,pheno=pd$Sample_Group, "450K")
DMRcateDMRs = data.frame("region"=rownames(DMRcateDMRs),DMRcateDMRs)

write.table(DMRcateDMRs, 
            file = paste('./figure/',prj.name,'_DMR.cate.txt',sep = ""), 
            sep="\t", quote=F, row.names=F)

# save(DMRcateDMRs,file = paste('./Rdata/',prj.name,'_step3_DMR.cate.Rdata',sep = ""))


sum(DMRcateDMRs$no.cpgs)
