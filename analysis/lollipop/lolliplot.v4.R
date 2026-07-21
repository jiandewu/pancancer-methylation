# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
library(trackViewer)
library(dplyr)
library(stringr)

options(dplyr.summarise.inform = FALSE)
setwd("./BRCA/figure/White.TNBC.vs.ctrl_lolliplot")
rm(list = ls())

humanmethylation450 <- read.csv(file = "references/humanmethylation450_15017482_v1-2.csv",
                                sep = ",", skip = 7)
M450 <- humanmethylation450[,c('Name','CHR','UCSC_RefGene_Name','MAPINFO','UCSC_RefGene_Group','UCSC_CpG_Islands_Name','Relation_to_UCSC_CpG_Island')]
M450 <- M450[M450$UCSC_RefGene_Name!='',]

### Split the same CpG site corresponding to different genes into different rows
### If the same CpG site corresponds to multiple groups in the same gene, only the first group is taken
removedup2 <- function (str1,str2) {
  s1 <- unlist(strsplit(str1,split = ';'))
  if (length(s1) == 1) return(str2)
  s2 <- unlist(strsplit(str2,split = ';'))
  s1.order <- order(s1)
  # r1 <- s1[1]
  r2 <- s2[s1.order[1]]
  for (i in 2:length(s1.order)) {
    if (s1[s1.order[i]] == s1[s1.order[i-1]]) next
    # r1<-paste(r1,s1[i],sep = ';')
    r2<-paste(r2,s2[s1.order[i]],sep = ';')
  }
  return(r2)
}
removedup1 <- function (str1) {
  s1 <- unlist(strsplit(str1,split = ';'))
  if (length(s1) == 1) return(str1)
  s1.order <- order(s1)
  r1 <- s1[s1.order[1]]
  for (i in 2:length(s1.order)) {
    if (s1[s1.order[i]] == s1[s1.order[i-1]]) next
    r1<-paste(r1,s1[s1.order[i]],sep = ';')
  }
  return(r1)
}
M450[c('UCSC_RefGene_Group')] <- mapply(removedup2, M450$UCSC_RefGene_Name,M450$UCSC_RefGene_Group)
M450[c('UCSC_RefGene_Name')] <- mapply(removedup1, M450$UCSC_RefGene_Name)

library(tidyr)
M450 <- M450 %>% separate_rows(UCSC_RefGene_Name,UCSC_RefGene_Group,sep = ';')
# library(tidyverse)
# M450 <- separate_rows(M450, c(UCSC_RefGene_Name,UCSC_RefGene_Group),sep = ';')
M450 <- unique(M450[order(M450$CHR,M450$UCSC_RefGene_Name,M450$MAPINFO),])
M450 <- transform(M450, MAPINFO = as.numeric(MAPINFO))
###---------------------------
generange.color <- data.frame(UCSC_RefGene_Group = unique(M450$UCSC_RefGene_Group))
generange.color$color <- c(1:length(generange.color$UCSC_RefGene_Group))+128

CpGs.dat <- read.csv(file = "./BRCA/figure/TNBC.w.ctrl_DMP_genes.txt", sep = "\t")
CpGs.dat <- transform(CpGs.dat, MAPINFO = as.numeric(MAPINFO))
# CpGs.dat <- na.omit(CpGs.dat)
colnames(CpGs.dat)[which(names(CpGs.dat) == "CHR")] <- "chr"
CpGs.dat  = data.frame("CpG"=rownames(CpGs.dat),CpGs.dat)

SNPs.dat <- read.csv(file = "data/external/breast-tumor/gdac.broadinstitute.org_BRCA.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/trip/SNPs.txt", 
                     sep = "\t")
SNPs.dat <- transform(SNPs.dat, position = as.numeric(position))
SNPs.dat <- na.omit(SNPs.dat)
SNPs.dat = SNPs.dat[,c('gene','chr','position')]%>%group_by_all%>%count

deltaBeta <- 0.3
beta.data <- CpGs.dat[abs(CpGs.dat$deltaBeta)>deltaBeta,]

nCpGs <- length(unique(beta.data$CpG))
gene_chr <- unique(beta.data[c('gene','chr')])
ngene <- nrow(gene_chr)

mincount.cpgs <- 5 ## plot genes at least has this number CpGs
maxcount.cpgs <- 0 ## 0 means no limit
savepath <- paste('deltaBeta',deltaBeta,'/',sep = '')
dir.create(savepath, showWarnings = FALSE)

#gname <- 'DPP6'
#chr <- 7
#gname <- 'ZNF709'
#chr <- 19
for(i in 1:nrow(gene_chr)){
  gname <- gene_chr[i,'gene']
  chr <- gene_chr[i,'chr']

  gdat = beta.data[beta.data$gene==gname&beta.data$chr==chr,]

  CpGs = gdat[,c('CpG','MAPINFO')]%>%group_by_all%>%count
  CpGs.count = sum(CpGs$n)
  if (CpGs.count<mincount.cpgs) next
  if (maxcount.cpgs>0 & CpGs.count>maxcount.cpgs) next
  CpGs = CpGs[order(CpGs$MAPINFO),]
  SNPs = SNPs.dat[SNPs.dat$gene==gname&SNPs.dat$chr==chr,]
  SNPs = SNPs[order(SNPs$position),]
  
  geneinfo = unique(M450[M450$UCSC_RefGene_Name==gname&M450$CHR==chr,])
  # gene.start.end = unique(generange[generange$gene==gname&generange$chr==chr,])
  # gene.start.end <- gene.start.end %>% group_by(gene,chr) %>% summarise(
  #   start = min(start, na.rm = T),
  #   end = max(end, na.rm = T))
  generange <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(generange) <- c('UCSC_RefGene_Group','start', 'end')
  generange[1,] <- c(geneinfo$UCSC_RefGene_Group[1],geneinfo$MAPINFO[1],'')
  for (i in 2:dim(geneinfo)[1]) {
    if (geneinfo$UCSC_RefGene_Group[i]==geneinfo$UCSC_RefGene_Group[i-1]) next
    n = nrow(generange)
    # endpoint = geneinfo$MAPINFO[i-1]
    # generange$end[n] <- endpoint
    # generange[n+1,] <- c(geneinfo$UCSC_RefGene_Group[i],geneinfo$MAPINFO[i],'')
    ## extend the coordinate to adjacent
    endpoint = floor((geneinfo$MAPINFO[i-1]+geneinfo$MAPINFO[i])/2)
    generange$end[n] <- endpoint
    generange[n+1,] <- c(geneinfo$UCSC_RefGene_Group[i],endpoint+1,'')
  }
  generange$end[nrow(generange)] <- geneinfo$MAPINFO[nrow(geneinfo)]
  generange <- transform(generange, start = as.numeric(start))
  generange <- transform(generange, end = as.numeric(end))
  generange <- merge(generange,generange.color,by=c('UCSC_RefGene_Group'))
  generange = generange[order(generange$start),]

  island <- str_split_fixed(geneinfo$UCSC_CpG_Islands_Name, ":", 2)
  island <- unique(str_split_fixed(island[,2], "-", 2))
  island <- data.frame(island)
  colnames(island) <- c('start','end')
  island <- island[!island$start==''|!island$end=='',]
  island <- transform(island, start = as.numeric(start))
  island <- transform(island, end = as.numeric(end))
  island = island[order(island$start),]
  
  ## try to scale the small range group ... not finish
  # feature.start <- min(CpGs$MAPINFO,SNPs$position,generange$start,island$start)
  # feature.end <- max(CpGs$MAPINFO,SNPs$position,generange$end,island$end)
  # feature.minwidth <- round((feature.end-feature.start)/100)
  
  features.gene <- GRanges(paste("chr",chr,sep=''), 
                           IRanges(generange$start,end=generange$end,names=generange$UCSC_RefGene_Group))
  features.gene$fill <- generange$color
  features.gene$featureLayerID <- rep('gene',nrow(generange))
  
  if (nrow(island) > 0) {
    features.island <- GRanges(paste("chr",chr,sep=''), 
                        IRanges(island$start,end=island$end,names=rep('island',nrow(island))))
    features.island$fill <- rep(80,nrow(island))
    features.island$featureLayerID <- rep('island',nrow(island))
    
    features <- c(features.gene,features.island)
  } else {
    features <- c(features.gene)
  }
  names(features) <- c(generange$UCSC_RefGene_Group,rep('CpG-Island',nrow(island)))
  # features$height <- c(unit(1/256, "inches"),rep(unit(1/128, "inches"),length(features)-1))
  features$height <- rep(0.02,length(features))
  
  if (dim(SNPs)[1]>0) {
    lolli <- c(CpGs$MAPINFO,SNPs$position)
    color <- c(rep(3,length(CpGs$MAPINFO)),rep(4,length(SNPs$position)))
    legend <- c(3,4)
    names(legend) <- c('CpG','SNP')
  } else {
    lolli <- CpGs$MAPINFO
    color <- rep(3,length(CpGs$MAPINFO))
    legend <- c(3)
    names(legend) <- c('CpG')
  }
  
  sample.gr <- GRanges(paste("chr",chr,sep=''), IRanges(lolli, width=1))
  sample.gr$color <- color
  sample.gr$score <- sample.int(7, length(sample.gr), replace=TRUE)
  sample.gr$score[1] <- 2.1
  
  png(filename=paste(savepath,'CpGs',CpGs.count,'_SNPs',dim(SNPs)[1],'_',gname,".png",sep = ''),width=800,height=200)
  lolliplot(sample.gr, features,yaxis=FALSE, legend=legend, ylab=gdat$gene[1],cex=1)
  # rescale <- c(.1, .1, .1, .1, .1, .1, .1, .1, .2, .1, .1, .1, .1, .1, .1, .4)
  # lolliplot(sample.gr, features,yaxis=FALSE, legend=legend,
  #           rescale = rescale, ylab=gdat$gene[1],cex=1)
  # lolliplot(sample.gr, features, ranges = GRanges(paste("chr",chr,sep=''), IRanges(min(generange$start), max(generange$end))),
  #           yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  dev.off()
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  # grid.text("title here", x=.5, y=.98, just="top", 
  #           gp=gpar(cex=1.5, fontface="bold"))
}

