# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
library(trackViewer)
library(dplyr)

options(dplyr.summarise.inform = FALSE)
setwd("./BRCA/figure/White.TNBC.vs.ctrl_lolliplot")
rm(list = ls())

humanmethylation450 <- read.csv(file = "references/humanmethylation450_15017482_v1-2.csv",
                                sep = ",", skip = 7)
M450 <- humanmethylation450[,c('Name','CHR','UCSC_RefGene_Name','MAPINFO','UCSC_RefGene_Group','UCSC_CpG_Islands_Name','Relation_to_UCSC_CpG_Island')]

M450 <- M450[M450$UCSC_RefGene_Name!='',]
library(tidyr)
M450.split <- M450 %>% separate_rows(UCSC_RefGene_Name,UCSC_RefGene_Group,sep = ';')
# library(tidyverse)
# M450.split <- separate_rows(M450, c(UCSC_RefGene_Name,UCSC_RefGene_Group),sep = ';')
M450.split <- unique(M450.split[order(M450.split$CHR,M450.split$UCSC_RefGene_Name,M450.split$MAPINFO),])


CpGs.dat <- read.csv(file = "./BRCA/figure/TNBC.w.ctrl_DMP_genes.txt", sep = "\t")
CpGs.dat <- transform(CpGs.dat, MAPINFO = as.numeric(MAPINFO))
# CpGs.dat <- na.omit(CpGs.dat)
colnames(CpGs.dat)[which(names(CpGs.dat) == "CHR")] <- "chr"
CpGs.dat  = data.frame("CpG"=rownames(CpGs.dat),CpGs.dat)


generange <- read.csv(file = "references/hg37ref/Homo_sapiens.GRCh37.75.range.txt", 
                      sep = "\t", skip = 5)
generange.color <- data.frame(type = unique(generange[,'type']))
generange.color$color <- c(1:length(generange.color$type))+128

SNPs.dat <- read.csv(file = "data/external/breast-tumor/gdac.broadinstitute.org_BRCA.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/trip/SNPs.txt", 
                     sep = "\t")
SNPs.dat <- transform(SNPs.dat, position = as.numeric(position))
SNPs.dat <- na.omit(SNPs.dat)
SNPs.dat = SNPs.dat[,c('gene','chr','position')]%>%group_by_all%>%count

CpGs.dat <- merge(CpGs.dat,generange[generange$type=='gene',],by=c('gene','chr'))
CpGs.dat <- subset(CpGs.dat,select=-c(type))

CpGs.dat.inside <- CpGs.dat[CpGs.dat$MAPINFO>CpGs.dat$start&CpGs.dat$MAPINFO<CpGs.dat$end,]

deltaBeta <- 0.3
beta.data <- CpGs.dat.inside[abs(CpGs.dat.inside$deltaBeta)>deltaBeta,]

nCpGs <- length(unique(beta.data$CpG))
gene_chr <- unique(beta.data[c('gene','chr')])
ngene <- nrow(gene_chr)

mincount <- 5
maxcount <- -1
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
  if (CpGs.count<mincount) next
  if (maxcount>-1 & CpGs.count>maxcount) next
  
  geneinfo = unique(generange[generange$gene==gname&generange$chr==chr,])
  geneinfo$width = geneinfo$end-geneinfo$start+1
  thisgene = geneinfo[geneinfo$type=='gene',]
  
  thisgene <- thisgene %>% group_by(gene,chr,type) %>% summarise(
    start = min(start, na.rm = T),
    end = max(end, na.rm = T))
  thisgene$width = thisgene$end-thisgene$start+1
  
  geneinfo = geneinfo[geneinfo$type!='gene',]
  geneinfo = geneinfo[geneinfo$type!='transcript',]
  geneinfo.type = unique(geneinfo$type)
  
  
  minmallength <- round(thisgene$width*0.01)
  geneinfo$width[geneinfo$width<minmallength] <- minmallength
  
  
  type = geneinfo.type[1]
  rng <- geneinfo[geneinfo$type==type,]
  # rng = rng[,c('start','width')]%>%group_by_all%>%count
  
  features <- GRanges(paste("chr",thisgene$chr,sep=''), 
                      IRanges(rng$start,width=rng$width,names=rep(type,nrow(rng))))
  features$fill <- rep(generange.color[generange.color$type==type,'color'],nrow(rng))
  featureLayerID <- rep(type,nrow(rng))
  
  geneinfo.type=geneinfo.type[-1]
  #type = 'exon'
  for (type in geneinfo.type) {
    if (type=='gene' || type=='transcript') next
    rng <- geneinfo[geneinfo$type==type,]
    feat <- GRanges(paste("chr",thisgene$chr,sep=''), 
                    IRanges(rng$start,
                            width=rng$width,
                            names=rep(type,nrow(rng))))
    
    feat$fill <- rep(generange.color[generange.color$type==type,'color'],nrow(rng))
    features <- c(features,feat)
    featureLayerID <- c(featureLayerID,rep(type,nrow(rng)))
  }
  
  features$featureLayerID <- featureLayerID
  names(features) <- features$featureLayerID
  
  # features$height <- c(unit(1/256, "inches"),rep(unit(1/128, "inches"),length(features)-1))
  features$height <- rep(0.02,length(features))
  
  SNPs = SNPs.dat[SNPs.dat$gene==gname&SNPs.dat$chr==chr,]
  
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
  
  sample.gr <- GRanges(paste("chr",thisgene$chr,sep=''), IRanges(lolli, width=1))
  sample.gr$color <- color
  sample.gr$score <- sample.int(7, length(sample.gr), replace=TRUE)
  sample.gr$score[1] <- 2.1
  
  png(filename=paste(savepath,'CpGs',CpGs.count,'_SNPs',dim(SNPs)[1],'_',gname,".png",sep = ''),width=800,height=200)
  lolliplot(sample.gr, features, ranges = GRanges(paste("chr",thisgene$chr,sep=''), IRanges(thisgene$start, thisgene$end)),
            yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  dev.off()
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  # grid.text("title here", x=.5, y=.98, just="top", 
  #           gp=gpar(cex=1.5, fontface="bold"))
}

