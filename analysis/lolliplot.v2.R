# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
library(trackViewer)
library(dplyr)

setwd("./BRCA/figure/TNBC.w.ctrl_DMP_genes_thr_Dbeta0")
rm(list = ls())

CpGs.dat <- read.csv(file = "./TNBC.w.ctrl_DMP_genes_thr_Dbeta0_loc.txt", sep = "\t")
CpGs.dat <- transform(CpGs.dat, CpG.loc = as.numeric(CpG.loc))
CpGs.dat <- na.omit(CpGs.dat)

generange <- read.csv(file = "references/hg37ref/Homo_sapiens.GRCh37.75.range.txt", sep = "\t", skip = 5)
generange.color <- data.frame(type = unique(generange[,'type']))
generange.color$color <- c(1:length(generange.color$type))+128

SNPs.dat <- read.csv(file = "data/external/breast-tumor/gdac.broadinstitute.org_BRCA.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/trip/SNPs.txt", sep = "\t")
SNPs.dat <- transform(SNPs.dat, position = as.numeric(position))
SNPs.dat <- na.omit(SNPs.dat)
SNPs.dat = SNPs.dat[,c('gene','chr','position')]%>%group_by_all%>%count

CpGs.dat <- merge(CpGs.dat,generange[generange$type=='gene',],by=c('gene','chr'))
CpGs.dat <- subset(CpGs.dat,select=-c(type))

CpGs.dat.inside <- CpGs.dat[CpGs.dat$CpG.loc>CpGs.dat$start&CpGs.dat$CpG.loc<CpGs.dat$end,]

deltaBeta <- 0.3
beta.data <- CpGs.dat.inside[abs(CpGs.dat.inside$deltaBeta)>deltaBeta,]

nCpGs <- length(unique(beta.data$CpG))
gene_chr <- unique(beta.data[c('gene','chr')])
ngene <- nrow(gene_chr)

mincount <- 1
savepath <- paste('deltaBeta',deltaBeta,'minCpGs',mincount,'/',sep = '')
dir.create(savepath, showWarnings = FALSE)

#gname <- 'DPP6'
#chr <- 7
#gname <- 'CNTNAP5'
#chr <- 2
for(i in 1:nrow(gene_chr)){
  gname <- gene_chr[i,'gene']
  chr <- gene_chr[i,'chr']

  gdat = beta.data[beta.data$gene==gname&beta.data$chr==chr,]

  CpGs = gdat[,c('CpG','CpG.loc')]%>%group_by_all%>%count
  CpGs.count = sum(CpGs$n)
  if (CpGs.count<mincount) next
  
  geneinfo = unique(generange[generange$gene==gname&generange$chr==chr,])
  geneinfo$width = geneinfo$end-geneinfo$start+1
  thisgene = geneinfo[geneinfo$type=='gene',]
  geneinfo = geneinfo[geneinfo$type!='gene',]
  geneinfo = geneinfo[geneinfo$type!='transcript',]
  geneinfo.type = unique(geneinfo$type)
  
  
  minmallength <- round(thisgene$width*0.01)
  geneinfo$width[geneinfo$width<minmallength] <- minmallength
  
  
  type = geneinfo.type[1]
  rng <- geneinfo[geneinfo$type==type,]
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
    lolli <- c(CpGs$CpG.loc,SNPs$position)
    color <- c(rep(3,length(CpGs$CpG.loc)),rep(4,length(SNPs$position)))
    legend <- c(3,4)
    names(legend) <- c('CpG','SNP')
  } else {
    lolli <- CpGs$CpG.loc
    color <- rep(3,length(CpGs$CpG.loc))
    legend <- c(3)
    names(legend) <- c('CpG')
  }
  
  sample.gr <- GRanges(paste("chr",thisgene$chr,sep=''), IRanges(lolli, width=1))
  sample.gr$color <- color
  sample.gr$score <- sample.int(7, length(sample.gr), replace=TRUE)
  sample.gr$score[1] <- 2.1
  
  png(filename=paste(savepath,'SNPs',dim(SNPs)[1],'_CpGs',CpGs.count,'_',gname,".png",sep = ''),width=800,height=200)
  #lolliplot(sample.gr, features,yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  lolliplot(sample.gr, features, ranges = GRanges(paste("chr",thisgene$chr,sep=''), IRanges(thisgene$start, thisgene$end)),
            yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  dev.off()
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  # grid.text("title here", x=.5, y=.98, just="top", 
  #           gp=gpar(cex=1.5, fontface="bold"))
}

