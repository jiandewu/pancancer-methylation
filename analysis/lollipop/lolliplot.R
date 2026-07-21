# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
library(trackViewer)
library(dplyr)

setwd("./BRCA/figure/TNBC.w.ctrl_DMP_genes_thr_Dbeta0")
rm(list = ls())

hdat <- read.csv(file = "./TNBC.w.ctrl_DMP_genes_thr_Dbeta0_Som_pos_loc_grange_inside_islandin.txt", sep = "\t")

generange <- read.csv(file = "references/hg37ref/Homo_sapiens.GRCh37.75.range.txt", sep = "\t", skip = 5)
generange.color <- data.frame(type = unique(generange[,'type']))
generange.color$color <- c(1:length(generange.color$type))+30

deltaBeta <- 0.3
savepath <- paste('range_deltaBeta',deltaBeta,'/',sep = '')
dir.create(savepath, showWarnings = FALSE)

beta.data <- hdat[abs(hdat$deltaBeta)>deltaBeta,]
nCpGs <- length(unique(beta.data$CpG))
genelist <- unique(beta.data$gene)
ngene <- length(genelist)
#gname <- 'ELFN1'
for (gname in genelist) {
  gdat = beta.data[beta.data$gene==gname,]
  
  geneinfo = unique(generange[generange$gene==gname,])
  geneinfo$width = geneinfo$end-geneinfo$start+1
  thisgene = geneinfo[geneinfo$type=='gene',]
  geneinfo = geneinfo[geneinfo$type!='gene',]
  geneinfo = geneinfo[geneinfo$type!='transcript',]
  geneinfo.type = unique(geneinfo$type)
  
  features <- GRanges(paste("chr",thisgene$chr,sep=''), IRanges(thisgene$start, thisgene$end))
  features$fill <- generange.color[generange.color$type=='transcript','color']
  featureLayerID <- c('transcript')
  
  #geneinfo.type=geneinfo.type[-1]
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
  features$height <- rep(0.01,length(features))
  
  CpGs = gdat[,c('CpG','CpG.loc')]%>%group_by_all%>%count
  SNPs = gdat[,c('gene','SNP')]%>%group_by_all%>%count
  
  lolli <- c(CpGs$CpG.loc,SNPs$SNP)
  color <- c(rep(3,length(CpGs$CpG.loc)),rep(4,length(SNPs$SNP)))
  
  sample.gr <- GRanges(paste("chr",thisgene$chr,sep=''), IRanges(lolli, width=1))
  sample.gr$color <- color
  sample.gr$score <- sample.int(11, length(sample.gr), replace=TRUE)
  sample.gr$score[1] <- 11
  legend <- c(3,4)
  names(legend) <- c('CpG','SNP')
  
  png(filename=paste(savepath,gname,"_lolli.png",sep = ''),width=800,height=200)
  lolliplot(sample.gr, features,yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  lolliplot(sample.gr, features, ranges = GRanges(paste("chr",thisgene$chr,sep=''), IRanges(thisgene$start, thisgene$end)),
            yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  dev.off()
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  # grid.text("title here", x=.5, y=.98, just="top", 
  #           gp=gpar(cex=1.5, fontface="bold"))
}
