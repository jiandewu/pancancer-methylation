# ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/trackViewer/inst/doc/lollipopPlot.html
library(trackViewer)
library(dplyr)

setwd("./BRCA/figure/TNBC.w.ctrl_DMP_genes_thr_Dbeta0")
rm(list = ls())

hdat <- read.csv(file = "./TNBC.w.ctrl_DMP_genes_thr_Dbeta0_Som_pos_loc_grange_inside_islandin.txt", sep = "\t")

generange <- read.csv(file = "references/hg37ref/Homo_sapiens.GRCh37.75.range.txt", sep = "\t", skip = 5)
generange.color <- data.frame(type = unique(generange[,'type']))
generange.color$color <- c(1:length(generange.color$type))+60

deltaBeta <- 0.3
savepath <- paste('range_deltaBeta',deltaBeta,'/',sep = '')
dir.create(savepath, showWarnings = FALSE)

beta.data <- hdat[abs(hdat$deltaBeta)>deltaBeta,]
nCpGs <- length(unique(beta.data$CpG))

gene_chr <- unique(beta.data[c('gene','chr')])
#genelist <- unique(beta.data$gene)
ngene <- nrow(gene_chr)
#gname <- 'PCDHA2'
#chr <- 5
for(i in 1:nrow(gene_chr)){
  gname <- gene_chr[i,'gene']
  chr <- gene_chr[i,'chr']

  gdat = beta.data[beta.data$gene==gname&beta.data$chr==chr,]
  
  geneinfo = unique(generange[generange$gene==gname&generange$chr==chr,])
  geneinfo$width = geneinfo$end-geneinfo$start+1
  thisgene = geneinfo[geneinfo$type=='gene',]
  geneinfo = geneinfo[geneinfo$type!='gene',]
  geneinfo = geneinfo[geneinfo$type!='transcript',]
  geneinfo.type = unique(geneinfo$type)
  
  
  minmallength <- round(thisgene$width*0.005)
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
  features$height <- rep(0.01,length(features))
  
  CpGs = gdat[,c('CpG','CpG.loc')]%>%group_by_all%>%count
  SNPs = gdat[,c('gene','SNP')]%>%group_by_all%>%count
  
  lolli <- c(CpGs$CpG.loc,SNPs$SNP)
  color <- c(rep(3,length(CpGs$CpG.loc)),rep(4,length(SNPs$SNP)))
  
  sample.gr <- GRanges(paste("chr",thisgene$chr,sep=''), IRanges(lolli, width=1))
  sample.gr$color <- color
  sample.gr$score <- sample.int(7, length(sample.gr), replace=TRUE)
  sample.gr$score[1] <- 2.1
  legend <- c(3,4)
  names(legend) <- c('CpG','SNP')
  
  png(filename=paste(savepath,gname,"_lolli.png",sep = ''),width=800,height=200)
  #lolliplot(sample.gr, features,yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  lolliplot(sample.gr, features, ranges = GRanges(paste("chr",thisgene$chr,sep=''), IRanges(thisgene$start, thisgene$end)),
            yaxis=FALSE, legend=legend, ylab=gdat$gene[1], rescale = TRUE, cex=1)
  dev.off()
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  # grid.text("title here", x=.5, y=.98, just="top", 
  #           gp=gpar(cex=1.5, fontface="bold"))
}
