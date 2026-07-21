## https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

library(data.table)
# library(rlist)

rm(list = ls())
cancernames <- c("BRCA","COAD","LUAD","PRAD")

prj.name <- "Tumor.vs.Ctrl"
path <- "./multi_cancer/manus"
setwd(path)

deltaBeta_t <- 0.3
P.Value_t <- 10^-15

genes <- list()
CpGs <- list()
genes.surv <- list()
CpGs.surv <- list()
# cancername = "BRCA"
cutoff <- 0.05
for (cancername in cancernames) {
  file = paste('../',cancername,'/figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_DMR0.05.txt',sep = "")
  DMPresult = data.table::fread(file, data.table = F)
   
  CpGs <- c(CpGs,list(DMPresult$CpG))
  names(CpGs)[length(CpGs)] <- cancername
  genes <- c(genes,list(unique(DMPresult$gene)))
  names(genes)[length(genes)] <- cancername
  
  df_DMP_cutoff <- DMPresult[DMPresult$log_rank_p<cutoff & DMPresult$cox_results<cutoff,]
  CpGs.surv <- c(CpGs.surv,list(df_DMP_cutoff$CpG))
  names(CpGs.surv)[length(CpGs.surv)] <- cancername
  genes.surv <- c(genes.surv,list(unique(df_DMP_cutoff$gene)))
  names(genes.surv)[length(genes.surv)] <- cancername
}

library(gplots)
genes.table <- venn(genes, show.plot=FALSE)
genes.table <- as.matrix(genes.table)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  # library(VennDiagram)
  # grid.newpage()
  venn.diagram(x,
               # category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4"),
               filename = NULL,
               # Circles
               lwd = 2,
               lty = 'blank',
               fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
               # Numbers
               cex = .9,
               fontface = "italic",
               # Set names
               cat.cex = 1,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.dist = c(0.055, 0.055, 0.1, 0.1)
  )
  # grid.draw(venn_object)
}

## https://statisticsglobe.com/add-title-venn-diagram-r
library("gridExtra") 

venn_genes <- display_venn(genes)
grid.arrange(gTree(children = venn_genes), # Add title & subtitle
             # bottom = "My subtitle",
             top = "Genes")

png(file=paste(prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_venn_genes.png',sep = ""))
dev.set(which = 2)
dev.copy(which = 4)
dev.off()

venn_CpGs <- display_venn(CpGs)
grid.arrange(gTree(children = venn_CpGs), # Add title & subtitle
             # bottom = "My subtitle",
             top = "CpGs")

png(file=paste(prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_venn_CpGs.png',sep = ""))
dev.set(which = 2)
dev.copy(which = 4)
dev.off()


venn_genes <- display_venn(genes.surv)
grid.arrange(gTree(children = venn_genes), # Add title & subtitle
             # bottom = "My subtitle",
             top = "Survival Genes")

png(file=paste(prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_venn_genes.surv.png',sep = ""))
dev.set(which = 2)
dev.copy(which = 4)
dev.off()

venn_CpGs <- display_venn(CpGs.surv)
grid.arrange(gTree(children = venn_CpGs), # Add title & subtitle
             # bottom = "My subtitle",
             top = "Survival CpGs")

png(file=paste(prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_venn_CpGs.surv.png',sep = ""))
dev.set(which = 2)
dev.copy(which = 4)
dev.off()

