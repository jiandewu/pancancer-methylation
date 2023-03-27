rm(list = ls())
# setwd("/external/prjs/Analysis/ColonRectum/MethySurvival")
load("./PRAD/Rdata/step2_filtered_pd_myNorm.Rdata")
load("./BRCA/Rdata/step2_filtered_pd_myNorm(2).Rdata")
load("./LUNG/Rdata/step2_filtered_pd_myNorm.Rdata")
load("../colon-rectum/Rdata/step2_filtered_pd_myNorm.Rdata")

load("./BRCA/Rdata/TNBC.w.ctrl_step2_filtered_pd_myNorm.Rdata")

# path <- "."
# cancername <- "PRAD"
# cancername <- "BRCA"
# prj.name <- "COAD"
setwd(paste(path,cancername,sep = "/"))
dim(myNorm)
dim(pd)

group_list <- pd$group_list
library(ChAMP)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)

head(myDMP$TNBC_to_Ctrl)
# head(myDMP$Tumor_to_Normal)
# head(myDMP$Normal_to_Tumor)

DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=group_list)

df_DMP <- myDMP$TNBC_to_Ctrl
#df_DMP <- myDMP$Tumor_to_Normal
# df_DMP <- myDMP$Normal_to_Tumor


df_DMP=df_DMP[df_DMP$gene!="",]
# logFC_t <- 0.45
deltaBeta_t <- 0
deltaBeta_t <- 0.5
P.Value_t <- 0.05 #10^-15
# df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$logFC) > logFC_t,
#                         ifelse(df_DMP$logFC > logFC_t ,'UP','DOWN'),'NOT') 
df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$deltaBeta) > deltaBeta_t,
                        ifelse(df_DMP$deltaBeta > deltaBeta_t ,'UP','DOWN'),'NOT') 
table(df_DMP$change) 
write.table(df_DMP, file = paste('./figure/',prj.name,'_DMP_genes.txt',sep = ""), sep="\t", quote = F, row.names=T)
df_DMP_threshoff <- df_DMP[df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$deltaBeta) > deltaBeta_t,]

write.table(df_DMP_threshoff, file = paste('./figure/',prj.name,'_DMP_genes_thr_Dbeta0.txt',sep = ""), sep="\t", quote = F, row.names=T)
write.table(df_DMP_threshoff, file = paste('./figure/',prj.name,'_DMP_genes_thr.txt',sep = ""), sep="\t", quote = F, row.names=T)

save(df_DMP,prj.name,path,cancername,file = paste('./Rdata/',prj.name,'_step3.df_DMP.Rdata',sep = ""))




rm(list = ls())
load("./BRCA/Rdata/TNBC.w.ctrl_step2_filtered_pd_myNorm.Rdata")
load("./BRCA/Rdata/TNBC.w.ctrl_step3.df_DMP.Rdata")
setwd(paste(path,cancername,sep = "/"))

library(tibble)
dat  = rownames_to_column(df_DMP)
library(ggplot2)
# p <- ggplot(data = dat, 
#             aes(x = logFC, 
#                 y = -log10(adj.P.Val))) +
#   geom_point(alpha=0.4, size=3.5, 
#              aes(color=change)) +
#   ylab("-log10(Pvalue)")+
#   scale_color_manual(values=c("blue", "grey","red"))+
#   geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
#   geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
#   theme_bw()
p <- ggplot(data = dat, 
            aes(x = deltaBeta, 
                y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=2.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-deltaBeta_t,deltaBeta_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

png(file=paste('./figure/',prj.name,'_volcano.png',sep = ""),width=900, height=900)
p
dev.off()

for_label <- dat%>% head(10)
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = rowname),
    data = for_label,
    color="black"
  )
volcano_plot

png(file=paste('./figure/',prj.name,'_volcano2.png',sep = ""),width=900, height=900)
volcano_plot
dev.off()

deltaBeta <- 0.3
cg <-  rownames(df_DMP[df_DMP$change != "NOT" ,])
cg <-  rownames(df_DMP[abs(df_DMP$deltaBeta)>deltaBeta ,])

plot_matrix <- myNorm[cg,]
annotation_col <- data.frame(Sample=pd$group_list) 
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Sample = c(Ctrl="#4DAF4A", TNBC="#E41A1C"))

library(pheatmap)
pheatmap(plot_matrix,show_colnames = F,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors,show_rownames = F)


png(file=paste('./figure/',prj.name,'_heatmap_Dbeta0.png',sep = ""),width=900, height=900)
pheatmap(plot_matrix,show_colnames = F,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors,show_rownames = F)
dev.off()
# DMRs
myDMR <- champ.DMR(beta=myNorm,pheno=group_list,method="Bumphunter")
head(myDMR$BumphunterDMR)
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=group_list)

library(dplyr)
