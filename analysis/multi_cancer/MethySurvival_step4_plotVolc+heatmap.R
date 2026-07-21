
rm(list=setdiff(ls(), c('prj.name','cancer.top','workingpath','cancername','pd','myNorm','DMP','DMRBumphunter')))
if (!exists('pd') | !exists('myNorm')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step2_filterNAClinNorm_pd.Rdata',sep = ''))
}
if (!exists('DMRBumphunter')) {
  DMRBumphunter = data.table::fread(paste(workingpath,'/figure/',prj.name,'_DMR.Bumphunter.txt',sep = ""), data.table = F)
}
if (!exists('DMP')) {
  DMP = data.table::fread(paste(workingpath,'/figure/',prj.name,'_DMP.DMR.txt',sep = ""), data.table = F)
}
rownames(DMP) <- DMP[, 1]
setwd(workingpath)

# library(ChAMP)
# DMP.GUI(DMP=DMP,beta=myNorm,pheno=pd$Sample_Group)
# DMR.GUI(DMR=DMR.Bumphunter,beta=myNorm,pheno=pd$Sample_Group)

DMP=DMP[DMP$gene!="",]

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','*********** DMP(gene): adj p-value < 0.05, |deltaBeta| > 0\n')
print(table(DMP$DMPrefCtrl))
cat('\n')
print(table(cut(DMP$adj.P.Val,breaks=c(0,10^-45,10^-35,10^-25,10^-15,10^-5,0.005,0.01,0.05,1))))
cat('\n','*********** DMP(gene) integrate DMR: adj p-value < 0.05, |deltaBeta| > 0; DMR.p-value < 0.05\n')
print(table(DMP$DMRDMPrefCtrl))
sink()

write.table(DMP[DMP$adj.P.Val<0.05,],
            file = paste('./figure/',prj.name,'_DMP_thr_p0.05_dB0.txt',sep = ""),
            sep="\t", quote=F, row.names=F)

# logFC_t <- 0.45
deltaBeta_t <- 0.3
P.Value_t <- 10^-15
P.Value_t <- 0.005 # TNBC
DMP$DMPrefCtrl <- ifelse(DMP$adj.P.Val<P.Value_t & abs(DMP$deltaBeta) > deltaBeta_t, 
                         ifelse(DMP$logFC<0,'Hypo','Hyper'),'NOT') 
DMP <- DMP %>% mutate(DMRDMPrefCtrl = case_when(is.na(DMRrefCtrl) & DMPrefCtrl == 'NOT' ~ 'NOT BOTH',
                                               is.na(DMRrefCtrl) ~ 'DMP ONLY',
                                               DMPrefCtrl == 'NOT' ~ 'DMR ONLY',
                                               DMPrefCtrl != DMRrefCtrl ~ 'INV MUTE',
                                               TRUE ~ paste(DMPrefCtrl,'BOTH')))

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','*********** DMP(gene): adj p-value < ',P.Value_t,' |deltaBeta| > ',deltaBeta_t,'\n')
print(table(DMP$DMPrefCtrl))
cat('\n','*********** DMP(gene) integrate DMR: adj p-value < ',P.Value_t,' |deltaBeta| > ',deltaBeta_t,'DMR.p-value < 0.05\n')
print(table(DMP$DMRDMPrefCtrl))
sink()

DMP_cutoff <- DMP[DMP$DMRDMPrefCtrl=='INV MUTE'|DMP$DMRDMPrefCtrl=='Hypo BOTH'|DMP$DMRDMPrefCtrl=='Hyper BOTH',]

write.table(DMP_cutoff,
            file = paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_DMR0.05.txt',sep = ""),
            sep="\t", quote = F, row.names=F)

save(deltaBeta_t,P.Value_t,file = paste('./Rdata/',prj.name,'_step4_DMPthreshoff.Rdata',sep = ""))


library(tibble)
library(ggplot2)
library(latex2exp)
p <- ggplot(data = DMP,aes(x = -deltaBeta,y = -log10(adj.P.Val))) + 
  geom_point(alpha=0.4, size=0.5,aes(color=DMPrefCtrl)) + 
  ylab(TeX("$-\\log_{10}(adj.Pvalue)$")) + labs(color = "Significance") + scale_x_continuous(name=TeX("$\\Delta\\beta$")) + 
  # scale_y_continuous(name=TeX("$-\\log_{10}(adj.Pvalue)$")) +
  scale_color_manual(values=c("red","blue", "grey")) + 
  geom_vline(xintercept=c(-deltaBeta_t,deltaBeta_t),lty=4,col="black",lwd=0.8) + 
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) + 
  theme_bw() + ggtitle(paste(cancername,prj.name,sep = ' - '))

  
png(file=paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_volcano.png',sep = ""))
# dev.list()
dev.set(which = 2)
p
dev.copy(which = 4)
# dev.cur()
dev.off()

p <- ggplot(data = DMP,aes(x = -deltaBeta,y = -log10(adj.P.Val))) + 
  geom_point(alpha=0.8, size=0.5,aes(color=DMRDMPrefCtrl)) + 
  ylab(TeX("$-\\log_{10}(adj.Pvalue)$")) + labs(color = "Significance") + scale_x_continuous(name=TeX("$\\Delta\\beta$")) + 
  # scale_y_continuous(name=TeX("$-\\log_{10}(adj,Pvalue)$")) +
  # scale_color_manual(values=c("blue","green", "grey","red"),labels=c('Hypo DMR', 'DMP x DMR', 'NOT', 'Hyper DMR')) + 
  scale_color_manual(values=c("azure4","darkorchid1","red","blue","green","grey")) + 
  geom_vline(xintercept=c(-deltaBeta_t,deltaBeta_t),lty=4,col="black",lwd=0.8) + 
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) + 
  theme_bw() + ggtitle(paste(cancername,prj.name,'DMP & DMR',sep = ' - '))
for_label <- DMP[DMP$DMRDMPrefCtrl=='INV MUTE',]
p <- p +
  geom_point(size = 2, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = CpG),
    data = for_label,
    color="black"
  )

png(file=paste('./figure/',prj.name,'_DMP.DMR_thr_p',P.Value_t,'_dB',deltaBeta_t,'_volcano.png',sep = ""))
dev.set(which = 2)
p
dev.copy(which = 4)
dev.off()


cg <- DMP_cutoff[,c('CpG')]
j <- c(which(pd$Sample_Group == "Tumor"), which(pd$Sample_Group == "Ctrl"))

plot_matrix <- myNorm[cg,j]
annotation_col <- data.frame(Sample=pd$Sample_Group[j]) 
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors = list(Sample = c(Ctrl="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
png(file=paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_heatmap.png',sep = ""))
dev.set(which = 2)
pheatmap(plot_matrix,show_colnames = F,
         annotation_col = annotation_col,cluster_cols = F,
         border_color=NA, main = paste(cancername,' - ',prj.name,' (',length(cg),' CpGs)',sep = ''),
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors,show_rownames = F)
dev.copy(which = 4)
dev.off()

plot(-log10(DMP$DMPfwer), -log10(DMP$adj.P.Val), main="Scatterplot Example", 
     xlab="DMR ", ylab="DMP", pch=19)
