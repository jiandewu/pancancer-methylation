
rm(list=setdiff(ls(), c('prj.name','workingpath','cancer.top','cancername','deltaBeta_t','P.Value_t','suv_dat','myNorm_signif','pd','myNorm','DMP_cutoff_survival')))
if (!exists('myNorm') | !exists('pd')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step2_filterNAClinNorm_pd.Rdata',sep = ''))
}
if (!exists('suv_dat') | !exists('myNorm_signif')) {
  load(paste('./Rdata/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_step5_DMP_cox_logrank.Rdata',sep = ""))
}
if (!exists('DMP_cutoff_survival')) {
  DMP_cutoff_survival = data.table::fread(paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_DMR0.05.txt',sep = ""), data.table = F)
}
rownames(DMP_cutoff_survival) <- DMP_cutoff_survival[, 1]

cutoff <- 0.05
df_DMP_cutoff <- DMP_cutoff_survival[DMP_cutoff_survival$log_rank_p<cutoff & DMP_cutoff_survival$cox_results<cutoff,]
df_DMP_cutoff <- DMP_cutoff_survival[DMP_cutoff_survival$log_rank_p<cutoff,]

df_DMP_cutoff <- df_DMP_cutoff[order(df_DMP_cutoff$log_rank_p),]
numplot <- min(dim(df_DMP_cutoff)[1],16)

surv_CpGs <- df_DMP_cutoff$CpG[1:numplot] 
choose_matrix <- myNorm[surv_CpGs,]
rownames(choose_matrix) <- paste(surv_CpGs,' (',df_DMP_cutoff[surv_CpGs,"gene"],')',sep='')
annotation_col <- data.frame(Sample=pd$Sample_Group) 
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Sample = c(Ctrl="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
png(file=paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_logrankcox',cutoff,'_',numplot,'_heatmap.png',sep = ""))
# dev.new()
# dev.list()
dev.set(which = 2)
pheatmap(choose_matrix,show_colnames = F,
         annotation_col = annotation_col,
         border_color=NA, cluster_cols = T,
         color = colorRampPalette(colors = c("white","navy"))(50),
         main = paste(cancername,prj.name,sep = ' - '),
         annotation_colors = ann_colors)
dev.copy(which = 4)
# dev.cur()
dev.off()

gs=head(surv_CpGs,numplot)
exprSet = myNorm_signif
meta = suv_dat
# BiocManager::install("survminer")
library(survminer)
splots <- lapply(gs, function(g){
  meta$group=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(months, OS)~group, data=meta)
  ggsurv <- ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE, pval.method = FALSE,
                       # main = paste(g,' (',DMP_cutoff_survival[g,"gene"],')',sep = ''),
                       # font.main = c(14, "bold", "darkblue"),
                       ylab = paste(g,' (',DMP_cutoff_survival[g,"gene"],')',sep = ''),
                       font.y = c(12, "bold", "black"),conf.int = TRUE,
                       xlab = 'Time (Month)',
                       font.x = c(12, "bold", "black"),
                       surv.median.line='hv',#risk.table.col="strata",
                       legend.title = '',
                       font.legend = c(12, "bold", "black"),
                       legend.labs = c("Hyper.", "Hypo."),)
  ggsurv %+% theme_survminer(
    font.legend = c(12, "plain", "black")
  )
  ggsurv$table <- ggsurv$table + theme_cleantable()
  ggsurv
}) 

ncol = min(4,numplot)
nrow = ceiling(numplot/ncol)
# Arrange and save into png file
res <- arrange_ggsurvplots(splots, print = FALSE,  
                           title = paste(cancername,' (',prj.name,') - Hyper. vs Hypo. (Top ',numplot,' of ',dim(df_DMP_cutoff)[1],' significance CpGs)',sep = ''),  
                           font.title = c(14, "bold", "black"),
                           ncol = ncol, nrow = nrow, risk.table.height = 0.3)
ggsave(paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_logrankcox',cutoff,'_',numplot,'_survival1.png',sep = ""), width = ncol*4, height = nrow*3.5, units = "in",scale = 1, res)

