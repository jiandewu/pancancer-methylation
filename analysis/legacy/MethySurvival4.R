rm(list = ls())
load("./BRCA/Rdata/TNBC.w.ctrl_step2_filtered_pd_myNorm.Rdata")
load("./BRCA/Rdata/TNBC.w.ctrl_step3.df_DMP.Rdata")
setwd(paste(path,cancername,sep = "/"))

cg <- rownames(df_DMP[df_DMP$change != "NOT",])
myNorm_tumor <- myNorm[cg,]

suv_dat <- data.table::fread("./raw_data/TCGA-BRCA.survival.tsv")
library(stringr)
suv_dat$sample = str_sub(suv_dat$sample,1,15)
suv_dat <- suv_dat[suv_dat$sample %in% colnames(myNorm_tumor),]
suv_dat <- suv_dat[str_sub(suv_dat$sample,14,15)=="01",] 
colnames(pd)[1] <- "sampleID"
suv_dat = merge(pd,suv_dat,by.x = "sampleID",by.y = "sample")
suv_dat$months <- suv_dat$OS.time/30
myNorm_tumor <- myNorm_tumor[,suv_dat$sample]
identical(colnames(myNorm_tumor),suv_dat$sample)

library(survival)
# cox method
cox_results <- apply(myNorm_tumor, 1, function(x){
  #x <- myNorm_tumor[1,]
  suv_dat$cox_group <- ifelse(x>mean(x),"High","Low")
  res <- coxph(Surv(months, OS)~cox_group, data=suv_dat)
  beta <- coef(res)
  se <- sqrt(diag(vcov(res)))
  p.val <- 1 - pchisq((beta/se)^2, 1)
})
table(cox_results<0.05) 

# log-rank method
log_rank_p <- apply(myNorm_tumor , 1 , function(gene){
  # gene=exprSet[1,]
  suv_dat$log_rank_p_group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(months, OS)~log_rank_p_group,data=suv_dat)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
# log_rank_p=sort(log_rank_p)
table(log_rank_p<0.05) 



surv_gene <- names(sort(log_rank_p))[1:2] 
choose_matrix <- myNorm[surv_gene,]
rownames(choose_matrix) <- df_DMP[surv_gene,"gene"]
annotation_col <- data.frame(Sample=pd$group_list) 
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
pheatmap(choose_matrix,show_colnames = F,
         annotation_col = annotation_col,
         border_color=NA, cluster_cols = T,
         color = colorRampPalette(colors = c("white","navy"))(50),
         #filename = "./figure/Survival_heatmap.png",
         annotation_colors = ann_colors)

gs=head(surv_gene,12)
exprSet = myNorm_tumor
meta = suv_dat
# BiocManager::install("survminer")
library(survminer)
splots <- lapply(gs, function(g){
  meta$group=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(months, OS)~group, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta, risk.table = TRUE, 
             legend.title = g, legend.labs = c("High", "Low"),)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 1, risk.table.height = 0.4)
# Arrange and save into png file
res <- arrange_ggsurvplots(splots, print = FALSE,  
                           ncol = 2, nrow = 1, risk.table.height = 0.33)
ggsave("./figure/survival_prob.png", width = 16, height = 12, units = "in",scale = 1, res)




library(data.table)
