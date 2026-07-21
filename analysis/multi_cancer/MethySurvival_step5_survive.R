
rm(list=setdiff(ls(), c('prj.name','cancer.top','workingpath','cancername','deltaBeta_t','P.Value_t','pd','myNorm','DMP_cutoff')))
if (!exists('pd') | !exists('myNorm')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step2_filterNAClinNorm_pd.Rdata',sep = ''))
}
if (!exists('deltaBeta_t') | !exists('P.Value_t')) {
  load(paste(workingpath,'/Rdata/',prj.name,'_step4_DMPthreshoff.Rdata',sep = ''))
}
if (!exists('DMP_cutoff')) {
  DMP_cutoff = data.table::fread(paste(workingpath,'/figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_DMR0.05.txt',sep = ""), data.table = F)
}
rownames(DMP_cutoff) <- DMP_cutoff[, 1]

setwd(workingpath)

cg <- DMP_cutoff[,c('CpG')]
myNorm_signif <- myNorm[cg,]

# suv_dat <- data.table::fread(paste("../../",cancername,"/raw_data/AllDataHub/TCGA-",cancername,".survival.tsv",sep = ""))
suv_dat <- data.table::fread(paste("../../",cancer.top,"/raw_data/survival_",cancer.top,"_survival.txt",sep = ""))

library(stringr)
suv_dat$sample = str_sub(suv_dat$sample,1,15)
suv_dat <- suv_dat[suv_dat$sample %in% colnames(myNorm_signif),]
suv_dat <- suv_dat[str_sub(suv_dat$sample,14,15)=="01",] 
colnames(pd)[1] <- "sampleID"
suv_dat = merge(pd,suv_dat,by.x = "sampleID",by.y = "sample")
suv_dat$months <- suv_dat$OS.time/30

myNorm_signif <- myNorm_signif[,suv_dat$sample]
identical(colnames(myNorm_signif),suv_dat$sample)

library(survival)
# cox method
cox_results <- apply(myNorm_signif, 1, function(x){
  #x <- myNorm_signif[1,]
  suv_dat$cox_group <- ifelse(x>mean(x),"High","Low")
  res <- coxph(Surv(months, OS)~cox_group, data=suv_dat)
  beta <- coef(res)
  se <- sqrt(diag(vcov(res)))
  p.val <- 1 - pchisq((beta/se)^2, 1)
})

# log-rank method
log_rank_p <- apply(myNorm_signif , 1 , function(gene){
  # gene=exprSet[1,]
  suv_dat$log_rank_p_group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(months, OS)~log_rank_p_group,data=suv_dat)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
# log_rank_p=sort(log_rank_p)

DMP_cutoff_survival <- merge(DMP_cutoff,log_rank_p,by.x = 'CpG',by.y='row.names')
colnames(DMP_cutoff_survival)[ncol(DMP_cutoff_survival)] <- 'log_rank_p'
DMP_cutoff_survival <- merge(DMP_cutoff_survival,cox_results,by.x = 'CpG',by.y='row.names')
colnames(DMP_cutoff_survival)[ncol(DMP_cutoff_survival)] <- 'cox_results'

sink(paste('./figure/',cancername,'_',prj.name,'_results.txt',sep = ""), append = TRUE)
cat('\n','*********** cox_results','\n')
print(table(cut(DMP_cutoff_survival$cox_results,breaks=c(0,0.001,0.005,0.01,0.05,1))))
cat('\n','*********** log_rank_p','\n')
print(table(cut(DMP_cutoff_survival$log_rank_p,breaks=c(0,0.001,0.005,0.01,0.05,1))))
cat('\n','*********** log_rank_p integrate cox_results','\n')
print(table(cut(pmax(DMP_cutoff_survival$log_rank_p,DMP_cutoff_survival$cox_results),breaks=c(0,0.001,0.005,0.01,0.05,1))))
sink()

write.table(DMP_cutoff_survival,
            file = paste('./figure/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_DMR0.05.txt',sep = ""),
            sep="\t", quote = F, row.names=F)

save(suv_dat,myNorm_signif,file=paste('./Rdata/',prj.name,'_DMP_thr_p',P.Value_t,'_dB',deltaBeta_t,'_step5_DMP_cox_logrank.Rdata',sep = ""))
