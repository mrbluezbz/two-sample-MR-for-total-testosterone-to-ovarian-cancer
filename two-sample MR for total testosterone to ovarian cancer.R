remove(list = ls())
########################################
# Date 2022-2-07
# Author 本本
# Data  TwoSampleMR
# TwoSampleMR process
######################################

#安装TwoSampleMR包
install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#设置工作路径
setwd("F:")
getwd()
t <-extract_instruments(outcomes='ieu-b-4864',access_token = NULL) 
head(t) 
t_exp = extract_instruments(
  outcomes='ieu-b-4864',
  clump=TRUE, r2=0.01,
  kb=5000,access_token = NULL
)
ov_out = extract_outcome_data(
  snps=t_exp$SNP,
  outcomes='ieu-a-1120',
  proxies = FALSE,
  maf_threshold = 0.01,
  access_token = NULL
)
mydata = harmonise_data(
  exposure_dat=t_exp,
  outcome_dat=ov_out,
  action= 2
)
res = mr(mydata)
het <- mr_heterogeneity(mydata)
mr(mydata,method_list=c('mr_ivw_mre'))#使用随机效应模型
pleio <- mr_pleiotropy_test(mydata)
single <- mr_leaveoneout(mydata)
single_plot = single[c(1:50),]
mr_leaveoneout_plot(single_plot)
mr_scatter_plot(res,mydata)
res_single <- mr_singlesnp(mydata)
mr_forest_plot(res_single)
mr_funnel_plot(res_single)
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
library(MRPRESSO)
cols = c("beta.exposure","beta.outcome",
         "se.outcome","se.exposure",
         "pval.outcome","pval.exposure")
summarystats = mydata[,colnames(mydata)%in%cols]
mrpresso = mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = summarystats, NbDistribution = 1000,  
                     SignifThreshold = 0.05)
gen_res = generate_odds_ratios(res)
write.csv(gen_res,file = "F:",
          row.names = T,quote = F)
write.csv(het,file = "F:",
          row.names = T,quote = F)
write.csv(pleio,file = "F:",
          row.names = T,quote = F)
dat = mydata
dat$samplesize.exposure = 199569
r2 = directionality_test(dat)
f = as.data.frame((199569-135-1)/135)*(r2$snp_r2.exposure/(1-r2$snp_r2.exposure))
write.csv(r2,file = "F:",
          row.names = T,quote = F)
write.csv(f,file = "F:",
          row.names = T,quote = F)









