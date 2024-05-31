#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
vcf_file_path <- "C:/Users/18002/Desktop/SA2ukb-d-G6_SLEEPAPNO.vcf/SA2ukb-d-G6_SLEEPAPNO.vcf"
vcf <- readVcf(vcf_file_path, genome="hg19")

devtools::install_github("MRCIEU/TwoSampleMR")
install.packages("ieugwasr")
library(TwoSampleMR)
library(ieugwasr)
#记得先检查一下工作空间
getwd()
#ebi-a-GCST005810
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
exp_dat <-extract_instruments(outcomes='ieu-a-1073',p1=5e-6) #获取暴露数据（2选1）
head(exp_dat) #查看暴露数据

#读取暴露数据（2选1）
exp_dat <- read_exposure_data(
  filename = "锰-GCST90100527_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  samplesize_col = "sample_size"
)


exp_dat <- exp_dat[exp_dat$pval.exposure < 5e-6,]#筛选符合的P值

exp_dat <- clump_data(exp_dat,
                      clump_kb = 10000,
                      clump_r2 = 0.001,
                      clump_p1 = 1,
                      clump_p2 = 1,
                      pop = "EUR")#剔除连锁不平衡

write.csv(exp_dat, file="exp_dat_Cu.csv")

#提取结局数据（二选一——1）
out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = 'ebi-a-GCST005810'
)

#读取本地数据（二选一——2）
out_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "硒-GCST90100532_buildGRCh37.tsv",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value")

write.csv(out_dat, file="out_dat_Cu.csv")

#协调效应，合并数据
dat <- harmonise_data(
  exposure_dat =  exp_dat, 
  outcome_dat = out_dat
)
write.csv(dat, file="dat_Cu.csv")
write.table(dat,file="dat.txt")

dat <- read.csv("C:\\Users\\18002\\Desktop\\课题——serum trace elements in patients with osteoarthritis\\Mendelian randomization\\结果\\铜\\dat_Cu.csv")

#Steiger过滤
steiger_results <- steiger_filtering(dat)
valid_ivs <- steiger_results[steiger_results$steiger_dir == TRUE,]
write.csv(valid_ivs, file="steiger_results_Se.csv")

#根据Steiger过滤结果更新数据
dat <- dat[dat$SNP %in% valid_ivs$SNP,]
#进行MR分析
res <- mr(dat)
res
write.csv(res, file="res_Cu.csv")

or<-generate_odds_ratios(res)
write.csv(or, file="or_Cu.csv")



#查看MR的方法
mr_method_list()

#用IVW，mr-egger两种方法
res1 <- mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_ivw_fe"))
res1

#敏感性分析
#异质性检验
mr_heterogeneity(dat)
heterogeneity <-mr_heterogeneity(dat)
write.csv(heterogeneity, file="heterogeneity_Cu.csv")

#跑离群值(需要的时间很久，且数值越大时间越久)
run_mr_presso(dat,NbDistribution = 10000)

#剔除离群值
ex<-c(58,60,69)
dat1<-dat[-ex,]

##水平多效性检验
mr_pleiotropy_test(dat)
duoxx <-mr_pleiotropy_test(dat)
write.csv(duoxx, file="duoxx_Cu.csv")


##获取单独的每个SNP的效应值
res_single <- mr_singlesnp(dat)
res_single
write.csv(res_single, file="res_single_Cu.csv")

##留一法分析
res_loo <- mr_leaveoneout(dat)
res_loo 
write.csv(res_loo , file="res_loo_Cu.csv")

##散点图 
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
library(ggplot2)
ggsave(p1[[1]], file="res.pdf", width=7, height=7)


#森林图
p2 <- mr_forest_plot(res_single)
p2[[1]]
##保存图片
ggsave(p2[[1]], file="p2.pdf", width=7, height=7)


##留一法图 
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
##保存图片
ggsave(p3[[1]], file="p3.pdf", width=7, height=7)


##漏斗图 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
##保存图片
ggsave(p4[[1]], file="p4.pdf", width=7, height=7)

