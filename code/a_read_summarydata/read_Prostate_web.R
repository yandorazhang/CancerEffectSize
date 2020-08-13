#-------------------------------------------------------------------
# Update Date: 03/16/2019
# Create Date: 03/16/2019
# Author: Yan (Dora) Zhang
# Email: yandorazhang@gmail.com
# data downloaded from http://practical.icr.ac.uk/blog/?page_id=8164
# sample size is given in GWAS_MetaAnalysis_results_summary_public.docx
#-------------------------------------------------------------------
rm(list=ls())
library(data.table)
trait_name = "Prostate_web"
n.case = 79148
n.control = 61106

#---------------------------------------
# Prostate Cancer
w_hm3.noMHC.snplist = fread("/dcl01/chatterj/data/yzhang/Data_ldscore/eur_w_ld_chr/w_hm3.noMHC.snplist")

df = NULL
for(iter in c(1:23)){
  dfchr = fread(file=paste0("/dcl01/chatterj/data/yzhang/CrossCancer/Prostate_web/meta/chromosome/meta_v3_onco_euro_overall_chr",iter,"_1_release.txt"))
  dfchr$chr = iter
  tem = merge(dfchr,w_hm3.noMHC.snplist,by="SNP",sort=F)
  df = rbind(df, tem)
  print(iter)
  gc()
}


#---------------------------------------#---------------------------------------
#1. get output file (change colnames of df) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
f = function(x,y){x*y/(x+y)}

output = data.frame(matrix(0,nrow=nrow(df), ncol=1)); colnames(output) ="chr"
output$chr = df$chr
output$bp = df$position
output$snp = df$SNP
# output$effAllele = df$Effect
# output$refAllele = df$Baseline
output$beta = df$Effect
output$se = df$StdErr
output$z = output$beta/output$se
output$n.cases=n.case
output$n.controls=n.control
output$n = f(output$n.cases, output$n.controls)
output$pvalue = 2*pnorm(-abs(output$z))
output$pvalue_raw = df$Pvalue

median(na.omit(output$z^2))/qchisq(0.5,1) #1.206964


summary(output$n.cases)
summary(output$n.controls)
# > summary(output$n.cases)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 79148   79148   79148   79148   79148   79148 
# > summary(output$n.controls)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 61106   61106   61106   61106   61106   61106 

#---------------------------------------#---------------------------------------
# 2. filter output file
# load the LD data SNPlist 
load(file="/dcl01/chatterj/data/yzhang/Data_ldscore/1000G_EUR_Phase3_MAF05_hm3only_plink/output_ld/ana_1000G_EUR_Phase3_MAF05_hm3only_cutoff01_SNPinfo.RData")

filter <- function(summarydata){
  #a. If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
  nn = as.numeric(as.character(summarydata$n))
  ikeep1 <- which(nn>=0.67*quantile(na.omit(nn), 0.9))
  summarydata <- summarydata[ikeep1,]
  
}


foutput = filter(output)
rawdata = merge(foutput,SNPinfo,by.x=c("chr","bp"), by.y=c("CHR","BP"),  sort=F)
inx = which(duplicated(rawdata$SNPname))
if(length(inx)>0) rawdata = rawdata[-inx,]


#b. Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs.
ikeep2 <- which(as.character(rawdata$SNPname) %in% w_hm3.noMHC.snplist$SNP)
rawdata <- rawdata[ikeep2,]

#---------------------------------------#---------------------------------------
# 3. Save the (rawdata.RData) for later use, and (summarydata.RData) used for GENESIS
outputdir = paste0("/dcl01/chatterj/data/yzhang/2018_02_cancer/data_new/",trait_name,"/")
if(!dir.exists(outputdir)){dir.create(outputdir)}
save(rawdata,file=paste0(outputdir,"rawdata.RData") )


summarydata <- cbind(rawdata$SNPname, rawdata$z, rawdata$n)
colnames(summarydata) <- c("snp","z","n")
summarydata <- data.frame(summarydata)
save(summarydata,file=paste0(outputdir,"summarydata.RData") )

