#-------------------------------------------------------------------
# Update Date: 03/16/2019
# Create Date: 03/16/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
temp <- commandArgs(TRUE)
iter <- as.numeric(temp[1])
startingpic <- as.numeric(temp[2])
LDwindow <- 1
LDcutoff <- 0.1

library(data.table)
library(GENESIS)
cores = 24

traitlist = c("bcac_onco_icogs","bcac_icogs2",
              "Prostate_web", "Prostate_web_Onco", 
              "Bladder",
              "bcac_onco_icogs_erpos","bcac_onco_icogs_erneg",
              "bcac_onco_icogs_gwas", "bcac_onco2", 
              "bcac_gwas_all", 
  "bcac_onco_icogs_gwas_erpos", "bcac_onco2_erpos", "bcac_icogs2_erpos", "bcac_gwas_erpos",
  "bcac_onco_icogs_gwas_erneg", "bcac_onco2_erneg", "bcac_icogs2_erneg", "bcac_gwas_erneg",
  "Lymphoma_CLL4","Lymphoma_DLBCL4", "Lymphoma_FL4",
              "Colorectal","Endometrial", "Glioma","Glioma_GBM","Glioma_nonGBM",
              "HNC", "Lung", "LungGadenocarcinoma","LungsquamousCell",
              "Melanoma","Ovarian", "Esophageal",
              "Ovarian_serous", "Ovarian_serous_hg", "Ovarian_ser_lg_lmp",
              "Ovarian_serous_lmp", "Ovarian_serouslowgrade", "Ovarian_mucinous_all",
              "Ovarian_mucinous", "Ovarian_mucinous_lmp", "Ovarian_lmp",
              "Ovarian_clearcell", "Ovarian_endometrioid", "Pancreas", "Renal","TECAC")


trait_name = traitlist[iter]
load(paste0("../data_new/",trait_name, "/summarydata.RData"))

output_path = paste0("../genesis_result_new/",trait_name)
if(!dir.exists(output_path)){dir.create(output_path)}


fit = genesis(summarydata, filter=FALSE, 
              modelcomponents=2, cores=cores, LDcutoff=LDcutoff,LDwindow=LDwindow,M=1070777,
              c0=10, print=TRUE, printfreq=10, starting=NA,startingpic=startingpic, tolerance=NA, 
              qqplot=TRUE, qqplotCI.coverage=0.8, qqplot.name=paste0(output_path,"/fit2_pic",startingpic), 
              summaryGWASdata.save=FALSE, qqplotdata.save=T, 
              herit.liability=FALSE)

save(fit, file=paste0(output_path,"/fit2_pic",startingpic,".RData"))



#---------------------------------#---------------------------------
# get the best fit2
llk = NULL
for(startingpic in c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)){
  load(paste0(output_path,"/fit2_pic",startingpic,".RData"))
  llk = c(llk, fit$estimates$`Composite log-likelihood of fitted model`)
}
bestpic = c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)[which.max(llk)]
load(paste0(output_path,"/fit2_pic",startingpic,".RData"))
save(fit,bestpic, file=paste0(output_path,"/bestfit2.RData"))


#---------------------------------#---------------------------------
# set initial value for 3-component model 
load(paste0(output_path,"/bestfit2.RData"))  
est <- fit$estimates$`Parameter (pic, sigmasq, a) estimates`
starting <- rep(0,5)
starting[1] <- est[1]
starting[2] <- 1/9
starting[3] <- est[2]*5
starting[4] <- starting[3]/10
starting[5] <- est[3]


fit = genesis(summarydata, filter=FALSE, 
              modelcomponents=3, cores=cores, LDcutoff=LDcutoff,LDwindow=LDwindow,M=1070777,
              c0=10, print=TRUE, printfreq=10, starting=starting, tolerance=NA, 
              qqplot=TRUE, qqplotCI.coverage=0.8, qqplot.name=paste0(output_path,"/fit3"), 
              summaryGWASdata.save=FALSE, qqplotdata.save=T, 
              herit.liability=FALSE)

save(fit, file=paste0(output_path,"/fit3.RData"))

