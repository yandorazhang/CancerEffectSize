#-------------------------------------------------------------------
# Update Date: 05/11/2019
# Create Date: 04/30/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")

library(data.table)
# tr0 = fread("~/Dropbox/2018_02_cancer/data_samplesize/cancer_sample_size.csv")
tr0 = fread("../../data_samplesize//cancer_sample_size.csv")


inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)
inx = c(inx1,inx2,inx3)

tr = tr0[inx,]

M = 1070777


traitlist = unlist(tr[,1])
traitlistplot = unlist(tr[,7])

tem = matrix(0, length(inx),4); 
#-------------------------------
for(iter in 1:nrow(tr)){

  trait_name = traitlist[iter]
  trait_name_plot = traitlistplot[iter]
  output_path = paste0("../../genesis_result_new/",trait_name)
  
  if(!file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
    load(paste0(output_path,"/fit3.RData"))
    print(c(trait_name,"no"))
    herit_OutlierIndep = 0;
    n_OutlierIndep = 0; 
  }
  
  if(file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
    load(paste0(output_path,"/fit3_RemoveOutlier.RData"))
    load(paste0("../../data_new/",trait_name,"/herit_OutlierIndep.RData"))
  }
  
  bic3  = fit$estimates$`Model selection related`$BIC
  
  
  #--------------------
  if(!file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
    load(paste0(output_path,"/bestfit2.RData"))
    print(c(trait_name,"no"))
    herit_OutlierIndep = 0;
    n_OutlierIndep = 0; 
  }
  
  if(file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
    load(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))
    load(paste0("../../data_new/",trait_name,"/herit_OutlierIndep.RData"))
  }
  bic2 = fit$estimates$`Model selection related`$BIC
  
  tem[iter,] = c(trait_name_plot,bic3,bic2,bic3<=bic2)
}

tem = data.frame(tem)
colnames(tem) = c("Trait", "BIC-3com", "BIC-2com", "BIC3<BIC2")

output_path = paste0("../../result_png_csv_new_clump//")
write.csv(tem,file= paste0(output_path,"/table_est_bic.csv"),row.names=F)


which(tem[,4]==F) # c(2,8)
