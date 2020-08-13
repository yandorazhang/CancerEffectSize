#-------------------------------------------------------------------
# Update Date: 07/21/2019
# Create Date: 04/22/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")

library(data.table)

M = 1070777
source("function_projection.R")


# tr0 = fread("~/Dropbox/2018_02_cancer/data_samplesize//cancer_sample_size.csv")
tr0 = fread("../../data_samplesize//cancer_sample_size.csv")


inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)

tr = tr0[inx,]
f = function(x,y){x*y/(x+y)}

inx2com = c(2,4)

output = matrix(0,length(inx), 5)
for(i in 1:nrow(tr)){

  trait_name = tr[i,1]
  trait_name_plot = tr[i,7]
  
  
  output_path = paste0("../../genesis_result_new/",trait_name)
  if(i %in% inx2com){
    if(!file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
      load(paste0(output_path,"/bestfit2.RData"))
      print(c(trait_name,"no"))
      herit_OutlierIndep = 0;
      n_OutlierIndep = 0; 
      summarydata = NULL
    }
    
    if(file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
      load(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))
      load(paste0("../../data_new/",trait_name,"/herit_OutlierIndep.RData"))
      load(paste0("../../data_new/",trait_name,"/summarydata_OutlierIndep.RData"))
    }
    est = fit$estimates$`Parameter (pic, sigmasq, a) estimates`
    v =fit$estimates$`Covariance matrix of parameter estimates`
    
  }
  
  if(!i %in% inx2com){
    if(!file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
      load(paste0(output_path,"/fit3.RData"))
      print(c(trait_name,"no"))
      herit_OutlierIndep = 0;
      n_OutlierIndep = 0; 
      summarydata = NULL
    }
    
    if(file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
      load(paste0(output_path,"/fit3_RemoveOutlier.RData"))
      load(paste0("../../data_new/",trait_name,"/herit_OutlierIndep.RData"))
      load(paste0("../../data_new/",trait_name,"/summarydata_OutlierIndep.RData"))
    }
    
    
    est = fit$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`
    v =fit$estimates$`Covariance matrix of parameter estimates`
  }
  
  te0 = projection(est, n= as.numeric(f(tr[i,5], tr[i,6])),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te1=projection(est,n= as.numeric(f(2*tr[i,5], 2*tr[i,6])),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te2=projection(est,n= as.numeric(f(4*tr[i,5], 4*tr[i,6])),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te3=projection(est,n= as.numeric(f(1e9,1e9)),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  
  herit0 = te0$heritability
  
  
  output[i, ] = unlist(c(trait_name_plot, te0$pheno.variance,te1$pheno.variance,te2$pheno.variance,te3$pheno.variance))
}
output = data.frame(output)
colnames(output) = c("cancer", "current sample size Chip heritability explained  (frailty scale)",
                     "2 folds sample size Chip heritability explained  (frailty scale)",
                     "4 folds sample size Chip heritability explained  (frailty scale)",
                     "Infinite sample size Chip heritability explained  (frailty scale)"
                     )

output_path = paste0("../../result_png_csv_new_clump//")
write.csv(output, file=paste0(output_path,"/table_projection_amber_updated.csv"),row.names = F)





