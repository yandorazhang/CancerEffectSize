#-------------------------------------------------------------------
# Update Date: 07/21/2019
# Create Date: 04/29/2019
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
tr = tr0[c(inx1,inx2,inx3),]



f = function(x,y){x*y/(x+y)}
output = matrix(0,nrow(tr)*4, 7)


ns = c(80,  #1
       180, #2
       60,  #3
       330, #4
       190, #5
       250, #6
       170, #7
       110, #8
       220, #9
       250, #10
       270, #11
       1000, #12
       370, #13
       780)*1e3 #14

inx2com = c(2,4)

for(i in 1:length(ns)){
  
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
    
    est2com = est; 
    est_jiade = c(est2com[1], 0.5, est2com[2], est2com[2], est2com[3])
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
    est_jiade = est
  }
  

  result = fit

  sibling = result$estimates$`Sibling risk (sd)`
  sibling = as.numeric(strsplit(sibling, " ")[[1]][1])
  herit0 = as.numeric(strsplit(result$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][1])
  
  print(c(i, projection(est,v=NULL, n= f(ns[i],ns[i]), herit_OutlierIndep = herit_OutlierIndep,
             summarydata_OutlierIndep = summarydata,M=M)$GVpercentage))
  
}
