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
output = matrix(0,nrow(tr)*4, 7)

iter = 1
inx2com = c(2,4)

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
  
  result = fit


  current.case = tr[i,5]
  current.control = tr[i,6]
  
  
  te0=projection(est,n= as.numeric(f(1*current.case, 1*current.control)),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te1=projection(est,n= as.numeric(f(2*current.case, 2*current.control)),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te2=projection(est,n= as.numeric(f(4*current.case, 4*current.control)),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  te3=projection(est,n= as.numeric(f(1e9,1e9)),M=M,herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata)
  
  herit0 = te0$heritability
  
  output[iter, ] = unlist(c(trait_name_plot,1*current.case, 1*current.control, 
                            te0$Numdiscoveries, sqrt(exp(te0$GVpercentage*herit0/100)), 
                            te0$pheno.variance, 
                            te0$GVpercentage
  ))
  iter = iter+1
  
  output[iter, ] = unlist(c(trait_name_plot, 2*current.case, 2*current.control, 
                            te1$Numdiscoveries, sqrt(exp(te1$GVpercentage*herit0/100)), 
                            te1$pheno.variance, 
                            te1$GVpercentage))
  iter = iter+1
  
  output[iter, ] = unlist(c(trait_name_plot, 4*current.case, 4*current.control, 
                            te2$Numdiscoveries, sqrt(exp(te2$GVpercentage*herit0/100)), 
                            te2$pheno.variance, 
                            te2$GVpercentage))
  iter = iter+1
  
  output[iter, ] = unlist(c(trait_name_plot, 1e9, 1e9,
                            te3$Numdiscoveries, sqrt(exp(te3$GVpercentage*herit0/100)), 
                            te3$pheno.variance, 
                            te3$GVpercentage))
  
  iter = iter+1
 
  
}
output = data.frame(output)
colnames(output) = c("cancer","# of cases","# of controls",
                     "Predicted GW discoveries (MAF>5)", 
                     "Sibling FRR",
                     "Chip heritability explained  (log-odds-ratio/frailty scale)",
                     "Chip heritability explained as a % of total GWAS heritability"
                     )

output_path = paste0("../../result_png_csv_new_clump//")
write.csv(output, file=paste0(output_path,"/table_projection_all.csv"),row.names=F)



