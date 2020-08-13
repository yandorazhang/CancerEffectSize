#-------------------------------------------------------------------
# Update Date: 04/30/2019
# Create Date: 04/22/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump/")
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

inx2com = c(2,4)

f = function(x,y){x*y/(x+y)}
output = matrix(0,nrow(tr)*4, 7)

output = NULL
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
  
  n.case = as.numeric(tr[i,2][[1]])
  n.control = as.numeric(tr[i,3][[1]])
  
  if(i!=4){
    load(paste0("../../data_new/",trait_name,"/genomewide_signif.RData"))
    te0 = projection(est, v=v,n= as.numeric(f(n.case,n.control)),M=M,
                     herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata,
                     CI=TRUE)
    tem = unlist(c(trait_name_plot, n_signif, herit_signif, te0$Numdiscoveries, te0$pheno.variance, te0$heritability))
    
  }
  
  
  if(i==4){
    n_signif=0; herit_signif=0
    te0 = projection(est, v=v,n= as.numeric(f(n.case,n.control)),M=M,
                     herit_OutlierIndep = herit_OutlierIndep,summarydata_OutlierIndep = summarydata,CI=F)
    tem = unlist(c(trait_name_plot, n_signif, herit_signif, te0$Numdiscoveries,NA,NA,
                   te0$pheno.variance, NA,NA,te0$heritability))
    }
  

  output = rbind(output,tem)
}

output = data.frame(output)
colnames(output) = c("cancer","n_signif", "herit_signif", 
                     "# of discoveries","# of discoveries (lowCI)", "# of discoveries (upCI)", 
                     "heritability explained", "heritability explained (lowCI)",
                     "heritability explained (upCI)", "total heritability")


output_path = paste0("../../result_png_csv_new_clump//")
write.csv(output, file=paste0(output_path,"/table_projection_identified.csv"),row.names=F)

