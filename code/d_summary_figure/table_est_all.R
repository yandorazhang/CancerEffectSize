#-------------------------------------------------------------------
# Update Date: 07/21/2019
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


tem = matrix(0, length(inx),35); 

inx2com = c(2,4)
#-------------------------------
for(iter in 1:nrow(tr)){

  trait_name = traitlist[iter]
  trait_name_plot = traitlistplot[iter]
  output_path = paste0("../../genesis_result_new/",trait_name)
  
  ncase = tr[iter,2][[1]]; ncontrol = tr[iter,3][[1]]
  current.case = tr[iter, 5][[1]]; current.control = tr[iter,6][[1]]
  
  if(iter %in% inx2com){
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
    
    result = fit
   
    est = result$estimates$`Parameter (pic, sigmasq, a) estimates`
    sd = result$estimates$`S.D. of parameter estimates`
    llk = result$estimates$`Composite log-likelihood of fitted model`
    n_SNP_ana = result$estimates$`Total number of SNPs in the GWAS study after quality control`
    
    bic  = result$estimates$`Model selection related`$BIC
    
    nssnp = as.numeric(strsplit(result$estimates$`Number of sSNPs (sd)`, " ")[[1]][1])
    sd_nssnp = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Number of sSNPs (sd)`, " ")[[1]][2]))
    nssnp = nssnp + n_OutlierIndep
    
    nssnp1 = NA
    sd_nssnp1 = NA
    
    herit1 = NA
    sd_herit1 = NA
    herit2 = NA
    sd_herit2 = NA
    
    herit = as.numeric(strsplit(result$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][1])
    sd_herit = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][2]))
    herit = herit+herit_OutlierIndep
    
    pic = est[1]; sd_pic = sd[1];
    p1 = NA; sd_p1 = NA
    sig1 = est[2]; sd_sig1 = sd[2]
    sig2 = NA; sd_sig2 = NA
    a = est[3]; sd_a = sd[3]
    
    effect_perSNP = sig1
    sd_effect_perSNP = sd_sig1
    
    auc<-pnorm(sqrt(herit/2))
    sd_auc <- sd_herit * dnorm(sqrt(herit/2))/(sqrt(2)*2*sqrt(herit))
    
    #-----
    # number of effective associated SNPs
    ebeta2 = sig1
    ebeta4 = 3*(sig1^2)
    kap = ebeta4/(ebeta2)^2
    effect.causal = 3*(M*pic+n_OutlierIndep)/kap
  }
  
  #-----------
  # 3com GENESIS
  
  if(!iter %in% inx2com){
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
    
    result = fit
    est = result$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`
    sd = result$estimates$`S.D. of parameter estimates`
    llk = result$estimates$`Composite log-likelihood of fitted model`
    n_SNP_ana = result$estimates$`Total number of SNPs in the GWAS study after quality control`
    
    bic  = result$estimates$`Model selection related`$BIC
    
    nssnp = as.numeric(strsplit(result$estimates$`Number of sSNPs (sd)`, " ")[[1]][1])
    sd_nssnp = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Number of sSNPs (sd)`, " ")[[1]][2]))
    nssnp = nssnp + n_OutlierIndep
    
    nssnp1 = as.numeric(strsplit(result$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`, " ")[[1]][1])
    sd_nssnp1 = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`, " ")[[1]][2]))
    nssnp1 = nssnp1 + n_OutlierIndep
    
    herit1 = as.numeric(strsplit(result$estimates$`Heritability explained by the cluster with larger variance component (sd)`, " ")[[1]][1])
    sd_herit1 = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Heritability explained by the cluster with larger variance component (sd)`, " ")[[1]][2]))
    herit1 = herit1 + herit_OutlierIndep
    
    herit2 = as.numeric(strsplit(result$estimates$`Heritability explained by the cluster with samller variance component`, " ")[[1]][1])
    sd_herit2 = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Heritability explained by the cluster with samller variance component`, " ")[[1]][2]))
    
    herit = as.numeric(strsplit(result$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][1])
    sd_herit = as.numeric(gsub(".*\\((.*)\\).*", "\\1",strsplit(result$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][2]))
    herit = herit+herit_OutlierIndep
    
   
    pic = est[1]; sd_pic = sd[1];
    p1 = est[2]; sd_p1 = sd[2]
    sig1 = est[3]; sd_sig1 = sd[3]
    sig2 = est[4]; sd_sig2 = sd[4]
    a = est[5]; sd_a = sd[5]
    effect_perSNP = (p1*sig1+(1-p1)*sig2)
    
    var_est = result$estimates$`Covariance matrix of parameter estimates`
    temtem = matrix(c(0, 
                      (est[3] - est[4]), 
                      (est[2]),
                      (1-est[2]),
                      0), ncol=1)
    sd_effect_perSNP = sqrt( t(temtem) %*% var_est %*% temtem) # standard error of effect per sSNP
    
    auc<-pnorm(sqrt(herit/2))
    sd_auc <- sd_herit * dnorm(sqrt(herit/2))/(sqrt(2)*2*sqrt(herit))
    
    #-----
    # number of effective associated SNPs
    p2 = 1-p1
    ebeta4 = 3*(p1*sig1^2+p2*sig2^2)
    ebeta2 = p1*sig1 + p2*sig2
    kap = ebeta4/(ebeta2)^2
    effect.causal = 3*(M*pic+n_OutlierIndep)/kap
  }
  

  
  tem[iter,] = c(trait_name_plot, n_SNP_ana, 
                 nssnp,sd_nssnp, nssnp1,sd_nssnp1, herit1,sd_herit1,herit2,sd_herit2,
                 herit,sd_herit, bic, auc, sd_auc,
                 pic, sd_pic,
                 p1,sd_p1,sig1,sd_sig1,sig2,sd_sig2,a,sd_a,llk,
                 effect_perSNP, sd_effect_perSNP, 
                 ncase, ncontrol, current.case, current.control,
                 herit_OutlierIndep,n_OutlierIndep, effect.causal)
}

tem = data.frame(tem)
colnames(tem) = c("Trait", 
                  "Total number of SNPs in the GWAS study after quality control", 
                  "Estimated # of susceptibility SNPs (contain outlier)", "sd_num_sSNP",
                  "Estimated # of susceptibility SNPs  in cluster with larger effects  (contain outlier)","sd_num_sSNP1",
                  "Heritability explained in cluster with larger effects  (contain outlier)", "sd_herit1",
                  "Heritability explained in cluster with smaller effects","sd_herit2",
                  "Estimate of total observed scale heritability  (contain outlier)", "sd_herit",
                  "BIC","AUC  (contain outlier)", "sd_AUC",
                  "pic", "sd_pic", 
                  "p1", "sd_p1", 
                  "sig1", "sd_sig1", "sig2", "sd_sig2", 
                  "a", "sd_a", "llk",
                  "Average heritability explained per sSNP (no outlier)", "sd_effect per sSNP", 
                  "# of cases", "#of controls","current cases", "current controls",
                  "herit_OutlierIndep", "n_OutlierIndep", "# of effective causal SNPs")

output_path = paste0("../../result_png_csv_new_clump//")
write.csv(tem,file= paste0(output_path,"/table_est_all.csv"),row.names=F)



