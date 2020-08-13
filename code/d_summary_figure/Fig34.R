#-------------------------------------------------------------------
# Update Date: 02/01/2020
# Create Date: 05/09/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")

library(data.table)
library(foreach)
library(ggplot2)


M = 1070777
source("function_polyriskpredict.R")
source("function_pheno_add.R")


tr0 = fread("../../data_samplesize//cancer_sample_size.csv")
inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)
tr = tr0[inx,]

tr$group = c(rep(1,5),rep(2,6),rep(3,3))
tr$trait_color =c(rep("black",length(inx)))



ff = function(x,y){x*y/(x+y)}
output = matrix(0,nrow(tr)*4, 7)


temp <- commandArgs(TRUE)
i <- as.numeric(temp[1])

for(i in 1:14){
  if(i%in%c(1,2)) {y_auclimit0 = 0.5; y_auclimit = 0.9; y_rrlimit0 = 1;  y_rrlimit=13}
  if(i==3) {y_auclimit0 = 0.5; y_auclimit = 0.9; y_rrlimit0 = 1; y_rrlimit = 13}
  if(!i%in%c(1,2,3)) {y_auclimit0 = 0.5; y_auclimit = 0.9; y_rrlimit0 = 1; y_rrlimit = 13}
  
  
  trait_name = tr[i,1]
  trait_name_plot = tr[i,7]
  trait_color = as.character(tr[i,"trait_color"])
  
  
  if(i%in%c(1,2,3,4,5)){
    xlimit_start = 10
    xlimit = 200
    n.plot_seq = seq(10e3,200e3,length.out =1000)
    n.effect_seq = n.plot_seq/4
    n.plot_seqK = n.plot_seq/1e3
  }
  
  
  if(i%in%c(6:11)){
    xlimit_start = 20
    xlimit = 400
    n.plot_seq = seq(20e3,400e3,length.out =1000)
    n.effect_seq = n.plot_seq/4
    n.plot_seqK = n.plot_seq/1e3
  }
  
  
  
  if(i%in%c(12:14)){
    xlimit_start = 40
    xlimit = 1000
    n.plot_seq = seq(40e3,1000e3,length.out =1000)
    n.effect_seq = n.plot_seq/4
    n.plot_seqK = n.plot_seq/1e3
  }
  
  output_path = paste0("../../genesis_result_new/",trait_name)
  inx2com  = c(2,4)
  
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
  
  
  num.cases = as.numeric(tr[i,5])
  num.controls = as.numeric(tr[i,6])
  
  
  #---------------------------------------------------------------------------------------#
  # 1. caclulate at current sample size---------------#
  n.effect = ff(num.cases, num.controls)
  n.plot = 4*n.effect
  
  obj = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="optimum")
  p.tem =  as.numeric(obj$alpha)
  pheno.tem1 = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.tem)
  p.best = p.tem
  pheno.best = 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1

  
  auc.best = pnorm(sqrt(pheno.best/2))
  r99.best = exp(-pheno.best/2 + qnorm(.99)*pheno.best^.5)
  
  #-----at genomewide significance levle---------------#
  objGWAS = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="GWAS")
  p.temGWAS =  as.numeric(objGWAS$alpha)
  pheno.tem1GWAS = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.temGWAS)
  p.bestGWAS = p.temGWAS
  pheno.bestGWAS = 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS

  auc.bestGWAS = pnorm(sqrt(pheno.bestGWAS/2))
  r99.bestGWAS = exp(-pheno.bestGWAS/2 + qnorm(.99)*pheno.bestGWAS^.5)
  
  stats1=as.data.frame(cbind(n.plot/1e3, p.best, pheno.best, auc.best, r99.best))
  stats1.GWAS=as.data.frame(cbind(n.plot/1e3,p.bestGWAS, pheno.bestGWAS, auc.bestGWAS, r99.bestGWAS))
  
  #---------------------------------------------------------------------------------------#
  #---------------------------------------------------------------------------------------#
  # 2 caclulate at 2*current sample size---------------#
  n.effect = 2*ff(num.cases, num.controls)
  n.plot = 4*n.effect
  
  obj = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="optimum")
  p.tem =  as.numeric(obj$alpha)
  pheno.tem1 = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.tem)
  p.best = p.tem
  pheno.best = 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  
  auc.best = pnorm(sqrt(pheno.best/2))
  r99.best = exp(-pheno.best/2 + qnorm(.99)*pheno.best^.5)
  
  #-----at genomewide significance levle---------------#
  objGWAS = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="GWAS")
  p.temGWAS =  as.numeric(objGWAS$alpha)
  pheno.tem1GWAS = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.temGWAS)
  p.bestGWAS = p.temGWAS
  pheno.bestGWAS = 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  
  auc.bestGWAS = pnorm(sqrt(pheno.bestGWAS/2))
  r99.bestGWAS = exp(-pheno.bestGWAS/2 + qnorm(.99)*pheno.bestGWAS^.5)
  
  
  
  stats2=as.data.frame(cbind(n.plot/1e3, p.best, pheno.best, auc.best, r99.best))
  stats2.GWAS=as.data.frame(cbind(n.plot/1e3, p.bestGWAS, pheno.bestGWAS, auc.bestGWAS, r99.bestGWAS))
  
  #---------------------------------------------------------------------------------------#
  #---------------------------------------------------------------------------------------#
  # 3. caclulate at 4*current sample size---------------#
  n.effect = 4*ff(num.cases, num.controls)
  n.plot = 4*n.effect
  
  obj = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="optimum")
  p.tem =  as.numeric(obj$alpha)
  pheno.tem1 = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.tem)
  p.best = p.tem
  pheno.best = 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  
  auc.best = pnorm(sqrt(pheno.best/2))
  r99.best = exp(-pheno.best/2 + qnorm(.99)*pheno.best^.5)
  
  #-----at genomewide significance levle---------------#
  objGWAS = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="GWAS")
  p.temGWAS =  as.numeric(objGWAS$alpha)
  pheno.tem1GWAS = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.temGWAS)
  p.bestGWAS = p.temGWAS
  pheno.bestGWAS = 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  
  auc.bestGWAS = pnorm(sqrt(pheno.bestGWAS/2))
  r99.bestGWAS = exp(-pheno.bestGWAS/2 + qnorm(.99)*pheno.bestGWAS^.5)
  
  
  
  stats4=as.data.frame(cbind(n.plot/1e3, p.best, pheno.best, auc.best, r99.best))
  stats4.GWAS=as.data.frame(cbind(n.plot/1e3, p.bestGWAS, pheno.bestGWAS, auc.bestGWAS, r99.bestGWAS))
  
  #---------------------------------------------------------------------------------------#
  #---------------------------------------------------------------------------------------#
  # 4. caclulate atinfitie sample size---------------#
  n.effect = 1e9
  n.plot = 4*n.effect
  
  obj = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="optimum")
  p.tem =  as.numeric(obj$alpha)
  pheno.tem1 = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.tem)
  p.best = p.tem
  pheno.best = 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  
  auc.best = pnorm(sqrt(pheno.best/2))
  r99.best = exp(-pheno.best/2 + qnorm(.99)*pheno.best^.5)
  
  #-----at genomewide significance levle---------------#
  objGWAS = polyriskpredict(N=n.effect, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="GWAS")
  p.temGWAS =  as.numeric(objGWAS$alpha)
  pheno.tem1GWAS = pheno_add(est,n.effect,summarydata_OutlierIndep = summarydata, gwas.significance=p.temGWAS)
  p.bestGWAS = p.temGWAS
  pheno.bestGWAS = 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  
  auc.bestGWAS = pnorm(sqrt(pheno.bestGWAS/2))
  r99.bestGWAS = exp(-pheno.bestGWAS/2 + qnorm(.99)*pheno.bestGWAS^.5)
  
  
  stats9=as.data.frame(cbind(n.plot/1e3, p.best, pheno.best, auc.best, r99.best))
  stats9.GWAS=as.data.frame(cbind(n.plot/1e3, p.bestGWAS, pheno.bestGWAS, auc.bestGWAS, r99.bestGWAS))
  colnames(stats9)[1] = "n.plot"
  colnames(stats9.GWAS)[1] = "n.plot"
  
  #---------------------------------------------------------------------------------------#
  df0=rbind(stats1,stats2,stats4)
  df0.GWAS = rbind(stats1.GWAS, stats2.GWAS, stats4.GWAS)
  
  df0$`Sample Size` <- c('Current',
                         'Double',
                         'Quadruple')
  
  df0$`Sample Size` <- factor(df0$`Sample Size`,levels=c('Current',
                                                         'Double',
                                                         'Quadruple'))
  colnames(df0)[1] = "n.plot"; 
  
  df0.GWAS$`Sample Size` <- c('Current',
                              'Double',
                              'Quadruple')
  
  df0.GWAS$`Sample Size` <- factor(df0.GWAS$`Sample Size`,levels=c('Current',
                                                                   'Double',
                                                                   'Quadruple'))
  colnames(df0.GWAS)[1] = "n.plot"; 
  
  
  
  # PRS is calculated with SNPs included at optimum p-value threshold
  p.best = rep(0, length(n.effect_seq)); p.bestGWAS = p.best
  pheno.best = rep(0, length(n.effect_seq)); pheno.bestGWAS = pheno.best
  
  for(k in 1:length(n.effect_seq)){
    nn = n.effect_seq[k]
    obj = polyriskpredict(N=nn, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="optimum")
    p.tem =  as.numeric(obj$alpha)
    pheno.tem1 = pheno_add(est,nn,summarydata_OutlierIndep = summarydata, gwas.significance=p.tem)
    p.best[k] = p.tem
    pheno.best[k] = 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
    
    #-----at genomewide significance levle
    objGWAS = polyriskpredict(N=nn, Ps=c(est_jiade[2],1-est_jiade[2]), Sig2s=c(est_jiade[3],est_jiade[4]), M=1070777, M1=1070777*est_jiade[1], type="GWAS")
    p.temGWAS =  as.numeric(objGWAS$alpha)
    pheno.tem1GWAS = pheno_add(est,nn,summarydata_OutlierIndep = summarydata, gwas.significance=p.temGWAS)
    p.bestGWAS[k] = p.temGWAS
    pheno.bestGWAS[k] = 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  }
  auc.best = pnorm(sqrt(pheno.best/2))
  r99.best = exp(-pheno.best/2 + qnorm(.99)*pheno.best^.5)
  df = data.frame(cbind(n.plot_seqK, p.best, pheno.best, auc.best, r99.best))
  
  #-------
  auc.bestGWAS = pnorm(sqrt(pheno.bestGWAS/2))
  r99.bestGWAS = exp(-pheno.bestGWAS/2 + qnorm(.99)*pheno.bestGWAS^.5)
  df.GWAS = data.frame(cbind(n.plot_seqK, p.bestGWAS, pheno.bestGWAS, auc.bestGWAS, r99.bestGWAS))
  
  colnames(df)[1] = "n.plot"; colnames(df.GWAS)[1] = "n.plot"
  
  df$p.best = -log10(df$p.best)
  df0$p.best = -log10(df0$p.best)
  df0.GWAS$p.bestGWAS = -log10(df0.GWAS$p.bestGWAS)
  df.GWAS$p.bestGWAS = -log10(df.GWAS$p.bestGWAS)
  
  if(i==1){
    
    p_auc = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=auc.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=auc.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=auc.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=auc.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$auc.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(y_auclimit0,y_auclimit), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position=c(0.7,0.3), #"none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    
    p_p = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=p.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=p.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=p.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=p.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_hline(yintercept = stats9$p.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(3,7),breaks=c(3,4,5,6,7))+
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position=c(0.7,0.65), #"none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    
    p_pheno = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=pheno.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=pheno.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=pheno.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=pheno.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$pheno.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(0,2), breaks = c(0, 0.5, 1, 1.5, 2)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position=c(0.7,0.3), #"none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    
    p_r99 = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=r99.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=r99.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=r99.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=r99.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$r99.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_log10(limits=c(1,y_rrlimit), breaks = c(1,2,3,4,6,8,10,12)) +
      # scale_y_continuous(limits=c(y_rrlimit0,y_rrlimit), breaks = c(0,2,4,6,8,10,12)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position=c(0.7,0.3), #"none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
  }
  
  if(i!=1){
    p_auc = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=auc.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=auc.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=auc.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=auc.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$auc.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(y_auclimit0,y_auclimit), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position="none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    p_p = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=p.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=p.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=p.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=p.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_hline(yintercept = stats9$p.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(3,7),breaks=c(3,4,5,6,7))+
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position = "none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    p_pheno = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=pheno.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=pheno.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=pheno.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=pheno.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$pheno.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_continuous(limits=c(0,2), breaks = c(0, 0.5, 1, 1.5, 2)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position="none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
    
    
    
    p_r99 = ggplot() +
      geom_line(data=df, aes(x=n.plot, y=r99.best), size=2.5,color="grey38") +
      geom_point(data=df0,aes(x=n.plot, y=r99.best,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      # geom_line(data=df.GWAS, aes(x=n.plot, y=r99.bestGWAS), size=2.5,color="grey38") +
      # geom_point(data=df0.GWAS,aes(x=n.plot, y=r99.bestGWAS,fill=`Sample Size`, color=`Sample Size`), size=6)+    
      geom_hline(yintercept = stats9$r99.best,linetype="dotted",size=2,col="red")+
      theme_bw() + scale_color_manual(values=c("royalblue1","palevioletred","lightgoldenrod3"))+
      labs(x = paste0("\n"),y = "\n") +
      ggtitle(trait_name_plot)+
      scale_x_continuous(limits=c(0,xlimit),breaks = seq(0,xlimit,length.out=5)) +
      scale_y_log10(limits=c(1,y_rrlimit), breaks = c(1,2,3,4,6,8,10,12)) +
      # scale_y_continuous(limits=c(y_rrlimit0,y_rrlimit), breaks = c(0,2,4,6,8,10,12)) +
      theme(text = element_text(size=24), 
            # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title=element_text(size=24,face="italic",hjust=0.5),
            axis.text = element_text(size=24), 
            legend.position="none",
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
      theme(legend.title=element_blank())
  }
  
  save(p_p, p_auc, p_pheno, p_r99, 
       file = paste0("../../result_png_csv_new_clump/RData/","Fig34",trait_name,".RData"))
  
}
  


