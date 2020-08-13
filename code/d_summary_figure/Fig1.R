#-------------------------------------------------------------------
# Update Date: 01/29/2020
# Create Date: 04/29/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")
library(data.table)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(grid)
tr = fread("../../data_samplesize/cancer_sample_size.csv")

inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)


M = 1070777

traitlist = unlist(tr[inx,1])
traitlistplot = unlist(tr[inx,7])


mixpdf <- function(x,est){
  
  if(length(est)==5){
    pic = est[1]
    p0 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return((p0 * dnorm(x/s1)/s1 + (1-p0)*dnorm(x/s2) /s2))}
  }
  
  if(length(est)==3){
    pic = est[1]
    s1 = sqrt(est[2])
    den <- function(x){return(dnorm(x/s1)/s1)}
  }
  return(den(x))
}


inx2com = c(2,4)

x_seq0 = seq(-0.05,0.05,length.out = 1000); 
#-------------------------------#-------------------------------#-------------------------------
result = NULL
for(iter in 1:length(inx)){
  trait_name = traitlist[iter]; trait_name_plot = traitlistplot[iter]
  output_path = paste0("../../genesis_result_new/",trait_name)
  
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
    est = fit$estimates$`Parameter (pic, sigmasq, a) estimates`
  }
  
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
    
    est = fit$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`; 
  }
  y_seq = apply(matrix(x_seq0,ncol=1),1,function(t) mixpdf(t,est))
  tem = data.frame(cbind(trait_name_plot, y_seq))
  result = rbind(result, tem)
}

result$y_seq = as.numeric(as.character(result$y_seq))


x_seq = exp(x_seq0)
df1 = result[1:(length(inx1)*1000),]
df2 = result[(length(inx1)*1000+1): ((length(c(inx1,inx2))*1000)),]
df3 = result[(length(c(inx1,inx2))*1000+1): nrow(result),]

df1 = cbind(df1, rep(x_seq,length(inx1))); colnames(df1) = c("Cancer","y","x")
df2 = cbind(df2, rep(x_seq,length(inx2))); colnames(df2) = c("Cancer","y","x")
df3 = cbind(df3, rep(x_seq,length(inx3))); colnames(df3) = c("Cancer","y","x")



p1 = ggplot(data=df1) +
  scale_color_brewer(palette="Set2") +
  geom_line(size=3,aes(x=x,y=y,color=Cancer,group=Cancer)) +
  theme_bw() + labs(x="Odds ratio",y="Probability density") +
  # scale_y_continuous(limits=c(0,55)) +
  ggtitle("Cancer sites with <10,000 cases")+
  theme(text = element_text(size=35), 
        axis.text = element_text(size=35), 
        legend.position=c(0.22,0.85),legend.text=element_text(size=28),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title=element_text(size=35,face="italic",hjust=0.5),
        legend.key.size = unit(1.2, "cm"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(legend.title=element_blank())



p2 = ggplot(data=df2) +
  scale_color_brewer(palette="Paired") +
  geom_line(size=3,aes(x=x,y=y,color=Cancer,group=Cancer)) +
  theme_bw() + labs(x="Odds ratio",y="Probability density") +
  # scale_y_continuous(limits=c(0,55)) +
  ggtitle("Cancer sites with 10,000-25,000 cases")+
  theme(text = element_text(size=35), 
        axis.text = element_text(size=35), 
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position=c(0.22,0.82),legend.text=element_text(size=28),
        plot.title=element_text(size=35,face="italic",hjust=0.5),
        legend.key.size = unit(1.2, "cm"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(legend.title=element_blank())


p3 = ggplot(data=df3) +
  scale_color_brewer(palette="Set2") +
  geom_line(size=3,aes(x=x,y=y,color=Cancer,group=Cancer)) +
  theme_bw() + labs(x="Odds ratio",y="Probability density") +
  # scale_y_continuous(limits=c(0,55)) +
  ggtitle("Cancer sites with >25,000 cases")+
  theme(text = element_text(size=35), 
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=35), 
        plot.title=element_text(size=35,face="italic",hjust=0.5),
        legend.position=c(0.22,0.82),legend.text=element_text(size=28),
        legend.key.size = unit(1.2, "cm"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(legend.title=element_blank())

pp = list(length=nrow(tr))
pp[[1]] = p1; pp[[2]] = p2; pp[[3]] = p3;

pdf(paste0("../../result_png_csv_new_clump//Fig1.pdf"), width=36, height=12)
do.call(grid.arrange,c(pp,ncol=3,as.table=T))
dev.off()

setEPS()
postscript(paste0("../../result_png_csv_new_clump//Fig1.eps"), width=36, height=12)
do.call(grid.arrange,c(pp,ncol=3,as.table=T))
dev.off()

