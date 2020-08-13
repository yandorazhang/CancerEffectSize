#-------------------------------------------------------------------
# Update Date: 02/01/2020
# Create Date: 05/09/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")

rm(list=ls())
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(grid)
library(doParallel)
library(foreach)

tr0 = fread("../../data_samplesize//cancer_sample_size.csv")
inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)
tr = tr0[inx,]


lp = list(length=nrow(tr))
lauc = list(length=nrow(tr))
lgv = list(length=nrow(tr))
lr99 = list(length=nrow(tr))

for(i in 1:14){
  
  trait_name = tr[i,1]
  trait_name_plot = tr[i,7]
  load(paste0("../../result_png_csv_new_clump/RData/","Fig2_SFig34",trait_name,".RData"))
  lp[[i]] = p_p
  lauc[[i]] = p_auc
  lgv[[i]] = p_gv
  lr99[[i]] = p_r99
}



#---------------
pdf(paste0("../../result_png_csv_new_clump//Fig2.pdf"), width=28, height=28)
# do.call(grid.arrange,c(lgv,ncol=4,as.table=T))
grid.arrange(
  arrangeGrob(grobs=lapply(lgv, function(p) p ), ncol=4,
              bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=-2,gp=gpar(fontsize=25)),
              left=textGrob("y: Percentage of genetic variance explained", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90)))

dev.off()



#---------------
pdf(paste0("../../result_png_csv_new_clump//SFig3.pdf"), width=26, height=26)
# do.call(grid.arrange,c(lauc,ncol=4,as.table=T))
grid.arrange(
  arrangeGrob(grobs=lapply(lauc, function(p) p ), ncol=4,
              bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=-2,gp=gpar(fontsize=25)),
              left=textGrob("y: AUC associated with the PRS", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90)))

dev.off()


#---------------
pdf(paste0("../../result_png_csv_new_clump//SFig4.pdf"), width=28, height=28)
# do.call(grid.arrange,c(lr99,ncol=4,as.table=T))
grid.arrange(
  arrangeGrob(grobs=lapply(lr99,  function(p) p ), ncol=4, 
              bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=-2,gp=gpar(fontsize=25)), 
              left=textGrob("y: Relative risk for people at 99th centile compared to average risk of the population", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90)))
dev.off()


