#-------------------------------------------------------------------
# Update Date: 01/03/2020
# Create Date: 04/29/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump//")
setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")

library(data.table)
tr = fread("../../data_samplesize//cancer_sample_size.csv")

inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)



M = 1070777
source("function_figure_qqplot1.R")
# -------------------------------------------------------------------
# -------------------------------------------------------------------
png(file = paste0("../../result_png_csv_new_clump//",'SFig1.png'), width = 3500, height =3500, res =300,type="cairo")
par(mfrow = c(4,4),mar=c(5,5,2.5,2.5))
for(iter in inx[1:4]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq2com(trait_name,trait_name_plot)
  if(iter==inx[1]) mtext("2-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}

for(iter in inx[1:4]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq3com(trait_name,trait_name_plot)
  if(iter==inx[1]) mtext("3-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}


#---------
for(iter in inx[5:8]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq2com(trait_name,trait_name_plot)
  if(iter==inx[5])mtext("2-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}

for(iter in inx[5:8]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq3com(trait_name,trait_name_plot)
  if(iter==inx[5]) mtext("3-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}
dev.off()

# -------------------------------------------------------------------
# -------------------------------------------------------------------
png(file = paste0("../../result_png_csv_new_clump//",'SFig2.png'), width = 3500, height =3500, res =300,type="cairo")
par(mfrow = c(4,4),mar=c(5,5,2.5,2.5))
for(iter in inx[9:12]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq2com(trait_name,trait_name_plot)
  if(iter==inx[9])mtext("2-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}

for(iter in inx[9:12]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq3com(trait_name,trait_name_plot)
  if(iter==inx[9])mtext("3-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}


#---------
for(iter in inx[13:14]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq2com(trait_name,trait_name_plot)
  if(iter==inx[13])mtext("2-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}
plot.new();plot.new(); 

for(iter in inx[13:14]){
  traitlist = unlist(tr[iter,1])
  traitlistplot = unlist(tr[iter,7])
  trait_name = traitlist; trait_name_plot = traitlistplot; figure_qq3com(trait_name,trait_name_plot)
  if(iter==inx[13])mtext("3-component", side=2, line=2.5, cex = par("cex.lab"),col="black",font=3)
}
plot.new();plot.new(); 
dev.off()
