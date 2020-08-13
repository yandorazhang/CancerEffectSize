#-------------------------------------------------------------------
# Update Date: 04/30/2019
# Create Date: 04/29/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
rm(list=ls())
# setwd("~/Dropbox/2018_02_cancer/code_new/results_summary_clump/")
# tr = fread("~/Dropbox/2018_02_cancer/data_samplesize/cancer_sample_size.csv")

setwd("~/OneDrive - The University Of Hong Kong/AAA/2018_02_cancer/CancerEffectSize/code/d_summary_figure")
library(data.table)
tr = fread("../../data_samplesize//cancer_sample_size.csv")

inx1 = c(23,15,43,19,39)
inx2 = c(42,16,26,13,14,27)
inx3 = c(20,40,4)

inx = c(inx1,inx2,inx3)


M = 1070777

f = function(x,y){x*y/(x+y)}

traitlist = unlist(tr[,1])
traitlistplot = unlist(tr[,7])

tem = matrix(0, length(inx),6); 
i=1
#-------------------------------
for(iter in inx){

  trait_name = traitlist[iter]
  trait_name_plot = traitlistplot[iter]

  ncase = tr[iter,2][[1]]; ncontrol = tr[iter,3][[1]]
  current.case = tr[iter, 5][[1]]; current.control = tr[iter,6][[1]]
  
  tem[i,] = c(trait_name_plot,ncase, ncontrol,
              f(as.numeric(ncase),as.numeric(ncontrol)),current.case, current.control)
  i = i+1
}

tem = data.frame(tem)
colnames(tem) = c("Trait", 
                  "# of cases", "# of controls",
                  "Effective sample size",
                  "# of cases in current GWAS", 
                  "# of controls in current GWAS")
output_path = paste0("../../result_png_csv_new_clump/")
write.csv(tem,file= paste0(output_path,"/table_samplesize.csv"),row.names=F)



