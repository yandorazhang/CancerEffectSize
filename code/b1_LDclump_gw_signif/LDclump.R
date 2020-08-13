#-------------------------------------------------------------------
# Update Date: 04/30/2019
# Create Date: 04/30/2019
# Author: Yan (Dora) Zhang
# Email: yandorazhang@gmail.com
# cases and controls come from email
#-------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)

temp <- commandArgs(TRUE)
iter <- as.numeric(temp[1])


pthr = 5e-8
r2thr = 0.1
kbpthr = 1000
plinkpath = paste0("/dcl01/chatterj/data/yzhang/software/plink-1.07-x86_64/")
LDrefpanel_plink_path = paste0("/dcl01/chatterj/data/yzhang//Data_ldscore/1000G_EUR_Phase3_hm3only1.2/1000G_EUR_Phase3_hm3only1.2.")


traitlist = c("bcac_onco_icogs",
              "Prostate_web", "Prostate_web_Onco", 
              "Lymphoma_CLL4",
              "Colorectal","Endometrial", "Glioma",
              "HNC", "Lung", "Melanoma","Ovarian", "Esophageal", "Pancreas", "Renal","TECAC")
              #----
              # "bcac_icogs2", "Bladder",
              # "bcac_onco_icogs_erpos","bcac_onco_icogs_erneg",
              # "bcac_onco_icogs_gwas", "bcac_onco2",
              # "bcac_gwas_all",
              # "bcac_onco_icogs_gwas_erpos", "bcac_onco2_erpos", "bcac_icogs2_erpos", "bcac_gwas_erpos",
              # "bcac_onco_icogs_gwas_erneg", "bcac_onco2_erneg", "bcac_icogs2_erneg", "bcac_gwas_erneg",
              # "Lymphoma_DLBCL4", "Lymphoma_FL4","Glioma_GBM","Glioma_nonGBM",
              # "LungGadenocarcinoma","LungsquamousCell",
              # "Ovarian_serous", "Ovarian_serous_hg", "Ovarian_ser_lg_lmp",
              # "Ovarian_serous_lmp", "Ovarian_serouslowgrade", "Ovarian_mucinous_all",
              # "Ovarian_mucinous", "Ovarian_mucinous_lmp", "Ovarian_lmp",
              # "Ovarian_clearcell", "Ovarian_endometrioid")



trait_name = traitlist[iter]
load(paste0("/dcl01/chatterj/data/yzhang/2018_02_cancer/data_new/", trait_name, "/rawdata.RData")); 
# dim(rawdata) # check

outpath = paste0("/dcl01/chatterj/data/yzhang/2018_02_cancer/code_new/LDclump_gw_signif/",trait_name,"/")
if(!dir.exists(outpath)){dir.create(outpath)}


#@@@@@@@@@@@@@@@@@@@------------------------------------------------------
# freeze below part
#@@@@@@@@@@@@@@@@@@@------------------------------------------------------
#---------------------------------------#---------------------------------------
# 0. read the summaryGWAS data, filtered it to 1.2 million SNP list, 
# should remove SNPs from the WHC region. 
#---------------------------------------#---------------------------------------
outdata = rawdata
colnames(outdata)[which(colnames(outdata)=="CHR")] = "chr"
colnames(outdata)[which(colnames(outdata)=="BP")] = "bp"

#---------------------------------------#---------------------------------------
# 1. format to LDclump, only need p-value
#---------------------------------------#---------------------------------------
qassocdata = data.frame(CHR = outdata$chr, 
                      SNP = outdata$snp, BP =outdata$bp, 
                      NMISS = NA, BETA = NA, SE = NA, R2 = NA, T = NA,P=outdata$pvalue)


for(chrnum in 1:22){
  #---------------------------------------#---------------------------------------
  clumpFILE_path = paste0(outpath,"/chr",chrnum,"_r2thr",r2thr, "_pthr",pthr, "_window", kbpthr/1000,"MB")
  # LD-clumping 
  tempdata_chr = qassocdata %>% filter((CHR==chrnum))
  if (nrow(tempdata_chr)>0){
    # 1. LD clumping 
    write.table(tempdata_chr, paste0(clumpFILE_path, ".qassoc"), row.names = FALSE, quote = FALSE, sep = "\t")
    
    plinkcode = paste(paste0(plinkpath,"/plink"),
                      "--bfile", paste0(LDrefpanel_plink_path,chrnum),
                      "--clump", paste0(clumpFILE_path, ".qassoc"),
                      "--clump-p1", pthr,
                      # "--clump-p2", 1.1,
                      "--clump-r2", r2thr, 
                      "--clump-kb", kbpthr,
                      "--noweb", 
                      "--out", clumpFILE_path)
    system(plinkcode)

    
    # linuxcode = paste0("rm ", clumpFILE_path, ".clumped"); system(linuxcode)
    linuxcode = paste0("rm ", clumpFILE_path, ".log"); system(linuxcode)
    linuxcode = paste0("rm ", clumpFILE_path, ".qassoc"); system(linuxcode)
    linuxcode = paste0("rm ", clumpFILE_path, ".nosex"); system(linuxcode)
  }
}
 

snpindep = NULL;
for(chrnum in 1:22){
  clumpFILE_path = paste0(outpath,"/chr",chrnum,"_r2thr",r2thr, "_pthr",pthr, "_window", kbpthr/1000,"MB")
  tem = fread(paste0(clumpFILE_path, ".clumped"))
  if(nrow(tem)>0){
    snpindep = c(snpindep, tem$SNP)
  }
}


if(length(snpindep)>0){
  inxindep = NULL;
  for(si in snpindep){
    inxindep = c(inxindep, which(rawdata$snp == si))
  }
  
  rawdata_OutlierIndep = rawdata[inxindep,]
  n_OutlierIndep = nrow(rawdata_OutlierIndep) 
  
  beta = rawdata_OutlierIndep$z/sqrt(rawdata_OutlierIndep$n)
  tau = 1/sqrt(rawdata_OutlierIndep$n)
  herit_OutlierIndep = sum( beta^2 - tau^2 )
  
  n_signif = n_OutlierIndep
  herit_signif = herit_OutlierIndep
  save(n_signif,herit_signif,file=paste0("/dcl01/chatterj/data/yzhang/2018_02_cancer/data_new/",trait_name,"/genomewide_signif.RData") )
}

