rm(list=ls())
library(iCARE)
library(cowplot)
library(openxlsx)
library(reshape2)
library(grid)
set.seed(7918)
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections')
#setwd('/Users/amberwilcox/NCI Google Drive (amber.wilcox@nih.gov)/Crosscancer Effect Size/Risk Projections')
############################
# Input File Coding
############################
# Create model formula
model.PRS=outcome ~ PRS + var

# Create covariate information
v1 = v2 = list()

v1$names = "PRS"
v1$type = "continuous"

v2$names = "var"
v2$type = "continuous"

cov_info_PRS = list(v1,v2)

# Create beta file
beta_PRS = c(1,0)
names(beta_PRS) = c("PRS", "var")

############################
# Load Rate Files
############################
breast_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Breast (overall)')
#bc_ERpos_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Breast (ERpos)')
#bc_ERneg_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Breast (ERneg)')
colorectal_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Colorectal')
endometrial_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Endometrial')
esophageal_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Esophageal')
glioma_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Glioma')
headneck_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Head and Neck')
lung_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Lung (Overall)')
#lungAC_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Lung (Adeno)')
#lungSCC_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Lung (SCC)')
lymphCLL_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='CLL')
#lymphFL_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Leukemia FL')
melanoma_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Melanoma')
ovarian_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Ovarian')
pancreatic_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Pancreatic')
prostate_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Prostate')
renal_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Renal')
testicular_inc=read.xlsx('IncRates_All_Cancers2015.xlsx', sheet ='Testicular')

setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/CDC Wonder Mortality Files')
mort_no_breast=read.xlsx('mort_no_breast_2016.xlsx', sheet = 1)
mort_no_colorectal=read.xlsx('mort_no_colorectal_2016.xlsx', sheet = 1)
mort_no_endometrial=read.xlsx('mort_no_endometrial_2016.xlsx', sheet=1)
mort_no_esophageal=read.xlsx('mort_no_esophageal_2016.xlsx', sheet=1)
mort_no_glioma=read.xlsx('mort_no_glioma_2016.xlsx', sheet=1)
mort_no_headneck=read.xlsx('mort_no_HNC_2016.xlsx', sheet=1)
mort_no_lung=read.xlsx('mort_no_lung_2016.xlsx', sheet=1)
mort_no_lymphCLL=read.xlsx('mort_no_CLL_2016.xlsx', sheet=1)
mort_no_lymphFL=read.xlsx('mort_no_FL_2016.xlsx', sheet=1)
mort_no_melanoma=read.xlsx('mort_no_melanoma_2016.xlsx', sheet=1)
mort_no_ovarian=read.xlsx('mort_no_ovarian_2016.xlsx', sheet=1)
mort_no_pancreatic=read.xlsx('mort_no_pancreatic_2016.xlsx', sheet=1)
mort_no_prostate=read.xlsx('mort_no_prostate_2016.xlsx', sheet=1)
mort_no_renal=read.xlsx('mort_no_renal_2016.xlsx', sheet=1)
mort_no_testicular=read.xlsx('mort_no_testicular_2016.xlsx', sheet=1)

#########################################
############################
# Calculating Absolute Risk
############################
#########################################
      ############################
      # Breast cancer
      ############################
var1=0.278144652
var2=0.337157763
var3=0.417057802
var4=0.602166453
##### CURRENT SAMPLE SIZE #####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
breast_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
BCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
BCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
BCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
BCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
BCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
BCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
BCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
BCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
BCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

BCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(BCA_PRS_30to34$risk), as.numeric(BCA_PRS_35to39$risk), as.numeric(BCA_PRS_40to44$risk), as.numeric(BCA_PRS_45to49$risk),
                                           as.numeric(BCA_PRS_50to54$risk), as.numeric(BCA_PRS_55to59$risk), as.numeric(BCA_PRS_60to64$risk), as.numeric(BCA_PRS_65to69$risk),
                                           as.numeric(BCA_PRS_70to74$risk)))

colnames(BCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

BCA_PRS_Average_Data1$risk_wtd_avg=(BCA_PRS_Average_Data1$risk_30to34*0.105319686 + BCA_PRS_Average_Data1$risk_35to39*0.102730061 + BCA_PRS_Average_Data1$risk_40to44*0.095906423 + 
                                       BCA_PRS_Average_Data1$risk_45to49*0.110520685 + BCA_PRS_Average_Data1$risk_50to54*0.120983294 + BCA_PRS_Average_Data1$risk_55to59*0.133676551 +
                                       BCA_PRS_Average_Data1$risk_60to64*0.12774873 + BCA_PRS_Average_Data1$risk_65to69*0.112737815 + BCA_PRS_Average_Data1$risk_70to74*0.090376754)

BCA_PRS_AR_CurrentSS=BCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
breast_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
BCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
BCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
BCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
BCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
BCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
BCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
BCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
BCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
BCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

BCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(BCA_PRS_30to34$risk), as.numeric(BCA_PRS_35to39$risk), as.numeric(BCA_PRS_40to44$risk), as.numeric(BCA_PRS_45to49$risk),
                                           as.numeric(BCA_PRS_50to54$risk), as.numeric(BCA_PRS_55to59$risk), as.numeric(BCA_PRS_60to64$risk), as.numeric(BCA_PRS_65to69$risk),
                                           as.numeric(BCA_PRS_70to74$risk)))

colnames(BCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

BCA_PRS_Average_Data2$risk_wtd_avg=(BCA_PRS_Average_Data2$risk_30to34*0.105319686 + BCA_PRS_Average_Data2$risk_35to39*0.102730061 + BCA_PRS_Average_Data2$risk_40to44*0.095906423 + 
                                       BCA_PRS_Average_Data2$risk_45to49*0.110520685 + BCA_PRS_Average_Data2$risk_50to54*0.120983294 + BCA_PRS_Average_Data2$risk_55to59*0.133676551 +
                                       BCA_PRS_Average_Data2$risk_60to64*0.12774873 + BCA_PRS_Average_Data2$risk_65to69*0.112737815 + BCA_PRS_Average_Data2$risk_70to74*0.090376754)

BCA_PRS_AR_2xSS=BCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
breast_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
BCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
BCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
BCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
BCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
BCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
BCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
BCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
BCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
BCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

BCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(BCA_PRS_30to34$risk), as.numeric(BCA_PRS_35to39$risk), as.numeric(BCA_PRS_40to44$risk), as.numeric(BCA_PRS_45to49$risk),
                                           as.numeric(BCA_PRS_50to54$risk), as.numeric(BCA_PRS_55to59$risk), as.numeric(BCA_PRS_60to64$risk), as.numeric(BCA_PRS_65to69$risk),
                                           as.numeric(BCA_PRS_70to74$risk)))

colnames(BCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

BCA_PRS_Average_Data3$risk_wtd_avg=(BCA_PRS_Average_Data3$risk_30to34*0.105319686 + BCA_PRS_Average_Data3$risk_35to39*0.102730061 + BCA_PRS_Average_Data3$risk_40to44*0.095906423 + 
                                       BCA_PRS_Average_Data3$risk_45to49*0.110520685 + BCA_PRS_Average_Data3$risk_50to54*0.120983294 + BCA_PRS_Average_Data3$risk_55to59*0.133676551 +
                                       BCA_PRS_Average_Data3$risk_60to64*0.12774873 + BCA_PRS_Average_Data3$risk_65to69*0.112737815 + BCA_PRS_Average_Data3$risk_70to74*0.090376754)

BCA_PRS_AR_4xSS=BCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
breast_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
BCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
BCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
BCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
BCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                                      model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
BCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
BCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
BCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
BCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
BCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = breast_ref,
                           model.disease.incidence.rates = breast_inc, model.competing.incidence.rates = mort_no_breast, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = breast_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

BCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(BCA_PRS_30to34$risk), as.numeric(BCA_PRS_35to39$risk), as.numeric(BCA_PRS_40to44$risk), as.numeric(BCA_PRS_45to49$risk),
                                           as.numeric(BCA_PRS_50to54$risk), as.numeric(BCA_PRS_55to59$risk), as.numeric(BCA_PRS_60to64$risk), as.numeric(BCA_PRS_65to69$risk),
                                           as.numeric(BCA_PRS_70to74$risk)))

colnames(BCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

BCA_PRS_Average_Data4$risk_wtd_avg=(BCA_PRS_Average_Data4$risk_30to34*0.105319686 + BCA_PRS_Average_Data4$risk_35to39*0.102730061 + BCA_PRS_Average_Data4$risk_40to44*0.095906423 + 
                                       BCA_PRS_Average_Data4$risk_45to49*0.110520685 + BCA_PRS_Average_Data4$risk_50to54*0.120983294 + BCA_PRS_Average_Data4$risk_55to59*0.133676551 +
                                       BCA_PRS_Average_Data4$risk_60to64*0.12774873 + BCA_PRS_Average_Data4$risk_65to69*0.112737815 + BCA_PRS_Average_Data4$risk_70to74*0.090376754)

BCA_PRS_AR_InfSS=BCA_PRS_Average_Data4$risk_wtd_avg

save(BCA_PRS_AR_CurrentSS,BCA_PRS_AR_2xSS,BCA_PRS_AR_4xSS,BCA_PRS_AR_InfSS, file='BCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
BCA_PRS_AR_30to34=as.data.frame(cbind(BCA_PRS_Average_Data1$risk_30to34, BCA_PRS_Average_Data2$risk_30to34,
                                      BCA_PRS_Average_Data3$risk_30to34, BCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(BCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

BCA_PRS_AR_40to44=as.data.frame(cbind(BCA_PRS_Average_Data1$risk_40to44, BCA_PRS_Average_Data2$risk_40to44,
                                      BCA_PRS_Average_Data3$risk_40to44, BCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(BCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

BCA_PRS_AR_50to54=as.data.frame(cbind(BCA_PRS_Average_Data1$risk_50to54, BCA_PRS_Average_Data2$risk_50to54,
                                      BCA_PRS_Average_Data3$risk_50to54, BCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(BCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

BCA_PRS_AR_60to64=as.data.frame(cbind(BCA_PRS_Average_Data1$risk_60to64, BCA_PRS_Average_Data2$risk_60to64,
                                      BCA_PRS_Average_Data3$risk_60to64, BCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(BCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(BCA_PRS_AR_30to34, BCA_PRS_AR_40to44, BCA_PRS_AR_50to54, BCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/BCA_PRS_AR_AgeStrat.Rda')
      ############################
      # Colorectal cancer
      ############################
var1=0.034061422
var2=0.090161616
var3=0.196584146
var4=0.428932631
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
colorectal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CRCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CRCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CRCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CRCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CRCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CRCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CRCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CRCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CRCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CRCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(CRCA_PRS_30to34$risk), as.numeric(CRCA_PRS_35to39$risk), as.numeric(CRCA_PRS_40to44$risk), as.numeric(CRCA_PRS_45to49$risk),
                                           as.numeric(CRCA_PRS_50to54$risk), as.numeric(CRCA_PRS_55to59$risk), as.numeric(CRCA_PRS_60to64$risk), as.numeric(CRCA_PRS_65to69$risk),
                                           as.numeric(CRCA_PRS_70to74$risk)))

colnames(CRCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CRCA_PRS_Average_Data1$risk_wtd_avg=(CRCA_PRS_Average_Data1$risk_30to34*0.108040659 + CRCA_PRS_Average_Data1$risk_35to39*0.104987889 + CRCA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                     CRCA_PRS_Average_Data1$risk_45to49*0.112278705 + CRCA_PRS_Average_Data1$risk_50to54*0.121739409 + CRCA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                     CRCA_PRS_Average_Data1$risk_60to64*0.125715739 + CRCA_PRS_Average_Data1$risk_65to69*0.109659337 + CRCA_PRS_Average_Data1$risk_70to74*0.086579466)

CRCA_PRS_AR_CurrentSS=CRCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
colorectal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CRCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CRCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CRCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CRCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CRCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CRCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CRCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CRCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CRCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CRCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(CRCA_PRS_30to34$risk), as.numeric(CRCA_PRS_35to39$risk), as.numeric(CRCA_PRS_40to44$risk), as.numeric(CRCA_PRS_45to49$risk),
                                           as.numeric(CRCA_PRS_50to54$risk), as.numeric(CRCA_PRS_55to59$risk), as.numeric(CRCA_PRS_60to64$risk), as.numeric(CRCA_PRS_65to69$risk),
                                           as.numeric(CRCA_PRS_70to74$risk)))

colnames(CRCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CRCA_PRS_Average_Data2$risk_wtd_avg=(CRCA_PRS_Average_Data2$risk_30to34*0.108040659 + CRCA_PRS_Average_Data2$risk_35to39*0.104987889 + CRCA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                     CRCA_PRS_Average_Data2$risk_45to49*0.112278705 + CRCA_PRS_Average_Data2$risk_50to54*0.121739409 + CRCA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                     CRCA_PRS_Average_Data2$risk_60to64*0.125715739 + CRCA_PRS_Average_Data2$risk_65to69*0.109659337 + CRCA_PRS_Average_Data2$risk_70to74*0.086579466)

CRCA_PRS_AR_2xSS=CRCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
colorectal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CRCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CRCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CRCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CRCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CRCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CRCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CRCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CRCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CRCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CRCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(CRCA_PRS_30to34$risk), as.numeric(CRCA_PRS_35to39$risk), as.numeric(CRCA_PRS_40to44$risk), as.numeric(CRCA_PRS_45to49$risk),
                                           as.numeric(CRCA_PRS_50to54$risk), as.numeric(CRCA_PRS_55to59$risk), as.numeric(CRCA_PRS_60to64$risk), as.numeric(CRCA_PRS_65to69$risk),
                                           as.numeric(CRCA_PRS_70to74$risk)))

colnames(CRCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CRCA_PRS_Average_Data3$risk_wtd_avg=(CRCA_PRS_Average_Data3$risk_30to34*0.108040659 + CRCA_PRS_Average_Data3$risk_35to39*0.104987889 + CRCA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                     CRCA_PRS_Average_Data3$risk_45to49*0.112278705 + CRCA_PRS_Average_Data3$risk_50to54*0.121739409 + CRCA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                     CRCA_PRS_Average_Data3$risk_60to64*0.125715739 + CRCA_PRS_Average_Data3$risk_65to69*0.109659337 + CRCA_PRS_Average_Data3$risk_70to74*0.086579466)

CRCA_PRS_AR_4xSS=CRCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
colorectal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CRCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CRCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CRCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CRCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CRCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CRCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CRCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CRCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CRCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = colorectal_ref,
                           model.disease.incidence.rates = colorectal_inc, model.competing.incidence.rates = mort_no_colorectal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = colorectal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CRCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(CRCA_PRS_30to34$risk), as.numeric(CRCA_PRS_35to39$risk), as.numeric(CRCA_PRS_40to44$risk), as.numeric(CRCA_PRS_45to49$risk),
                                           as.numeric(CRCA_PRS_50to54$risk), as.numeric(CRCA_PRS_55to59$risk), as.numeric(CRCA_PRS_60to64$risk), as.numeric(CRCA_PRS_65to69$risk),
                                           as.numeric(CRCA_PRS_70to74$risk)))

colnames(CRCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CRCA_PRS_Average_Data4$risk_wtd_avg=(CRCA_PRS_Average_Data4$risk_30to34*0.108040659 + CRCA_PRS_Average_Data4$risk_35to39*0.104987889 + CRCA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                     CRCA_PRS_Average_Data4$risk_45to49*0.112278705 + CRCA_PRS_Average_Data4$risk_50to54*0.121739409 + CRCA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                     CRCA_PRS_Average_Data4$risk_60to64*0.125715739 + CRCA_PRS_Average_Data4$risk_65to69*0.109659337 + CRCA_PRS_Average_Data4$risk_70to74*0.086579466)

CRCA_PRS_AR_InfSS=CRCA_PRS_Average_Data4$risk_wtd_avg

save(CRCA_PRS_AR_CurrentSS,CRCA_PRS_AR_2xSS,CRCA_PRS_AR_4xSS,CRCA_PRS_AR_InfSS, 
     file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/CRCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
CRCA_PRS_AR_30to34=as.data.frame(cbind(CRCA_PRS_Average_Data1$risk_30to34, CRCA_PRS_Average_Data2$risk_30to34,
                                      CRCA_PRS_Average_Data3$risk_30to34, CRCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(CRCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CRCA_PRS_AR_40to44=as.data.frame(cbind(CRCA_PRS_Average_Data1$risk_40to44, CRCA_PRS_Average_Data2$risk_40to44,
                                      CRCA_PRS_Average_Data3$risk_40to44, CRCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(CRCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CRCA_PRS_AR_50to54=as.data.frame(cbind(CRCA_PRS_Average_Data1$risk_50to54, CRCA_PRS_Average_Data2$risk_50to54,
                                      CRCA_PRS_Average_Data3$risk_50to54, CRCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(CRCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CRCA_PRS_AR_60to64=as.data.frame(cbind(CRCA_PRS_Average_Data1$risk_60to64, CRCA_PRS_Average_Data2$risk_60to64,
                                      CRCA_PRS_Average_Data3$risk_60to64, CRCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(CRCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(CRCA_PRS_AR_30to34, CRCA_PRS_AR_40to44, CRCA_PRS_AR_50to54, CRCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/CRCA_PRS_AR_AgeStrat.Rda')
      ############################
      # Endometrial cancer
      ############################
var1=0.050261388
var2=0.084606019
var3=0.140728816
var4=0.271009918
##### CURRENT SAMPLE SIZE #####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
endometrial_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
UECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
UECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
UECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
UECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
UECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
UECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
UECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
UECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
UECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

UECA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(UECA_PRS_30to34$risk), as.numeric(UECA_PRS_35to39$risk), as.numeric(UECA_PRS_40to44$risk), as.numeric(UECA_PRS_45to49$risk),
                                           as.numeric(UECA_PRS_50to54$risk), as.numeric(UECA_PRS_55to59$risk), as.numeric(UECA_PRS_60to64$risk), as.numeric(UECA_PRS_65to69$risk),
                                           as.numeric(UECA_PRS_70to74$risk)))

colnames(UECA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

UECA_PRS_Average_Data1$risk_wtd_avg=(UECA_PRS_Average_Data1$risk_30to34*0.105319686 + UECA_PRS_Average_Data1$risk_35to39*0.102730061 + UECA_PRS_Average_Data1$risk_40to44*0.095906423 + 
                                       UECA_PRS_Average_Data1$risk_45to49*0.110520685 + UECA_PRS_Average_Data1$risk_50to54*0.120983294 + UECA_PRS_Average_Data1$risk_55to59*0.133676551 +
                                       UECA_PRS_Average_Data1$risk_60to64*0.12774873 + UECA_PRS_Average_Data1$risk_65to69*0.112737815 + UECA_PRS_Average_Data1$risk_70to74*0.090376754)

UECA_PRS_AR_CurrentSS=UECA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
endometrial_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
UECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
UECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
UECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
UECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
UECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
UECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
UECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
UECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
UECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

UECA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(UECA_PRS_30to34$risk), as.numeric(UECA_PRS_35to39$risk), as.numeric(UECA_PRS_40to44$risk), as.numeric(UECA_PRS_45to49$risk),
                                           as.numeric(UECA_PRS_50to54$risk), as.numeric(UECA_PRS_55to59$risk), as.numeric(UECA_PRS_60to64$risk), as.numeric(UECA_PRS_65to69$risk),
                                           as.numeric(UECA_PRS_70to74$risk)))

colnames(UECA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

UECA_PRS_Average_Data2$risk_wtd_avg=(UECA_PRS_Average_Data2$risk_30to34*0.105319686 + UECA_PRS_Average_Data2$risk_35to39*0.102730061 + UECA_PRS_Average_Data2$risk_40to44*0.095906423 + 
                                       UECA_PRS_Average_Data2$risk_45to49*0.110520685 + UECA_PRS_Average_Data2$risk_50to54*0.120983294 + UECA_PRS_Average_Data2$risk_55to59*0.133676551 +
                                       UECA_PRS_Average_Data2$risk_60to64*0.12774873 + UECA_PRS_Average_Data2$risk_65to69*0.112737815 + UECA_PRS_Average_Data2$risk_70to74*0.090376754)

UECA_PRS_AR_2xSS=UECA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
endometrial_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
UECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
UECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
UECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
UECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
UECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
UECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
UECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
UECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
UECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

UECA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(UECA_PRS_30to34$risk), as.numeric(UECA_PRS_35to39$risk), as.numeric(UECA_PRS_40to44$risk), as.numeric(UECA_PRS_45to49$risk),
                                           as.numeric(UECA_PRS_50to54$risk), as.numeric(UECA_PRS_55to59$risk), as.numeric(UECA_PRS_60to64$risk), as.numeric(UECA_PRS_65to69$risk),
                                           as.numeric(UECA_PRS_70to74$risk)))

colnames(UECA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

UECA_PRS_Average_Data3$risk_wtd_avg=(UECA_PRS_Average_Data3$risk_30to34*0.105319686 + UECA_PRS_Average_Data3$risk_35to39*0.102730061 + UECA_PRS_Average_Data3$risk_40to44*0.095906423 + 
                                       UECA_PRS_Average_Data3$risk_45to49*0.110520685 + UECA_PRS_Average_Data3$risk_50to54*0.120983294 + UECA_PRS_Average_Data3$risk_55to59*0.133676551 +
                                       UECA_PRS_Average_Data3$risk_60to64*0.12774873 + UECA_PRS_Average_Data3$risk_65to69*0.112737815 + UECA_PRS_Average_Data3$risk_70to74*0.090376754)

UECA_PRS_AR_4xSS=UECA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
endometrial_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
UECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
UECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
UECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
UECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                                      model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
UECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
UECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
UECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
UECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
UECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = endometrial_ref,
                           model.disease.incidence.rates = endometrial_inc, model.competing.incidence.rates = mort_no_endometrial, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = endometrial_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

UECA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(UECA_PRS_30to34$risk), as.numeric(UECA_PRS_35to39$risk), as.numeric(UECA_PRS_40to44$risk), as.numeric(UECA_PRS_45to49$risk),
                                           as.numeric(UECA_PRS_50to54$risk), as.numeric(UECA_PRS_55to59$risk), as.numeric(UECA_PRS_60to64$risk), as.numeric(UECA_PRS_65to69$risk),
                                           as.numeric(UECA_PRS_70to74$risk)))

colnames(UECA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

UECA_PRS_Average_Data4$risk_wtd_avg=(UECA_PRS_Average_Data4$risk_30to34*0.105319686 + UECA_PRS_Average_Data4$risk_35to39*0.102730061 + UECA_PRS_Average_Data4$risk_40to44*0.095906423 + 
                                       UECA_PRS_Average_Data4$risk_45to49*0.110520685 + UECA_PRS_Average_Data4$risk_50to54*0.120983294 + UECA_PRS_Average_Data4$risk_55to59*0.133676551 +
                                       UECA_PRS_Average_Data4$risk_60to64*0.12774873 + UECA_PRS_Average_Data4$risk_65to69*0.112737815 + UECA_PRS_Average_Data4$risk_70to74*0.090376754)

UECA_PRS_AR_InfSS=UECA_PRS_Average_Data4$risk_wtd_avg

save(UECA_PRS_AR_CurrentSS,UECA_PRS_AR_2xSS,UECA_PRS_AR_4xSS,UECA_PRS_AR_InfSS, file='UECA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
UECA_PRS_AR_30to34=as.data.frame(cbind(UECA_PRS_Average_Data1$risk_30to34, UECA_PRS_Average_Data2$risk_30to34,
                                      UECA_PRS_Average_Data3$risk_30to34, UECA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(UECA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

UECA_PRS_AR_40to44=as.data.frame(cbind(UECA_PRS_Average_Data1$risk_40to44, UECA_PRS_Average_Data2$risk_40to44,
                                      UECA_PRS_Average_Data3$risk_40to44, UECA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(UECA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

UECA_PRS_AR_50to54=as.data.frame(cbind(UECA_PRS_Average_Data1$risk_50to54, UECA_PRS_Average_Data2$risk_50to54,
                                      UECA_PRS_Average_Data3$risk_50to54, UECA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(UECA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

UECA_PRS_AR_60to64=as.data.frame(cbind(UECA_PRS_Average_Data1$risk_60to64, UECA_PRS_Average_Data2$risk_60to64,
                                      UECA_PRS_Average_Data3$risk_60to64, UECA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(UECA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(UECA_PRS_AR_30to34, UECA_PRS_AR_40to44, UECA_PRS_AR_50to54, UECA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/UECA_PRS_AR_AgeStrat.Rda')
      ############################
      # Esophageal cancer
      ############################
var1=0.02267308
var2=0.14376533
var3=0.440513759
var4=1.242902083
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
esophageal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
ECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
ECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
ECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
ECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
ECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
ECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
ECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
ECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
ECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

ECA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(ECA_PRS_30to34$risk), as.numeric(ECA_PRS_35to39$risk), as.numeric(ECA_PRS_40to44$risk), as.numeric(ECA_PRS_45to49$risk), 
                                          as.numeric(ECA_PRS_50to54$risk), as.numeric(ECA_PRS_55to59$risk), as.numeric(ECA_PRS_60to64$risk), as.numeric(ECA_PRS_65to69$risk),
                                          as.numeric(ECA_PRS_70to74$risk)))

colnames(ECA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

ECA_PRS_Average_Data1$risk_wtd_avg=(ECA_PRS_Average_Data1$risk_30to34*0.108040659 + ECA_PRS_Average_Data1$risk_35to39*0.104987889 + ECA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                    ECA_PRS_Average_Data1$risk_45to49*0.112278705 + ECA_PRS_Average_Data1$risk_50to54*0.121739409 + ECA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                    ECA_PRS_Average_Data1$risk_60to64*0.125715739 + ECA_PRS_Average_Data1$risk_65to69*0.109659337 + ECA_PRS_Average_Data1$risk_70to74*0.086579466)

ECA_PRS_AR_CurrentSS=ECA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
esophageal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
ECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
ECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
ECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
ECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
ECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
ECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
ECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
ECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
ECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

ECA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(ECA_PRS_30to34$risk), as.numeric(ECA_PRS_35to39$risk), as.numeric(ECA_PRS_40to44$risk), as.numeric(ECA_PRS_45to49$risk), 
                                          as.numeric(ECA_PRS_50to54$risk), as.numeric(ECA_PRS_55to59$risk), as.numeric(ECA_PRS_60to64$risk), as.numeric(ECA_PRS_65to69$risk),
                                          as.numeric(ECA_PRS_70to74$risk)))

colnames(ECA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

ECA_PRS_Average_Data2$risk_wtd_avg=(ECA_PRS_Average_Data2$risk_30to34*0.108040659 + ECA_PRS_Average_Data2$risk_35to39*0.104987889 + ECA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                    ECA_PRS_Average_Data2$risk_45to49*0.112278705 + ECA_PRS_Average_Data2$risk_50to54*0.121739409 + ECA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                    ECA_PRS_Average_Data2$risk_60to64*0.125715739 + ECA_PRS_Average_Data2$risk_65to69*0.109659337 + ECA_PRS_Average_Data2$risk_70to74*0.086579466)

ECA_PRS_AR_2xSS=ECA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
esophageal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
ECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
ECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
ECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
ECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
ECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
ECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
ECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
ECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
ECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

ECA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(ECA_PRS_30to34$risk), as.numeric(ECA_PRS_35to39$risk), as.numeric(ECA_PRS_40to44$risk), as.numeric(ECA_PRS_45to49$risk), 
                                          as.numeric(ECA_PRS_50to54$risk), as.numeric(ECA_PRS_55to59$risk), as.numeric(ECA_PRS_60to64$risk), as.numeric(ECA_PRS_65to69$risk),
                                          as.numeric(ECA_PRS_70to74$risk)))

colnames(ECA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

ECA_PRS_Average_Data3$risk_wtd_avg=(ECA_PRS_Average_Data3$risk_30to34*0.108040659 + ECA_PRS_Average_Data3$risk_35to39*0.104987889 + ECA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                    ECA_PRS_Average_Data3$risk_45to49*0.112278705 + ECA_PRS_Average_Data3$risk_50to54*0.121739409 + ECA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                    ECA_PRS_Average_Data3$risk_60to64*0.125715739 + ECA_PRS_Average_Data3$risk_65to69*0.109659337 + ECA_PRS_Average_Data3$risk_70to74*0.086579466)

ECA_PRS_AR_4xSS=ECA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
esophageal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
ECA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
ECA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
ECA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
ECA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                                      model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
ECA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
ECA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
ECA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
ECA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
ECA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = esophageal_ref,
                           model.disease.incidence.rates = esophageal_inc, model.competing.incidence.rates = mort_no_esophageal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = esophageal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

ECA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(ECA_PRS_30to34$risk), as.numeric(ECA_PRS_35to39$risk), as.numeric(ECA_PRS_40to44$risk), as.numeric(ECA_PRS_45to49$risk), 
                                          as.numeric(ECA_PRS_50to54$risk), as.numeric(ECA_PRS_55to59$risk), as.numeric(ECA_PRS_60to64$risk), as.numeric(ECA_PRS_65to69$risk),
                                          as.numeric(ECA_PRS_70to74$risk)))

colnames(ECA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

ECA_PRS_Average_Data4$risk_wtd_avg=(ECA_PRS_Average_Data4$risk_30to34*0.108040659 + ECA_PRS_Average_Data4$risk_35to39*0.104987889 + ECA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                    ECA_PRS_Average_Data4$risk_45to49*0.112278705 + ECA_PRS_Average_Data4$risk_50to54*0.121739409 + ECA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                    ECA_PRS_Average_Data4$risk_60to64*0.125715739 + ECA_PRS_Average_Data4$risk_65to69*0.109659337 + ECA_PRS_Average_Data4$risk_70to74*0.086579466)

ECA_PRS_AR_InfSS=ECA_PRS_Average_Data4$risk_wtd_avg

save(ECA_PRS_AR_CurrentSS,ECA_PRS_AR_2xSS,ECA_PRS_AR_4xSS,ECA_PRS_AR_InfSS, file='ECA_PRS_AR.Rda')


##### AGE-STRATIFIED PLOTTING DATA #####
ECA_PRS_AR_30to34=as.data.frame(cbind(ECA_PRS_Average_Data1$risk_30to34, ECA_PRS_Average_Data2$risk_30to34,
                                      ECA_PRS_Average_Data3$risk_30to34, ECA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(ECA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

ECA_PRS_AR_40to44=as.data.frame(cbind(ECA_PRS_Average_Data1$risk_40to44, ECA_PRS_Average_Data2$risk_40to44,
                                      ECA_PRS_Average_Data3$risk_40to44, ECA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(ECA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

ECA_PRS_AR_50to54=as.data.frame(cbind(ECA_PRS_Average_Data1$risk_50to54, ECA_PRS_Average_Data2$risk_50to54,
                                      ECA_PRS_Average_Data3$risk_50to54, ECA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(ECA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

ECA_PRS_AR_60to64=as.data.frame(cbind(ECA_PRS_Average_Data1$risk_60to64, ECA_PRS_Average_Data2$risk_60to64,
                                      ECA_PRS_Average_Data3$risk_60to64, ECA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(ECA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(ECA_PRS_AR_30to34, ECA_PRS_AR_40to44, ECA_PRS_AR_50to54, ECA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/ECA_PRS_AR_AgeStrat.Rda')

      ############################
      # Glioma cancer
      ############################
var1=0.438219546
var2=0.472372237
var3=0.533605282
var4=0.874251941
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
glioma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
GCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
GCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
GCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
GCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
GCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
GCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
GCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
GCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
GCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

GCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(GCA_PRS_30to34$risk), as.numeric(GCA_PRS_35to39$risk), as.numeric(GCA_PRS_40to44$risk), as.numeric(GCA_PRS_45to49$risk),
                                           as.numeric(GCA_PRS_50to54$risk), as.numeric(GCA_PRS_55to59$risk), as.numeric(GCA_PRS_60to64$risk), as.numeric(GCA_PRS_65to69$risk),
                                           as.numeric(GCA_PRS_70to74$risk)))

colnames(GCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

GCA_PRS_Average_Data1$risk_wtd_avg=(GCA_PRS_Average_Data1$risk_30to34*0.108040659 + GCA_PRS_Average_Data1$risk_35to39*0.104987889 + GCA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                     GCA_PRS_Average_Data1$risk_45to49*0.112278705 + GCA_PRS_Average_Data1$risk_50to54*0.121739409 + GCA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                     GCA_PRS_Average_Data1$risk_60to64*0.125715739 + GCA_PRS_Average_Data1$risk_65to69*0.109659337 + GCA_PRS_Average_Data1$risk_70to74*0.086579466)

GCA_PRS_AR_CurrentSS=GCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
glioma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
GCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
GCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
GCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
GCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
GCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
GCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
GCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
GCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
GCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

GCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(GCA_PRS_30to34$risk), as.numeric(GCA_PRS_35to39$risk), as.numeric(GCA_PRS_40to44$risk), as.numeric(GCA_PRS_45to49$risk),
                                           as.numeric(GCA_PRS_50to54$risk), as.numeric(GCA_PRS_55to59$risk), as.numeric(GCA_PRS_60to64$risk), as.numeric(GCA_PRS_65to69$risk),
                                           as.numeric(GCA_PRS_70to74$risk)))

colnames(GCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

GCA_PRS_Average_Data2$risk_wtd_avg=(GCA_PRS_Average_Data2$risk_30to34*0.108040659 + GCA_PRS_Average_Data2$risk_35to39*0.104987889 + GCA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                     GCA_PRS_Average_Data2$risk_45to49*0.112278705 + GCA_PRS_Average_Data2$risk_50to54*0.121739409 + GCA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                     GCA_PRS_Average_Data2$risk_60to64*0.125715739 + GCA_PRS_Average_Data2$risk_65to69*0.109659337 + GCA_PRS_Average_Data2$risk_70to74*0.086579466)

GCA_PRS_AR_2xSS=GCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
glioma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
GCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
GCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
GCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
GCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
GCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
GCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
GCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
GCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
GCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

GCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(GCA_PRS_30to34$risk), as.numeric(GCA_PRS_35to39$risk), as.numeric(GCA_PRS_40to44$risk), as.numeric(GCA_PRS_45to49$risk),
                                           as.numeric(GCA_PRS_50to54$risk), as.numeric(GCA_PRS_55to59$risk), as.numeric(GCA_PRS_60to64$risk), as.numeric(GCA_PRS_65to69$risk),
                                           as.numeric(GCA_PRS_70to74$risk)))

colnames(GCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

GCA_PRS_Average_Data3$risk_wtd_avg=(GCA_PRS_Average_Data3$risk_30to34*0.108040659 + GCA_PRS_Average_Data3$risk_35to39*0.104987889 + GCA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                     GCA_PRS_Average_Data3$risk_45to49*0.112278705 + GCA_PRS_Average_Data3$risk_50to54*0.121739409 + GCA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                     GCA_PRS_Average_Data3$risk_60to64*0.125715739 + GCA_PRS_Average_Data3$risk_65to69*0.109659337 + GCA_PRS_Average_Data3$risk_70to74*0.086579466)

GCA_PRS_AR_4xSS=GCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
glioma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
GCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
GCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
GCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
GCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
GCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
GCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
GCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
GCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
GCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = glioma_ref,
                           model.disease.incidence.rates = glioma_inc, model.competing.incidence.rates = mort_no_glioma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = glioma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

GCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(GCA_PRS_30to34$risk), as.numeric(GCA_PRS_35to39$risk), as.numeric(GCA_PRS_40to44$risk), as.numeric(GCA_PRS_45to49$risk),
                                           as.numeric(GCA_PRS_50to54$risk), as.numeric(GCA_PRS_55to59$risk), as.numeric(GCA_PRS_60to64$risk), as.numeric(GCA_PRS_65to69$risk),
                                           as.numeric(GCA_PRS_70to74$risk)))

colnames(GCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

GCA_PRS_Average_Data4$risk_wtd_avg=(GCA_PRS_Average_Data4$risk_30to34*0.108040659 + GCA_PRS_Average_Data4$risk_35to39*0.104987889 + GCA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                     GCA_PRS_Average_Data4$risk_45to49*0.112278705 + GCA_PRS_Average_Data4$risk_50to54*0.121739409 + GCA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                     GCA_PRS_Average_Data4$risk_60to64*0.125715739 + GCA_PRS_Average_Data4$risk_65to69*0.109659337 + GCA_PRS_Average_Data4$risk_70to74*0.086579466)

GCA_PRS_AR_InfSS=GCA_PRS_Average_Data4$risk_wtd_avg

save(GCA_PRS_AR_CurrentSS,GCA_PRS_AR_2xSS,GCA_PRS_AR_4xSS,GCA_PRS_AR_InfSS, file='GCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
GCA_PRS_AR_30to34=as.data.frame(cbind(GCA_PRS_Average_Data1$risk_30to34, GCA_PRS_Average_Data2$risk_30to34,
                                      GCA_PRS_Average_Data3$risk_30to34, GCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(GCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

GCA_PRS_AR_40to44=as.data.frame(cbind(GCA_PRS_Average_Data1$risk_40to44, GCA_PRS_Average_Data2$risk_40to44,
                                      GCA_PRS_Average_Data3$risk_40to44, GCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(GCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

GCA_PRS_AR_50to54=as.data.frame(cbind(GCA_PRS_Average_Data1$risk_50to54, GCA_PRS_Average_Data2$risk_50to54,
                                      GCA_PRS_Average_Data3$risk_50to54, GCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(GCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

GCA_PRS_AR_60to64=as.data.frame(cbind(GCA_PRS_Average_Data1$risk_60to64, GCA_PRS_Average_Data2$risk_60to64,
                                      GCA_PRS_Average_Data3$risk_60to64, GCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(GCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(GCA_PRS_AR_30to34, GCA_PRS_AR_40to44, GCA_PRS_AR_50to54, GCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/GCA_PRS_AR_AgeStrat.Rda')

      ############################
      # Oral & Pharyngeal cancer
      ############################
var1=0.015616731
var2=0.091552096
var3=0.262511398
var4=0.683084386
##### 50K SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
headneck_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
HNC_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
HNC_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
HNC_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
HNC_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
HNC_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=52,
                                     apply.age.interval.length = 23, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
HNC_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=57,
                                     apply.age.interval.length = 18, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
HNC_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=62,
                                     apply.age.interval.length = 13, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
HNC_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=67,
                                     apply.age.interval.length = 8, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
HNC_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=72,
                                     apply.age.interval.length = 3, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

HNC_PRS_Average_Data1=as.data.frame(cbind(as.numeric(HNC_PRS_30to34$risk), as.numeric(HNC_PRS_35to39$risk), as.numeric(HNC_PRS_40to44$risk), as.numeric(HNC_PRS_45to49$risk),
                                           as.numeric(HNC_PRS_50to54$risk), as.numeric(HNC_PRS_55to59$risk), as.numeric(HNC_PRS_60to64$risk), as.numeric(HNC_PRS_65to69$risk),
                                           as.numeric(HNC_PRS_70to74$risk)))

colnames(HNC_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

HNC_PRS_Average_Data1$risk_wtd_avg=(HNC_PRS_Average_Data1$risk_30to34*0.108040659 + HNC_PRS_Average_Data1$risk_35to39*0.104987889 + HNC_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                     HNC_PRS_Average_Data1$risk_45to49*0.112278705 + HNC_PRS_Average_Data1$risk_50to54*0.121739409 + HNC_PRS_Average_Data1$risk_55to59*0.133239483 +
                                     HNC_PRS_Average_Data1$risk_60to64*0.125715739 + HNC_PRS_Average_Data1$risk_65to69*0.109659337 + HNC_PRS_Average_Data1$risk_70to74*0.086579466)

HNC_PRS_AR_50K=HNC_PRS_Average_Data1$risk_wtd_avg

##### 100K SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
headneck_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
HNC_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
HNC_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
HNC_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
HNC_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
HNC_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=52,
                                     apply.age.interval.length = 23, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
HNC_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=57,
                                     apply.age.interval.length = 18, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
HNC_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=62,
                                     apply.age.interval.length = 13, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
HNC_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=67,
                                     apply.age.interval.length = 8, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
HNC_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=72,
                                     apply.age.interval.length = 3, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

HNC_PRS_Average_Data2=as.data.frame(cbind(as.numeric(HNC_PRS_30to34$risk), as.numeric(HNC_PRS_35to39$risk), as.numeric(HNC_PRS_40to44$risk), as.numeric(HNC_PRS_45to49$risk),
                                           as.numeric(HNC_PRS_50to54$risk), as.numeric(HNC_PRS_55to59$risk), as.numeric(HNC_PRS_60to64$risk), as.numeric(HNC_PRS_65to69$risk),
                                           as.numeric(HNC_PRS_70to74$risk)))

colnames(HNC_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

HNC_PRS_Average_Data2$risk_wtd_avg=(HNC_PRS_Average_Data2$risk_30to34*0.108040659 + HNC_PRS_Average_Data2$risk_35to39*0.104987889 + HNC_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                     HNC_PRS_Average_Data2$risk_45to49*0.112278705 + HNC_PRS_Average_Data2$risk_50to54*0.121739409 + HNC_PRS_Average_Data2$risk_55to59*0.133239483 +
                                     HNC_PRS_Average_Data2$risk_60to64*0.125715739 + HNC_PRS_Average_Data2$risk_65to69*0.109659337 + HNC_PRS_Average_Data2$risk_70to74*0.086579466)

HNC_PRS_AR_100K=HNC_PRS_Average_Data2$risk_wtd_avg

##### 200K SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
headneck_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
HNC_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
HNC_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
HNC_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
HNC_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
HNC_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=52,
                                     apply.age.interval.length = 23, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
HNC_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=57,
                                     apply.age.interval.length = 18, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
HNC_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=62,
                                     apply.age.interval.length = 13, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
HNC_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=67,
                                     apply.age.interval.length = 8, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
HNC_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=72,
                                     apply.age.interval.length = 3, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

HNC_PRS_Average_Data3=as.data.frame(cbind(as.numeric(HNC_PRS_30to34$risk), as.numeric(HNC_PRS_35to39$risk), as.numeric(HNC_PRS_40to44$risk), as.numeric(HNC_PRS_45to49$risk),
                                           as.numeric(HNC_PRS_50to54$risk), as.numeric(HNC_PRS_55to59$risk), as.numeric(HNC_PRS_60to64$risk), as.numeric(HNC_PRS_65to69$risk),
                                           as.numeric(HNC_PRS_70to74$risk)))

colnames(HNC_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

HNC_PRS_Average_Data3$risk_wtd_avg=(HNC_PRS_Average_Data3$risk_30to34*0.108040659 + HNC_PRS_Average_Data3$risk_35to39*0.104987889 + HNC_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                     HNC_PRS_Average_Data3$risk_45to49*0.112278705 + HNC_PRS_Average_Data3$risk_50to54*0.121739409 + HNC_PRS_Average_Data3$risk_55to59*0.133239483 +
                                     HNC_PRS_Average_Data3$risk_60to64*0.125715739 + HNC_PRS_Average_Data3$risk_65to69*0.109659337 + HNC_PRS_Average_Data3$risk_70to74*0.086579466)

HNC_PRS_AR_200K=HNC_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
headneck_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
HNC_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
HNC_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
HNC_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
HNC_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                           model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
HNC_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=52,
                                     apply.age.interval.length = 23, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
HNC_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=57,
                                     apply.age.interval.length = 18, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
HNC_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=62,
                                     apply.age.interval.length = 13, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
HNC_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=67,
                                     apply.age.interval.length = 8, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
HNC_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = headneck_ref,
                                     model.disease.incidence.rates = headneck_inc, model.competing.incidence.rates = mort_no_headneck, model.bin.fh.name=NULL, apply.age.start=72,
                                     apply.age.interval.length = 3, apply.cov.profile = headneck_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

HNC_PRS_Average_Data4=as.data.frame(cbind(as.numeric(HNC_PRS_30to34$risk), as.numeric(HNC_PRS_35to39$risk), as.numeric(HNC_PRS_40to44$risk), as.numeric(HNC_PRS_45to49$risk),
                                           as.numeric(HNC_PRS_50to54$risk), as.numeric(HNC_PRS_55to59$risk), as.numeric(HNC_PRS_60to64$risk), as.numeric(HNC_PRS_65to69$risk),
                                           as.numeric(HNC_PRS_70to74$risk)))

colnames(HNC_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

HNC_PRS_Average_Data4$risk_wtd_avg=(HNC_PRS_Average_Data4$risk_30to34*0.108040659 + HNC_PRS_Average_Data4$risk_35to39*0.104987889 + HNC_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                     HNC_PRS_Average_Data4$risk_45to49*0.112278705 + HNC_PRS_Average_Data4$risk_50to54*0.121739409 + HNC_PRS_Average_Data4$risk_55to59*0.133239483 +
                                     HNC_PRS_Average_Data4$risk_60to64*0.125715739 + HNC_PRS_Average_Data4$risk_65to69*0.109659337 + HNC_PRS_Average_Data4$risk_70to74*0.086579466)

HNC_PRS_AR_InfSS=HNC_PRS_Average_Data4$risk_wtd_avg

save(HNC_PRS_AR_50K,HNC_PRS_AR_100K,HNC_PRS_AR_200K,HNC_PRS_AR_InfSS, file='HNC_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
HNC_PRS_AR_30to34=as.data.frame(cbind(HNC_PRS_Average_Data1$risk_30to34, HNC_PRS_Average_Data2$risk_30to34,
                                      HNC_PRS_Average_Data3$risk_30to34, HNC_PRS_Average_Data4$risk_30to34)*100)
                                colnames(HNC_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

HNC_PRS_AR_40to44=as.data.frame(cbind(HNC_PRS_Average_Data1$risk_40to44, HNC_PRS_Average_Data2$risk_40to44,
                                      HNC_PRS_Average_Data3$risk_40to44, HNC_PRS_Average_Data4$risk_40to44)*100)
                                colnames(HNC_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

HNC_PRS_AR_50to54=as.data.frame(cbind(HNC_PRS_Average_Data1$risk_50to54, HNC_PRS_Average_Data2$risk_50to54,
                                      HNC_PRS_Average_Data3$risk_50to54, HNC_PRS_Average_Data4$risk_50to54)*100)
                                colnames(HNC_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

HNC_PRS_AR_60to64=as.data.frame(cbind(HNC_PRS_Average_Data1$risk_60to64, HNC_PRS_Average_Data2$risk_60to64,
                                      HNC_PRS_Average_Data3$risk_60to64, HNC_PRS_Average_Data4$risk_60to64)*100)
                                colnames(HNC_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(HNC_PRS_AR_30to34, HNC_PRS_AR_40to44, HNC_PRS_AR_50to54, HNC_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/HNC_PRS_AR_AgeStrat.Rda')

      ############################
      # Lung cancer
      ############################
var1=0.056350792
var2=0.063400325
var3=0.09409151
var4=0.391382459
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
lung_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
LCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
LCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
LCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
LCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
LCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
LCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
LCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
LCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
LCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

LCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(LCA_PRS_30to34$risk), as.numeric(LCA_PRS_35to39$risk), as.numeric(LCA_PRS_40to44$risk), as.numeric(LCA_PRS_45to49$risk),
                                          as.numeric(LCA_PRS_50to54$risk), as.numeric(LCA_PRS_55to59$risk), as.numeric(LCA_PRS_60to64$risk), as.numeric(LCA_PRS_65to69$risk),
                                          as.numeric(LCA_PRS_70to74$risk)))

colnames(LCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

LCA_PRS_Average_Data1$risk_wtd_avg=(LCA_PRS_Average_Data1$risk_30to34*0.108040659 + LCA_PRS_Average_Data1$risk_35to39*0.104987889 + LCA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                      LCA_PRS_Average_Data1$risk_45to49*0.112278705 + LCA_PRS_Average_Data1$risk_50to54*0.121739409 + LCA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                      LCA_PRS_Average_Data1$risk_60to64*0.125715739 + LCA_PRS_Average_Data1$risk_65to69*0.109659337 + LCA_PRS_Average_Data1$risk_70to74*0.086579466)

LCA_PRS_AR_CurrentSS=LCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
lung_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
LCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
LCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
LCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
LCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
LCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
LCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
LCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
LCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
LCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

LCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(LCA_PRS_30to34$risk), as.numeric(LCA_PRS_35to39$risk), as.numeric(LCA_PRS_40to44$risk), as.numeric(LCA_PRS_45to49$risk),
                                          as.numeric(LCA_PRS_50to54$risk), as.numeric(LCA_PRS_55to59$risk), as.numeric(LCA_PRS_60to64$risk), as.numeric(LCA_PRS_65to69$risk),
                                          as.numeric(LCA_PRS_70to74$risk)))

colnames(LCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

LCA_PRS_Average_Data2$risk_wtd_avg=(LCA_PRS_Average_Data2$risk_30to34*0.108040659 + LCA_PRS_Average_Data2$risk_35to39*0.104987889 + LCA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                      LCA_PRS_Average_Data2$risk_45to49*0.112278705 + LCA_PRS_Average_Data2$risk_50to54*0.121739409 + LCA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                      LCA_PRS_Average_Data2$risk_60to64*0.125715739 + LCA_PRS_Average_Data2$risk_65to69*0.109659337 + LCA_PRS_Average_Data2$risk_70to74*0.086579466)

LCA_PRS_AR_2xSS=LCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
lung_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
LCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
LCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
LCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
LCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
LCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
LCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
LCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
LCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
LCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

LCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(LCA_PRS_30to34$risk), as.numeric(LCA_PRS_35to39$risk), as.numeric(LCA_PRS_40to44$risk), as.numeric(LCA_PRS_45to49$risk),
                                          as.numeric(LCA_PRS_50to54$risk), as.numeric(LCA_PRS_55to59$risk), as.numeric(LCA_PRS_60to64$risk), as.numeric(LCA_PRS_65to69$risk),
                                          as.numeric(LCA_PRS_70to74$risk)))

colnames(LCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

LCA_PRS_Average_Data3$risk_wtd_avg=(LCA_PRS_Average_Data3$risk_30to34*0.108040659 + LCA_PRS_Average_Data3$risk_35to39*0.104987889 + LCA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                      LCA_PRS_Average_Data3$risk_45to49*0.112278705 + LCA_PRS_Average_Data3$risk_50to54*0.121739409 + LCA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                      LCA_PRS_Average_Data3$risk_60to64*0.125715739 + LCA_PRS_Average_Data3$risk_65to69*0.109659337 + LCA_PRS_Average_Data3$risk_70to74*0.086579466)

LCA_PRS_AR_4xSS=LCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
lung_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
LCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
LCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
LCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
LCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                                     model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
LCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
LCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
LCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
LCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
LCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lung_ref,
                           model.disease.incidence.rates = lung_inc, model.competing.incidence.rates = mort_no_lung, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lung_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

LCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(LCA_PRS_30to34$risk), as.numeric(LCA_PRS_35to39$risk), as.numeric(LCA_PRS_40to44$risk), as.numeric(LCA_PRS_45to49$risk),
                                          as.numeric(LCA_PRS_50to54$risk), as.numeric(LCA_PRS_55to59$risk), as.numeric(LCA_PRS_60to64$risk), as.numeric(LCA_PRS_65to69$risk),
                                          as.numeric(LCA_PRS_70to74$risk)))

colnames(LCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

LCA_PRS_Average_Data4$risk_wtd_avg=(LCA_PRS_Average_Data4$risk_30to34*0.108040659 + LCA_PRS_Average_Data4$risk_35to39*0.104987889 + LCA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                      LCA_PRS_Average_Data4$risk_45to49*0.112278705 + LCA_PRS_Average_Data4$risk_50to54*0.121739409 + LCA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                      LCA_PRS_Average_Data4$risk_60to64*0.125715739 + LCA_PRS_Average_Data4$risk_65to69*0.109659337 + LCA_PRS_Average_Data4$risk_70to74*0.086579466)

LCA_PRS_AR_InfSS=LCA_PRS_Average_Data4$risk_wtd_avg

save(LCA_PRS_AR_CurrentSS,LCA_PRS_AR_2xSS,LCA_PRS_AR_4xSS,LCA_PRS_AR_InfSS, file='LCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
LCA_PRS_AR_30to34=as.data.frame(cbind(LCA_PRS_Average_Data1$risk_30to34, LCA_PRS_Average_Data2$risk_30to34,
                                      LCA_PRS_Average_Data3$risk_30to34, LCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(LCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

LCA_PRS_AR_40to44=as.data.frame(cbind(LCA_PRS_Average_Data1$risk_40to44, LCA_PRS_Average_Data2$risk_40to44,
                                      LCA_PRS_Average_Data3$risk_40to44, LCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(LCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

LCA_PRS_AR_50to54=as.data.frame(cbind(LCA_PRS_Average_Data1$risk_50to54, LCA_PRS_Average_Data2$risk_50to54,
                                      LCA_PRS_Average_Data3$risk_50to54, LCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(LCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

LCA_PRS_AR_60to64=as.data.frame(cbind(LCA_PRS_Average_Data1$risk_60to64, LCA_PRS_Average_Data2$risk_60to64,
                                      LCA_PRS_Average_Data3$risk_60to64, LCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(LCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(LCA_PRS_AR_30to34, LCA_PRS_AR_40to44, LCA_PRS_AR_50to54, LCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/LCA_PRS_AR_AgeStrat.Rda')

      ############################
      # Lymphoma (CLL)
      ############################
var1=0.527055781
var2=0.620917628
var3=0.747588098
var4=1.625033989
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
lymphCLL_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CLL_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CLL_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CLL_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CLL_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CLL_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CLL_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CLL_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CLL_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CLL_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CLL_PRS_Average_Data1=as.data.frame(cbind(as.numeric(CLL_PRS_30to34$risk), as.numeric(CLL_PRS_35to39$risk), as.numeric(CLL_PRS_40to44$risk), as.numeric(CLL_PRS_45to49$risk),
                                          as.numeric(CLL_PRS_50to54$risk), as.numeric(CLL_PRS_55to59$risk), as.numeric(CLL_PRS_60to64$risk), as.numeric(CLL_PRS_65to69$risk),
                                          as.numeric(CLL_PRS_70to74$risk)))

colnames(CLL_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CLL_PRS_Average_Data1$risk_wtd_avg=(CLL_PRS_Average_Data1$risk_30to34*0.108040659 + CLL_PRS_Average_Data1$risk_35to39*0.104987889 + CLL_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                      CLL_PRS_Average_Data1$risk_45to49*0.112278705 + CLL_PRS_Average_Data1$risk_50to54*0.121739409 + CLL_PRS_Average_Data1$risk_55to59*0.133239483 +
                                      CLL_PRS_Average_Data1$risk_60to64*0.125715739 + CLL_PRS_Average_Data1$risk_65to69*0.109659337 + CLL_PRS_Average_Data1$risk_70to74*0.086579466)

CLL_PRS_AR_CurrentSS=CLL_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
lymphCLL_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CLL_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CLL_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CLL_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CLL_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CLL_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CLL_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CLL_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CLL_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CLL_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CLL_PRS_Average_Data2=as.data.frame(cbind(as.numeric(CLL_PRS_30to34$risk), as.numeric(CLL_PRS_35to39$risk), as.numeric(CLL_PRS_40to44$risk), as.numeric(CLL_PRS_45to49$risk),
                                          as.numeric(CLL_PRS_50to54$risk), as.numeric(CLL_PRS_55to59$risk), as.numeric(CLL_PRS_60to64$risk), as.numeric(CLL_PRS_65to69$risk),
                                          as.numeric(CLL_PRS_70to74$risk)))

colnames(CLL_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CLL_PRS_Average_Data2$risk_wtd_avg=(CLL_PRS_Average_Data2$risk_30to34*0.108040659 + CLL_PRS_Average_Data2$risk_35to39*0.104987889 + CLL_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                      CLL_PRS_Average_Data2$risk_45to49*0.112278705 + CLL_PRS_Average_Data2$risk_50to54*0.121739409 + CLL_PRS_Average_Data2$risk_55to59*0.133239483 +
                                      CLL_PRS_Average_Data2$risk_60to64*0.125715739 + CLL_PRS_Average_Data2$risk_65to69*0.109659337 + CLL_PRS_Average_Data2$risk_70to74*0.086579466)

CLL_PRS_AR_2xSS=CLL_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
lymphCLL_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CLL_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CLL_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CLL_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CLL_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CLL_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CLL_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CLL_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CLL_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CLL_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CLL_PRS_Average_Data3=as.data.frame(cbind(as.numeric(CLL_PRS_30to34$risk), as.numeric(CLL_PRS_35to39$risk), as.numeric(CLL_PRS_40to44$risk), as.numeric(CLL_PRS_45to49$risk),
                                          as.numeric(CLL_PRS_50to54$risk), as.numeric(CLL_PRS_55to59$risk), as.numeric(CLL_PRS_60to64$risk), as.numeric(CLL_PRS_65to69$risk),
                                          as.numeric(CLL_PRS_70to74$risk)))

colnames(CLL_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CLL_PRS_Average_Data3$risk_wtd_avg=(CLL_PRS_Average_Data3$risk_30to34*0.108040659 + CLL_PRS_Average_Data3$risk_35to39*0.104987889 + CLL_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                      CLL_PRS_Average_Data3$risk_45to49*0.112278705 + CLL_PRS_Average_Data3$risk_50to54*0.121739409 + CLL_PRS_Average_Data3$risk_55to59*0.133239483 +
                                      CLL_PRS_Average_Data3$risk_60to64*0.125715739 + CLL_PRS_Average_Data3$risk_65to69*0.109659337 + CLL_PRS_Average_Data3$risk_70to74*0.086579466)

CLL_PRS_AR_4xSS=CLL_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
lymphCLL_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
CLL_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
CLL_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
CLL_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
CLL_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                                     model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
CLL_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
CLL_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
CLL_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
CLL_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
CLL_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = lymphCLL_ref,
                           model.disease.incidence.rates = lymphCLL_inc, model.competing.incidence.rates = mort_no_lymphCLL, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = lymphCLL_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

CLL_PRS_Average_Data4=as.data.frame(cbind(as.numeric(CLL_PRS_30to34$risk), as.numeric(CLL_PRS_35to39$risk), as.numeric(CLL_PRS_40to44$risk), as.numeric(CLL_PRS_45to49$risk),
                                          as.numeric(CLL_PRS_50to54$risk), as.numeric(CLL_PRS_55to59$risk), as.numeric(CLL_PRS_60to64$risk), as.numeric(CLL_PRS_65to69$risk),
                                          as.numeric(CLL_PRS_70to74$risk)))

colnames(CLL_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

CLL_PRS_Average_Data4$risk_wtd_avg=(CLL_PRS_Average_Data4$risk_30to34*0.108040659 + CLL_PRS_Average_Data4$risk_35to39*0.104987889 + CLL_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                      CLL_PRS_Average_Data4$risk_45to49*0.112278705 + CLL_PRS_Average_Data4$risk_50to54*0.121739409 + CLL_PRS_Average_Data4$risk_55to59*0.133239483 +
                                      CLL_PRS_Average_Data4$risk_60to64*0.125715739 + CLL_PRS_Average_Data4$risk_65to69*0.109659337 + CLL_PRS_Average_Data4$risk_70to74*0.086579466)

CLL_PRS_AR_InfSS=CLL_PRS_Average_Data4$risk_wtd_avg

save(CLL_PRS_AR_CurrentSS,CLL_PRS_AR_2xSS,CLL_PRS_AR_4xSS,CLL_PRS_AR_InfSS, file='CLL_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
CLL_PRS_AR_30to34=as.data.frame(cbind(CLL_PRS_Average_Data1$risk_30to34, CLL_PRS_Average_Data2$risk_30to34,
                                      CLL_PRS_Average_Data3$risk_30to34, CLL_PRS_Average_Data4$risk_30to34)*100)
                                colnames(CLL_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CLL_PRS_AR_40to44=as.data.frame(cbind(CLL_PRS_Average_Data1$risk_40to44, CLL_PRS_Average_Data2$risk_40to44,
                                      CLL_PRS_Average_Data3$risk_40to44, CLL_PRS_Average_Data4$risk_40to44)*100)
                                colnames(CLL_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CLL_PRS_AR_50to54=as.data.frame(cbind(CLL_PRS_Average_Data1$risk_50to54, CLL_PRS_Average_Data2$risk_50to54,
                                      CLL_PRS_Average_Data3$risk_50to54, CLL_PRS_Average_Data4$risk_50to54)*100)
                                colnames(CLL_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

CLL_PRS_AR_60to64=as.data.frame(cbind(CLL_PRS_Average_Data1$risk_60to64, CLL_PRS_Average_Data2$risk_60to64,
                                      CLL_PRS_Average_Data3$risk_60to64, CLL_PRS_Average_Data4$risk_60to64)*100)
                                colnames(CLL_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(CLL_PRS_AR_30to34, CLL_PRS_AR_40to44, CLL_PRS_AR_50to54, CLL_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/CLL_PRS_AR_AgeStrat.Rda')

      ############################
      # Melanoma skin cancer
      ############################
var1=0.271701543
var2=0.351471353
var3=0.455525056
var4=0.650217714
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
melanoma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
MSCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
MSCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
MSCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
MSCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
MSCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
MSCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
MSCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
MSCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
MSCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

MSCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(MSCA_PRS_30to34$risk), as.numeric(MSCA_PRS_35to39$risk), as.numeric(MSCA_PRS_40to44$risk), as.numeric(MSCA_PRS_45to49$risk),
                                          as.numeric(MSCA_PRS_50to54$risk), as.numeric(MSCA_PRS_55to59$risk), as.numeric(MSCA_PRS_60to64$risk), as.numeric(MSCA_PRS_65to69$risk),
                                          as.numeric(MSCA_PRS_70to74$risk)))

colnames(MSCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

MSCA_PRS_Average_Data1$risk_wtd_avg=(MSCA_PRS_Average_Data1$risk_30to34*0.108040659 + MSCA_PRS_Average_Data1$risk_35to39*0.104987889 + MSCA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                      MSCA_PRS_Average_Data1$risk_45to49*0.112278705 + MSCA_PRS_Average_Data1$risk_50to54*0.121739409 + MSCA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                      MSCA_PRS_Average_Data1$risk_60to64*0.125715739 + MSCA_PRS_Average_Data1$risk_65to69*0.109659337 + MSCA_PRS_Average_Data1$risk_70to74*0.086579466)

MSCA_PRS_AR_CurrentSS=MSCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
melanoma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
MSCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
MSCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
MSCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
MSCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
MSCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
MSCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
MSCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
MSCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
MSCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

MSCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(MSCA_PRS_30to34$risk), as.numeric(MSCA_PRS_35to39$risk), as.numeric(MSCA_PRS_40to44$risk), as.numeric(MSCA_PRS_45to49$risk),
                                          as.numeric(MSCA_PRS_50to54$risk), as.numeric(MSCA_PRS_55to59$risk), as.numeric(MSCA_PRS_60to64$risk), as.numeric(MSCA_PRS_65to69$risk),
                                          as.numeric(MSCA_PRS_70to74$risk)))

colnames(MSCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

MSCA_PRS_Average_Data2$risk_wtd_avg=(MSCA_PRS_Average_Data2$risk_30to34*0.108040659 + MSCA_PRS_Average_Data2$risk_35to39*0.104987889 + MSCA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                      MSCA_PRS_Average_Data2$risk_45to49*0.112278705 + MSCA_PRS_Average_Data2$risk_50to54*0.121739409 + MSCA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                      MSCA_PRS_Average_Data2$risk_60to64*0.125715739 + MSCA_PRS_Average_Data2$risk_65to69*0.109659337 + MSCA_PRS_Average_Data2$risk_70to74*0.086579466)

MSCA_PRS_AR_2xSS=MSCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
melanoma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
MSCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
MSCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
MSCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
MSCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
MSCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
MSCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
MSCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
MSCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
MSCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

MSCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(MSCA_PRS_30to34$risk), as.numeric(MSCA_PRS_35to39$risk), as.numeric(MSCA_PRS_40to44$risk), as.numeric(MSCA_PRS_45to49$risk),
                                          as.numeric(MSCA_PRS_50to54$risk), as.numeric(MSCA_PRS_55to59$risk), as.numeric(MSCA_PRS_60to64$risk), as.numeric(MSCA_PRS_65to69$risk),
                                          as.numeric(MSCA_PRS_70to74$risk)))

colnames(MSCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

MSCA_PRS_Average_Data3$risk_wtd_avg=(MSCA_PRS_Average_Data3$risk_30to34*0.108040659 + MSCA_PRS_Average_Data3$risk_35to39*0.104987889 + MSCA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                      MSCA_PRS_Average_Data3$risk_45to49*0.112278705 + MSCA_PRS_Average_Data3$risk_50to54*0.121739409 + MSCA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                      MSCA_PRS_Average_Data3$risk_60to64*0.125715739 + MSCA_PRS_Average_Data3$risk_65to69*0.109659337 + MSCA_PRS_Average_Data3$risk_70to74*0.086579466)

MSCA_PRS_AR_4xSS=MSCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
melanoma_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
MSCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=32,
                                     apply.age.interval.length = 43, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
MSCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=37,
                                     apply.age.interval.length = 38, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
MSCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=42,
                                     apply.age.interval.length = 33, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
MSCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                                     model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=47,
                                     apply.age.interval.length = 28, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
MSCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
MSCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
MSCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
MSCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
MSCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = melanoma_ref,
                           model.disease.incidence.rates = melanoma_inc, model.competing.incidence.rates = mort_no_melanoma, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = melanoma_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

MSCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(MSCA_PRS_30to34$risk), as.numeric(MSCA_PRS_35to39$risk), as.numeric(MSCA_PRS_40to44$risk), as.numeric(MSCA_PRS_45to49$risk),
                                          as.numeric(MSCA_PRS_50to54$risk), as.numeric(MSCA_PRS_55to59$risk), as.numeric(MSCA_PRS_60to64$risk), as.numeric(MSCA_PRS_65to69$risk),
                                          as.numeric(MSCA_PRS_70to74$risk)))

colnames(MSCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

MSCA_PRS_Average_Data4$risk_wtd_avg=(MSCA_PRS_Average_Data4$risk_30to34*0.108040659 + MSCA_PRS_Average_Data4$risk_35to39*0.104987889 + MSCA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                      MSCA_PRS_Average_Data4$risk_45to49*0.112278705 + MSCA_PRS_Average_Data4$risk_50to54*0.121739409 + MSCA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                      MSCA_PRS_Average_Data4$risk_60to64*0.125715739 + MSCA_PRS_Average_Data4$risk_65to69*0.109659337 + MSCA_PRS_Average_Data4$risk_70to74*0.086579466)

MSCA_PRS_AR_InfSS=MSCA_PRS_Average_Data4$risk_wtd_avg

save(MSCA_PRS_AR_CurrentSS,MSCA_PRS_AR_2xSS,MSCA_PRS_AR_4xSS,MSCA_PRS_AR_InfSS, file='MSCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
MSCA_PRS_AR_30to34=as.data.frame(cbind(MSCA_PRS_Average_Data1$risk_30to34, MSCA_PRS_Average_Data2$risk_30to34,
                                      MSCA_PRS_Average_Data3$risk_30to34, MSCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(MSCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

MSCA_PRS_AR_40to44=as.data.frame(cbind(MSCA_PRS_Average_Data1$risk_40to44, MSCA_PRS_Average_Data2$risk_40to44,
                                      MSCA_PRS_Average_Data3$risk_40to44, MSCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(MSCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

MSCA_PRS_AR_50to54=as.data.frame(cbind(MSCA_PRS_Average_Data1$risk_50to54, MSCA_PRS_Average_Data2$risk_50to54,
                                      MSCA_PRS_Average_Data3$risk_50to54, MSCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(MSCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

MSCA_PRS_AR_60to64=as.data.frame(cbind(MSCA_PRS_Average_Data1$risk_60to64, MSCA_PRS_Average_Data2$risk_60to64,
                                      MSCA_PRS_Average_Data3$risk_60to64, MSCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(MSCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(MSCA_PRS_AR_30to34, MSCA_PRS_AR_40to44, MSCA_PRS_AR_50to54, MSCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/MSCA_PRS_AR_AgeStrat.Rda')
      ############################
      # Ovarian cancer
      ############################
var1=0.0721132
var2=0.101750951
var3=0.146720099
var4=0.240963042
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
ovarian_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
OCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
OCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
OCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
OCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
OCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
OCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
OCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
OCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
OCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

OCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(OCA_PRS_30to34$risk), as.numeric(OCA_PRS_35to39$risk), as.numeric(OCA_PRS_40to44$risk), as.numeric(OCA_PRS_45to49$risk),
                                           as.numeric(OCA_PRS_50to54$risk), as.numeric(OCA_PRS_55to59$risk), as.numeric(OCA_PRS_60to64$risk), as.numeric(OCA_PRS_65to69$risk),
                                           as.numeric(OCA_PRS_70to74$risk)))

colnames(OCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

OCA_PRS_Average_Data1$risk_wtd_avg=(OCA_PRS_Average_Data1$risk_30to34*0.105319686 + OCA_PRS_Average_Data1$risk_35to39*0.102730061 + OCA_PRS_Average_Data1$risk_40to44*0.095906423 + 
                                       OCA_PRS_Average_Data1$risk_45to49*0.110520685 + OCA_PRS_Average_Data1$risk_50to54*0.120983294 + OCA_PRS_Average_Data1$risk_55to59*0.133676551 +
                                       OCA_PRS_Average_Data1$risk_60to64*0.12774873 + OCA_PRS_Average_Data1$risk_65to69*0.112737815 + OCA_PRS_Average_Data1$risk_70to74*0.090376754)

OCA_PRS_AR_CurrentSS=OCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
ovarian_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
OCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
OCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
OCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
OCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
OCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
OCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
OCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
OCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
OCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

OCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(OCA_PRS_30to34$risk), as.numeric(OCA_PRS_35to39$risk), as.numeric(OCA_PRS_40to44$risk), as.numeric(OCA_PRS_45to49$risk),
                                           as.numeric(OCA_PRS_50to54$risk), as.numeric(OCA_PRS_55to59$risk), as.numeric(OCA_PRS_60to64$risk), as.numeric(OCA_PRS_65to69$risk),
                                           as.numeric(OCA_PRS_70to74$risk)))

colnames(OCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

OCA_PRS_Average_Data2$risk_wtd_avg=(OCA_PRS_Average_Data2$risk_30to34*0.105319686 + OCA_PRS_Average_Data2$risk_35to39*0.102730061 + OCA_PRS_Average_Data2$risk_40to44*0.095906423 + 
                                       OCA_PRS_Average_Data2$risk_45to49*0.110520685 + OCA_PRS_Average_Data2$risk_50to54*0.120983294 + OCA_PRS_Average_Data2$risk_55to59*0.133676551 +
                                       OCA_PRS_Average_Data2$risk_60to64*0.12774873 + OCA_PRS_Average_Data2$risk_65to69*0.112737815 + OCA_PRS_Average_Data2$risk_70to74*0.090376754)

OCA_PRS_AR_2xSS=OCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
ovarian_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
OCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
OCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
OCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
OCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
OCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
OCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
OCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
OCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
OCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

OCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(OCA_PRS_30to34$risk), as.numeric(OCA_PRS_35to39$risk), as.numeric(OCA_PRS_40to44$risk), as.numeric(OCA_PRS_45to49$risk),
                                           as.numeric(OCA_PRS_50to54$risk), as.numeric(OCA_PRS_55to59$risk), as.numeric(OCA_PRS_60to64$risk), as.numeric(OCA_PRS_65to69$risk),
                                           as.numeric(OCA_PRS_70to74$risk)))

colnames(OCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

OCA_PRS_Average_Data3$risk_wtd_avg=(OCA_PRS_Average_Data3$risk_30to34*0.105319686 + OCA_PRS_Average_Data3$risk_35to39*0.102730061 + OCA_PRS_Average_Data3$risk_40to44*0.095906423 + 
                                       OCA_PRS_Average_Data3$risk_45to49*0.110520685 + OCA_PRS_Average_Data3$risk_50to54*0.120983294 + OCA_PRS_Average_Data3$risk_55to59*0.133676551 +
                                       OCA_PRS_Average_Data3$risk_60to64*0.12774873 + OCA_PRS_Average_Data3$risk_65to69*0.112737815 + OCA_PRS_Average_Data3$risk_70to74*0.090376754)

OCA_PRS_AR_4xSS=OCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
ovarian_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
OCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
OCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
OCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
OCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                                      model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
OCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
OCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
OCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
OCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
OCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = ovarian_ref,
                           model.disease.incidence.rates = ovarian_inc, model.competing.incidence.rates = mort_no_ovarian, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = ovarian_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

OCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(OCA_PRS_30to34$risk), as.numeric(OCA_PRS_35to39$risk), as.numeric(OCA_PRS_40to44$risk), as.numeric(OCA_PRS_45to49$risk),
                                           as.numeric(OCA_PRS_50to54$risk), as.numeric(OCA_PRS_55to59$risk), as.numeric(OCA_PRS_60to64$risk), as.numeric(OCA_PRS_65to69$risk),
                                           as.numeric(OCA_PRS_70to74$risk)))

colnames(OCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

OCA_PRS_Average_Data4$risk_wtd_avg=(OCA_PRS_Average_Data4$risk_30to34*0.105319686 + OCA_PRS_Average_Data4$risk_35to39*0.102730061 + OCA_PRS_Average_Data4$risk_40to44*0.095906423 + 
                                       OCA_PRS_Average_Data4$risk_45to49*0.110520685 + OCA_PRS_Average_Data4$risk_50to54*0.120983294 + OCA_PRS_Average_Data4$risk_55to59*0.133676551 +
                                       OCA_PRS_Average_Data4$risk_60to64*0.12774873 + OCA_PRS_Average_Data4$risk_65to69*0.112737815 + OCA_PRS_Average_Data4$risk_70to74*0.090376754)

OCA_PRS_AR_InfSS=OCA_PRS_Average_Data4$risk_wtd_avg

save(OCA_PRS_AR_CurrentSS,OCA_PRS_AR_2xSS,OCA_PRS_AR_4xSS,OCA_PRS_AR_InfSS, file='OCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
OCA_PRS_AR_30to34=as.data.frame(cbind(OCA_PRS_Average_Data1$risk_30to34, OCA_PRS_Average_Data2$risk_30to34,
                                      OCA_PRS_Average_Data3$risk_30to34, OCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(OCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

OCA_PRS_AR_40to44=as.data.frame(cbind(OCA_PRS_Average_Data1$risk_40to44, OCA_PRS_Average_Data2$risk_40to44,
                                      OCA_PRS_Average_Data3$risk_40to44, OCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(OCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

OCA_PRS_AR_50to54=as.data.frame(cbind(OCA_PRS_Average_Data1$risk_50to54, OCA_PRS_Average_Data2$risk_50to54,
                                      OCA_PRS_Average_Data3$risk_50to54, OCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(OCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

OCA_PRS_AR_60to64=as.data.frame(cbind(OCA_PRS_Average_Data1$risk_60to64, OCA_PRS_Average_Data2$risk_60to64,
                                      OCA_PRS_Average_Data3$risk_60to64, OCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(OCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(OCA_PRS_AR_30to34, OCA_PRS_AR_40to44, OCA_PRS_AR_50to54, OCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/OCA_PRS_AR_AgeStrat.Rda')

      ############################
      # Pancreatic cancer
      ############################
var1=0.132623955
var2=0.172315576
var3=0.238292227
var4=0.596304467
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var1))
var=rep(0, 100000)
pancreatic_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PACA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PACA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PACA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PACA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PACA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PACA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PACA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PACA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PACA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PACA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(PACA_PRS_30to34$risk), as.numeric(PACA_PRS_35to39$risk), as.numeric(PACA_PRS_40to44$risk), as.numeric(PACA_PRS_45to49$risk),
                                           as.numeric(PACA_PRS_50to54$risk), as.numeric(PACA_PRS_55to59$risk), as.numeric(PACA_PRS_60to64$risk), as.numeric(PACA_PRS_65to69$risk),
                                           as.numeric(PACA_PRS_70to74$risk)))

colnames(PACA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PACA_PRS_Average_Data1$risk_wtd_avg=(PACA_PRS_Average_Data1$risk_30to34*0.108040659 + PACA_PRS_Average_Data1$risk_35to39*0.104987889 + PACA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                       PACA_PRS_Average_Data1$risk_45to49*0.112278705 + PACA_PRS_Average_Data1$risk_50to54*0.121739409 + PACA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                       PACA_PRS_Average_Data1$risk_60to64*0.125715739 + PACA_PRS_Average_Data1$risk_65to69*0.109659337 + PACA_PRS_Average_Data1$risk_70to74*0.086579466)

PACA_PRS_AR_CurrentSS=PACA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var2))
var=rep(0, 100000)
pancreatic_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PACA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PACA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PACA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PACA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PACA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PACA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PACA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PACA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PACA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PACA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(PACA_PRS_30to34$risk), as.numeric(PACA_PRS_35to39$risk), as.numeric(PACA_PRS_40to44$risk), as.numeric(PACA_PRS_45to49$risk),
                                           as.numeric(PACA_PRS_50to54$risk), as.numeric(PACA_PRS_55to59$risk), as.numeric(PACA_PRS_60to64$risk), as.numeric(PACA_PRS_65to69$risk),
                                           as.numeric(PACA_PRS_70to74$risk)))

colnames(PACA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PACA_PRS_Average_Data2$risk_wtd_avg=(PACA_PRS_Average_Data2$risk_30to34*0.108040659 + PACA_PRS_Average_Data2$risk_35to39*0.104987889 + PACA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                       PACA_PRS_Average_Data2$risk_45to49*0.112278705 + PACA_PRS_Average_Data2$risk_50to54*0.121739409 + PACA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                       PACA_PRS_Average_Data2$risk_60to64*0.125715739 + PACA_PRS_Average_Data2$risk_65to69*0.109659337 + PACA_PRS_Average_Data2$risk_70to74*0.086579466)

PACA_PRS_AR_2xSS=PACA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var3))
var=rep(0, 100000)
pancreatic_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PACA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PACA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PACA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PACA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PACA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PACA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PACA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PACA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PACA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PACA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(PACA_PRS_30to34$risk), as.numeric(PACA_PRS_35to39$risk), as.numeric(PACA_PRS_40to44$risk), as.numeric(PACA_PRS_45to49$risk),
                                           as.numeric(PACA_PRS_50to54$risk), as.numeric(PACA_PRS_55to59$risk), as.numeric(PACA_PRS_60to64$risk), as.numeric(PACA_PRS_65to69$risk),
                                           as.numeric(PACA_PRS_70to74$risk)))

colnames(PACA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PACA_PRS_Average_Data3$risk_wtd_avg=(PACA_PRS_Average_Data3$risk_30to34*0.108040659 + PACA_PRS_Average_Data3$risk_35to39*0.104987889 + PACA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                       PACA_PRS_Average_Data3$risk_45to49*0.112278705 + PACA_PRS_Average_Data3$risk_50to54*0.121739409 + PACA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                       PACA_PRS_Average_Data3$risk_60to64*0.125715739 + PACA_PRS_Average_Data3$risk_65to69*0.109659337 + PACA_PRS_Average_Data3$risk_70to74*0.086579466)

PACA_PRS_AR_4xSS=PACA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(var4))
var=rep(0, 100000)
pancreatic_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PACA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PACA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PACA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PACA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                                      model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PACA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PACA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PACA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PACA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PACA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = pancreatic_ref,
                           model.disease.incidence.rates = pancreatic_inc, model.competing.incidence.rates = mort_no_pancreatic, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = pancreatic_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PACA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(PACA_PRS_30to34$risk), as.numeric(PACA_PRS_35to39$risk), as.numeric(PACA_PRS_40to44$risk), as.numeric(PACA_PRS_45to49$risk),
                                           as.numeric(PACA_PRS_50to54$risk), as.numeric(PACA_PRS_55to59$risk), as.numeric(PACA_PRS_60to64$risk), as.numeric(PACA_PRS_65to69$risk),
                                           as.numeric(PACA_PRS_70to74$risk)))

colnames(PACA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PACA_PRS_Average_Data4$risk_wtd_avg=(PACA_PRS_Average_Data4$risk_30to34*0.108040659 + PACA_PRS_Average_Data4$risk_35to39*0.104987889 + PACA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                       PACA_PRS_Average_Data4$risk_45to49*0.112278705 + PACA_PRS_Average_Data4$risk_50to54*0.121739409 + PACA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                       PACA_PRS_Average_Data4$risk_60to64*0.125715739 + PACA_PRS_Average_Data4$risk_65to69*0.109659337 + PACA_PRS_Average_Data4$risk_70to74*0.086579466)

PACA_PRS_AR_InfSS=PACA_PRS_Average_Data4$risk_wtd_avg

save(PACA_PRS_AR_CurrentSS,PACA_PRS_AR_2xSS,PACA_PRS_AR_4xSS,PACA_PRS_AR_InfSS, file='PACA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
PACA_PRS_AR_30to34=as.data.frame(cbind(PACA_PRS_Average_Data1$risk_30to34, PACA_PRS_Average_Data2$risk_30to34,
                                      PACA_PRS_Average_Data3$risk_30to34, PACA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(PACA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PACA_PRS_AR_40to44=as.data.frame(cbind(PACA_PRS_Average_Data1$risk_40to44, PACA_PRS_Average_Data2$risk_40to44,
                                      PACA_PRS_Average_Data3$risk_40to44, PACA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(PACA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PACA_PRS_AR_50to54=as.data.frame(cbind(PACA_PRS_Average_Data1$risk_50to54, PACA_PRS_Average_Data2$risk_50to54,
                                      PACA_PRS_Average_Data3$risk_50to54, PACA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(PACA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PACA_PRS_AR_60to64=as.data.frame(cbind(PACA_PRS_Average_Data1$risk_60to64, PACA_PRS_Average_Data2$risk_60to64,
                                      PACA_PRS_Average_Data3$risk_60to64, PACA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(PACA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(PACA_PRS_AR_30to34, PACA_PRS_AR_40to44, PACA_PRS_AR_50to54, PACA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/PACA_PRS_AR_AgeStrat.Rda')

      ############################
      # Prostate cancer
      ############################
var1=0.367091679
var2=0.454411072
var3=0.57244161
var4=0.771651228
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.275318202))
var=rep(0, 100000)
prostate_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(PCA_PRS_30to34$risk), as.numeric(PCA_PRS_35to39$risk), as.numeric(PCA_PRS_40to44$risk), as.numeric(PCA_PRS_45to49$risk),
                                           as.numeric(PCA_PRS_50to54$risk), as.numeric(PCA_PRS_55to59$risk), as.numeric(PCA_PRS_60to64$risk), as.numeric(PCA_PRS_65to69$risk),
                                           as.numeric(PCA_PRS_70to74$risk)))

colnames(PCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PCA_PRS_Average_Data1$risk_wtd_avg=(PCA_PRS_Average_Data1$risk_30to34*0.11083811 + PCA_PRS_Average_Data1$risk_35to39*0.107309178 + PCA_PRS_Average_Data1$risk_40to44*0.099664279 + 
                                       PCA_PRS_Average_Data1$risk_45to49*0.114086138 + PCA_PRS_Average_Data1$risk_50to54*0.122516776 + PCA_PRS_Average_Data1$risk_55to59*0.13279013 +
                                       PCA_PRS_Average_Data1$risk_60to64*0.123625607 + PCA_PRS_Average_Data1$risk_65to69*0.106494333 + PCA_PRS_Average_Data1$risk_70to74*0.082675449)

PCA_PRS_AR_CurrentSS=PCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.38622653))
var=rep(0, 100000)
prostate_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(PCA_PRS_30to34$risk), as.numeric(PCA_PRS_35to39$risk), as.numeric(PCA_PRS_40to44$risk), as.numeric(PCA_PRS_45to49$risk),
                                           as.numeric(PCA_PRS_50to54$risk), as.numeric(PCA_PRS_55to59$risk), as.numeric(PCA_PRS_60to64$risk), as.numeric(PCA_PRS_65to69$risk),
                                           as.numeric(PCA_PRS_70to74$risk)))

colnames(PCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PCA_PRS_Average_Data2$risk_wtd_avg=(PCA_PRS_Average_Data2$risk_30to34*0.11083811 + PCA_PRS_Average_Data2$risk_35to39*0.107309178 + PCA_PRS_Average_Data2$risk_40to44*0.099664279 + 
                                      PCA_PRS_Average_Data2$risk_45to49*0.114086138 + PCA_PRS_Average_Data2$risk_50to54*0.122516776 + PCA_PRS_Average_Data2$risk_55to59*0.13279013 +
                                      PCA_PRS_Average_Data2$risk_60to64*0.123625607 + PCA_PRS_Average_Data2$risk_65to69*0.106494333 + PCA_PRS_Average_Data2$risk_70to74*0.082675449)

PCA_PRS_AR_2xSS=PCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.51383947))
var=rep(0, 100000)
prostate_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(PCA_PRS_30to34$risk), as.numeric(PCA_PRS_35to39$risk), as.numeric(PCA_PRS_40to44$risk), as.numeric(PCA_PRS_45to49$risk),
                                           as.numeric(PCA_PRS_50to54$risk), as.numeric(PCA_PRS_55to59$risk), as.numeric(PCA_PRS_60to64$risk), as.numeric(PCA_PRS_65to69$risk),
                                           as.numeric(PCA_PRS_70to74$risk)))

colnames(PCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PCA_PRS_Average_Data3$risk_wtd_avg=(PCA_PRS_Average_Data3$risk_30to34*0.11083811 + PCA_PRS_Average_Data3$risk_35to39*0.107309178 + PCA_PRS_Average_Data3$risk_40to44*0.099664279 + 
                                      PCA_PRS_Average_Data3$risk_45to49*0.114086138 + PCA_PRS_Average_Data3$risk_50to54*0.122516776 + PCA_PRS_Average_Data3$risk_55to59*0.13279013 +
                                      PCA_PRS_Average_Data3$risk_60to64*0.123625607 + PCA_PRS_Average_Data3$risk_65to69*0.106494333 + PCA_PRS_Average_Data3$risk_70to74*0.082675449)

PCA_PRS_AR_4xSS=PCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.68099864))
var=rep(0, 100000)
prostate_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
PCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
PCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
PCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
PCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                                      model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
PCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
PCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
PCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
PCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
PCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = prostate_ref,
                           model.disease.incidence.rates = prostate_inc, model.competing.incidence.rates = mort_no_prostate, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = prostate_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

PCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(PCA_PRS_30to34$risk), as.numeric(PCA_PRS_35to39$risk), as.numeric(PCA_PRS_40to44$risk), as.numeric(PCA_PRS_45to49$risk),
                                           as.numeric(PCA_PRS_50to54$risk), as.numeric(PCA_PRS_55to59$risk), as.numeric(PCA_PRS_60to64$risk), as.numeric(PCA_PRS_65to69$risk),
                                           as.numeric(PCA_PRS_70to74$risk)))

colnames(PCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

PCA_PRS_Average_Data4$risk_wtd_avg=(PCA_PRS_Average_Data4$risk_30to34*0.11083811 + PCA_PRS_Average_Data4$risk_35to39*0.107309178 + PCA_PRS_Average_Data4$risk_40to44*0.099664279 + 
                                      PCA_PRS_Average_Data4$risk_45to49*0.114086138 + PCA_PRS_Average_Data4$risk_50to54*0.122516776 + PCA_PRS_Average_Data4$risk_55to59*0.13279013 +
                                      PCA_PRS_Average_Data4$risk_60to64*0.123625607 + PCA_PRS_Average_Data4$risk_65to69*0.106494333 + PCA_PRS_Average_Data4$risk_70to74*0.082675449)

PCA_PRS_AR_InfSS=PCA_PRS_Average_Data4$risk_wtd_avg

save(PCA_PRS_AR_CurrentSS,PCA_PRS_AR_2xSS,PCA_PRS_AR_4xSS,PCA_PRS_AR_InfSS, file='PCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
PCA_PRS_AR_30to34=as.data.frame(cbind(PCA_PRS_Average_Data1$risk_30to34, PCA_PRS_Average_Data2$risk_30to34,
                                      PCA_PRS_Average_Data3$risk_30to34, PCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(PCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PCA_PRS_AR_40to44=as.data.frame(cbind(PCA_PRS_Average_Data1$risk_40to44, PCA_PRS_Average_Data2$risk_40to44,
                                      PCA_PRS_Average_Data3$risk_40to44, PCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(PCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PCA_PRS_AR_50to54=as.data.frame(cbind(PCA_PRS_Average_Data1$risk_50to54, PCA_PRS_Average_Data2$risk_50to54,
                                      PCA_PRS_Average_Data3$risk_50to54, PCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(PCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

PCA_PRS_AR_60to64=as.data.frame(cbind(PCA_PRS_Average_Data1$risk_60to64, PCA_PRS_Average_Data2$risk_60to64,
                                      PCA_PRS_Average_Data3$risk_60to64, PCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(PCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(PCA_PRS_AR_30to34, PCA_PRS_AR_40to44, PCA_PRS_AR_50to54, PCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/PCA_PRS_AR_AgeStrat.Rda')
      ############################
      # Renal cancer
      ############################
var1=0.096621914
var2=0.132677739
var3=0.206765492
var4=0.565688661
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.075253584))
var=rep(0, 100000)
renal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
RCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
RCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
RCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
RCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
RCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
RCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
RCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
RCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
RCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

RCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(RCA_PRS_30to34$risk), as.numeric(RCA_PRS_35to39$risk), as.numeric(RCA_PRS_40to44$risk), as.numeric(RCA_PRS_45to49$risk),
                                           as.numeric(RCA_PRS_50to54$risk), as.numeric(RCA_PRS_55to59$risk), as.numeric(RCA_PRS_60to64$risk), as.numeric(RCA_PRS_65to69$risk),
                                           as.numeric(RCA_PRS_70to74$risk)))

colnames(RCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

RCA_PRS_Average_Data1$risk_wtd_avg=(RCA_PRS_Average_Data1$risk_30to34*0.108040659 + RCA_PRS_Average_Data1$risk_35to39*0.104987889 + RCA_PRS_Average_Data1$risk_40to44*0.097759312 + 
                                       RCA_PRS_Average_Data1$risk_45to49*0.112278705 + RCA_PRS_Average_Data1$risk_50to54*0.121739409 + RCA_PRS_Average_Data1$risk_55to59*0.133239483 +
                                       RCA_PRS_Average_Data1$risk_60to64*0.125715739 + RCA_PRS_Average_Data1$risk_65to69*0.109659337 + RCA_PRS_Average_Data1$risk_70to74*0.086579466)

RCA_PRS_AR_CurrentSS=RCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.110832963))
var=rep(0, 100000)
renal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
RCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
RCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
RCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
RCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
RCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
RCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
RCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
RCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
RCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

RCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(RCA_PRS_30to34$risk), as.numeric(RCA_PRS_35to39$risk), as.numeric(RCA_PRS_40to44$risk), as.numeric(RCA_PRS_45to49$risk),
                                           as.numeric(RCA_PRS_50to54$risk), as.numeric(RCA_PRS_55to59$risk), as.numeric(RCA_PRS_60to64$risk), as.numeric(RCA_PRS_65to69$risk),
                                           as.numeric(RCA_PRS_70to74$risk)))

colnames(RCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

RCA_PRS_Average_Data2$risk_wtd_avg=(RCA_PRS_Average_Data2$risk_30to34*0.108040659 + RCA_PRS_Average_Data2$risk_35to39*0.104987889 + RCA_PRS_Average_Data2$risk_40to44*0.097759312 + 
                                       RCA_PRS_Average_Data2$risk_45to49*0.112278705 + RCA_PRS_Average_Data2$risk_50to54*0.121739409 + RCA_PRS_Average_Data2$risk_55to59*0.133239483 +
                                       RCA_PRS_Average_Data2$risk_60to64*0.125715739 + RCA_PRS_Average_Data2$risk_65to69*0.109659337 + RCA_PRS_Average_Data2$risk_70to74*0.086579466)

RCA_PRS_AR_2xSS=RCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.191012429))
var=rep(0, 100000)
renal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
RCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
RCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
RCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
RCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
RCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
RCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
RCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
RCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
RCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

RCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(RCA_PRS_30to34$risk), as.numeric(RCA_PRS_35to39$risk), as.numeric(RCA_PRS_40to44$risk), as.numeric(RCA_PRS_45to49$risk),
                                           as.numeric(RCA_PRS_50to54$risk), as.numeric(RCA_PRS_55to59$risk), as.numeric(RCA_PRS_60to64$risk), as.numeric(RCA_PRS_65to69$risk),
                                           as.numeric(RCA_PRS_70to74$risk)))

colnames(RCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

RCA_PRS_Average_Data3$risk_wtd_avg=(RCA_PRS_Average_Data3$risk_30to34*0.108040659 + RCA_PRS_Average_Data3$risk_35to39*0.104987889 + RCA_PRS_Average_Data3$risk_40to44*0.097759312 + 
                                       RCA_PRS_Average_Data3$risk_45to49*0.112278705 + RCA_PRS_Average_Data3$risk_50to54*0.121739409 + RCA_PRS_Average_Data3$risk_55to59*0.133239483 +
                                       RCA_PRS_Average_Data3$risk_60to64*0.125715739 + RCA_PRS_Average_Data3$risk_65to69*0.109659337 + RCA_PRS_Average_Data3$risk_70to74*0.086579466)

RCA_PRS_AR_4xSS=RCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(0.545999396))
var=rep(0, 100000)
renal_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
RCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=32,
                                      apply.age.interval.length = 43, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
RCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=37,
                                      apply.age.interval.length = 38, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
RCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=42,
                                      apply.age.interval.length = 33, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
RCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                                      model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=47,
                                      apply.age.interval.length = 28, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
RCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
RCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
RCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
RCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
RCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = renal_ref,
                           model.disease.incidence.rates = renal_inc, model.competing.incidence.rates = mort_no_renal, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = renal_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

RCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(RCA_PRS_30to34$risk), as.numeric(RCA_PRS_35to39$risk), as.numeric(RCA_PRS_40to44$risk), as.numeric(RCA_PRS_45to49$risk),
                                           as.numeric(RCA_PRS_50to54$risk), as.numeric(RCA_PRS_55to59$risk), as.numeric(RCA_PRS_60to64$risk), as.numeric(RCA_PRS_65to69$risk),
                                           as.numeric(RCA_PRS_70to74$risk)))

colnames(RCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

RCA_PRS_Average_Data4$risk_wtd_avg=(RCA_PRS_Average_Data4$risk_30to34*0.108040659 + RCA_PRS_Average_Data4$risk_35to39*0.104987889 + RCA_PRS_Average_Data4$risk_40to44*0.097759312 + 
                                       RCA_PRS_Average_Data4$risk_45to49*0.112278705 + RCA_PRS_Average_Data4$risk_50to54*0.121739409 + RCA_PRS_Average_Data4$risk_55to59*0.133239483 +
                                       RCA_PRS_Average_Data4$risk_60to64*0.125715739 + RCA_PRS_Average_Data4$risk_65to69*0.109659337 + RCA_PRS_Average_Data4$risk_70to74*0.086579466)

RCA_PRS_AR_InfSS=RCA_PRS_Average_Data4$risk_wtd_avg

save(RCA_PRS_AR_CurrentSS,RCA_PRS_AR_2xSS,RCA_PRS_AR_4xSS,RCA_PRS_AR_InfSS, file='RCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
RCA_PRS_AR_30to34=as.data.frame(cbind(RCA_PRS_Average_Data1$risk_30to34, RCA_PRS_Average_Data2$risk_30to34,
                                      RCA_PRS_Average_Data3$risk_30to34, RCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(RCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

RCA_PRS_AR_40to44=as.data.frame(cbind(RCA_PRS_Average_Data1$risk_40to44, RCA_PRS_Average_Data2$risk_40to44,
                                      RCA_PRS_Average_Data3$risk_40to44, RCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(RCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

RCA_PRS_AR_50to54=as.data.frame(cbind(RCA_PRS_Average_Data1$risk_50to54, RCA_PRS_Average_Data2$risk_50to54,
                                      RCA_PRS_Average_Data3$risk_50to54, RCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(RCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

RCA_PRS_AR_60to64=as.data.frame(cbind(RCA_PRS_Average_Data1$risk_60to64, RCA_PRS_Average_Data2$risk_60to64,
                                      RCA_PRS_Average_Data3$risk_60to64, RCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(RCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(RCA_PRS_AR_30to34, RCA_PRS_AR_40to44, RCA_PRS_AR_50to54, RCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/RCA_PRS_AR_AgeStrat.Rda')

      ############################
      # Testicular cancer
      ############################
var1=1.102214479
var2=1.38300417
var3=1.710356629
var4=2.803422789
##### CURRENT SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(1.013321552))
var=rep(0, 100000)
testicular_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
TCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
TCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
TCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
TCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
TCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
TCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
TCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
TCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
TCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

TCA_PRS_Average_Data1=as.data.frame(cbind(as.numeric(TCA_PRS_30to34$risk), as.numeric(TCA_PRS_35to39$risk), as.numeric(TCA_PRS_40to44$risk), as.numeric(TCA_PRS_45to49$risk),
                                          as.numeric(TCA_PRS_50to54$risk), as.numeric(TCA_PRS_55to59$risk), as.numeric(TCA_PRS_60to64$risk), as.numeric(TCA_PRS_65to69$risk),
                                          as.numeric(TCA_PRS_70to74$risk)))
colnames(TCA_PRS_Average_Data1)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

TCA_PRS_Average_Data1$risk_wtd_avg=(TCA_PRS_Average_Data1$risk_30to34*0.11083811 + TCA_PRS_Average_Data1$risk_35to39*0.107309178 + TCA_PRS_Average_Data1$risk_40to44*0.099664279 + 
                                    TCA_PRS_Average_Data1$risk_45to49*0.114086138 + TCA_PRS_Average_Data1$risk_50to54*0.122516776 + TCA_PRS_Average_Data1$risk_55to59*0.13279013 +
                                    TCA_PRS_Average_Data1$risk_60to64*0.123625607 + TCA_PRS_Average_Data1$risk_65to69*0.106494333 + TCA_PRS_Average_Data1$risk_70to74*0.082675449)

TCA_PRS_AR_CurrentSS=TCA_PRS_Average_Data1$risk_wtd_avg

##### 2X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(1.447451055))
var=rep(0, 100000)
testicular_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
TCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
TCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
TCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
TCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
TCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
TCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
TCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
TCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
TCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

TCA_PRS_Average_Data2=as.data.frame(cbind(as.numeric(TCA_PRS_30to34$risk), as.numeric(TCA_PRS_35to39$risk), as.numeric(TCA_PRS_40to44$risk), as.numeric(TCA_PRS_45to49$risk),
                                          as.numeric(TCA_PRS_50to54$risk), as.numeric(TCA_PRS_55to59$risk), as.numeric(TCA_PRS_60to64$risk), as.numeric(TCA_PRS_65to69$risk),
                                          as.numeric(TCA_PRS_70to74$risk)))
colnames(TCA_PRS_Average_Data2)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

TCA_PRS_Average_Data2$risk_wtd_avg=(TCA_PRS_Average_Data2$risk_30to34*0.11083811 + TCA_PRS_Average_Data2$risk_35to39*0.107309178 + TCA_PRS_Average_Data2$risk_40to44*0.099664279 + 
                                    TCA_PRS_Average_Data2$risk_45to49*0.114086138 + TCA_PRS_Average_Data2$risk_50to54*0.122516776 + TCA_PRS_Average_Data2$risk_55to59*0.13279013 +
                                    TCA_PRS_Average_Data2$risk_60to64*0.123625607 + TCA_PRS_Average_Data2$risk_65to69*0.106494333 + TCA_PRS_Average_Data2$risk_70to74*0.082675449)

TCA_PRS_AR_2xSS=TCA_PRS_Average_Data2$risk_wtd_avg

##### 4X SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(1.957205694))
var=rep(0, 100000)
testicular_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
TCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
TCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
TCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
TCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
TCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
TCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
TCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
TCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
TCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

TCA_PRS_Average_Data3=as.data.frame(cbind(as.numeric(TCA_PRS_30to34$risk), as.numeric(TCA_PRS_35to39$risk), as.numeric(TCA_PRS_40to44$risk), as.numeric(TCA_PRS_45to49$risk),
                                          as.numeric(TCA_PRS_50to54$risk), as.numeric(TCA_PRS_55to59$risk), as.numeric(TCA_PRS_60to64$risk), as.numeric(TCA_PRS_65to69$risk),
                                          as.numeric(TCA_PRS_70to74$risk)))
colnames(TCA_PRS_Average_Data3)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

TCA_PRS_Average_Data3$risk_wtd_avg=(TCA_PRS_Average_Data3$risk_30to34*0.11083811 + TCA_PRS_Average_Data3$risk_35to39*0.107309178 + TCA_PRS_Average_Data3$risk_40to44*0.099664279 + 
                                    TCA_PRS_Average_Data3$risk_45to49*0.114086138 + TCA_PRS_Average_Data3$risk_50to54*0.122516776 + TCA_PRS_Average_Data3$risk_55to59*0.13279013 +
                                    TCA_PRS_Average_Data3$risk_60to64*0.123625607 + TCA_PRS_Average_Data3$risk_65to69*0.106494333 + TCA_PRS_Average_Data3$risk_70to74*0.082675449)

TCA_PRS_AR_4xSS=TCA_PRS_Average_Data3$risk_wtd_avg

##### INFINITE SAMPLE SIZE #####
PRS=rnorm(n=100000, mean=1, sd=sqrt(2.649999863))
var=rep(0, 100000)
testicular_ref=as.data.frame(cbind(PRS, var))

# AGES 30-34 YEARS #
TCA_PRS_30to34 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=32,
                           apply.age.interval.length = 43, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 35-39 YEARS #
TCA_PRS_35to39 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=37,
                           apply.age.interval.length = 38, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 40-44 YEARS #
TCA_PRS_40to44 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=42,
                           apply.age.interval.length = 33, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 45-49 YEARS #
TCA_PRS_45to49 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=47,
                           apply.age.interval.length = 28, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 50-54 YEARS #
TCA_PRS_50to54 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=52,
                           apply.age.interval.length = 23, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 55-59 YEARS #
TCA_PRS_55to59 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=57,
                           apply.age.interval.length = 18, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 60-64 YEARS #
TCA_PRS_60to64 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=62,
                           apply.age.interval.length = 13, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 65-69 YEARS # 
TCA_PRS_65to69 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=67,
                           apply.age.interval.length = 8, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

# AGES 70-74 YEARS #
TCA_PRS_70to74 = computeAbsoluteRisk(model.formula = model.PRS, model.cov.info = cov_info_PRS, model.snp.info=NULL, model.log.RR =beta_PRS, model.ref.dataset = testicular_ref,
                           model.disease.incidence.rates = testicular_inc, model.competing.incidence.rates = mort_no_testicular, model.bin.fh.name=NULL, apply.age.start=72,
                           apply.age.interval.length = 3, apply.cov.profile = testicular_ref, apply.snp.profile = NULL, use.c.code=1, return.lp = T)

TCA_PRS_Average_Data4=as.data.frame(cbind(as.numeric(TCA_PRS_30to34$risk), as.numeric(TCA_PRS_35to39$risk), as.numeric(TCA_PRS_40to44$risk), as.numeric(TCA_PRS_45to49$risk),
                                          as.numeric(TCA_PRS_50to54$risk), as.numeric(TCA_PRS_55to59$risk), as.numeric(TCA_PRS_60to64$risk), as.numeric(TCA_PRS_65to69$risk),
                                          as.numeric(TCA_PRS_70to74$risk)))
colnames(TCA_PRS_Average_Data4)=c('risk_30to34','risk_35to39','risk_40to44','risk_45to49','risk_50to54','risk_55to59','risk_60to64','risk_65to69','risk_70to74')

TCA_PRS_Average_Data4$risk_wtd_avg=(TCA_PRS_Average_Data4$risk_30to34*0.11083811 + TCA_PRS_Average_Data4$risk_35to39*0.107309178 + TCA_PRS_Average_Data4$risk_40to44*0.099664279 + 
                                    TCA_PRS_Average_Data4$risk_45to49*0.114086138 + TCA_PRS_Average_Data4$risk_50to54*0.122516776 + TCA_PRS_Average_Data4$risk_55to59*0.13279013 +
                                    TCA_PRS_Average_Data4$risk_60to64*0.123625607 + TCA_PRS_Average_Data4$risk_65to69*0.106494333 + TCA_PRS_Average_Data4$risk_70to74*0.082675449)

TCA_PRS_AR_InfSS=TCA_PRS_Average_Data4$risk_wtd_avg

save(TCA_PRS_AR_CurrentSS,TCA_PRS_AR_2xSS,TCA_PRS_AR_4xSS,TCA_PRS_AR_InfSS, file='TCA_PRS_AR.Rda')

##### AGE-STRATIFIED PLOTTING DATA #####
TCA_PRS_AR_30to34=as.data.frame(cbind(TCA_PRS_Average_Data1$risk_30to34, TCA_PRS_Average_Data2$risk_30to34,
                                      TCA_PRS_Average_Data3$risk_30to34, TCA_PRS_Average_Data4$risk_30to34)*100)
                                colnames(TCA_PRS_AR_30to34)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

TCA_PRS_AR_40to44=as.data.frame(cbind(TCA_PRS_Average_Data1$risk_40to44, TCA_PRS_Average_Data2$risk_40to44,
                                      TCA_PRS_Average_Data3$risk_40to44, TCA_PRS_Average_Data4$risk_40to44)*100)
                                colnames(TCA_PRS_AR_40to44)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

TCA_PRS_AR_50to54=as.data.frame(cbind(TCA_PRS_Average_Data1$risk_50to54, TCA_PRS_Average_Data2$risk_50to54,
                                      TCA_PRS_Average_Data3$risk_50to54, TCA_PRS_Average_Data4$risk_50to54)*100)
                                colnames(TCA_PRS_AR_50to54)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

TCA_PRS_AR_60to64=as.data.frame(cbind(TCA_PRS_Average_Data1$risk_60to64, TCA_PRS_Average_Data2$risk_60to64,
                                      TCA_PRS_Average_Data3$risk_60to64, TCA_PRS_Average_Data4$risk_60to64)*100)
                                colnames(TCA_PRS_AR_60to64)=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS')

save(TCA_PRS_AR_30to34, TCA_PRS_AR_40to44, TCA_PRS_AR_50to54, TCA_PRS_AR_60to64, 
           file='/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files/TCA_PRS_AR_AgeStrat.Rda')




#########################################
############################
# Plotting Absolute Risk
############################
#########################################
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
theme_set(theme_cowplot())

AR_dens_plot<- function(dataset, site=NULL, x.lims=NULL, x.breaks=waiver(), y.lims=NULL, y.breaks=waiver(), with_key=F) {
  
if(with_key==T){

  ggplot(dataset, aes(x=AB_risk_pct, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) + ggtitle(site) +
  scale_x_continuous("\nResidual absolute risk (%)", limits=x.lims, breaks=x.breaks, expand = c(0, 0)) +
  scale_y_continuous("Density", limits=y.lims, breaks=y.breaks, expand = c(0, 0)) +
  theme(plot.title = element_text(size = 18, colour = "black", vjust=3, hjust=0.5), axis.text.x=element_text(size=14),axis.text.y=element_text(size=14), 
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), legend.direction="vertical", 
        legend.position=c(.15,.8), legend.title=element_blank(), legend.text=element_text(size=14), 
        legend.key.size=unit(0.7,'cm'), plot.margin=unit(c(0.5,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size"))

}

  else {
    
  ggplot(dataset, aes(x=AB_risk_pct, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) + ggtitle(site) +
  scale_x_continuous("\nResidual absolute risk (%)", limits=x.lims, breaks=x.breaks, expand = c(0, 0)) +
  scale_y_continuous("Density", limits=y.lims, breaks=y.breaks, expand = c(0, 0)) +
  theme(plot.title = element_text(size = 18, colour = "black", vjust=3, hjust=0.5), axis.text.x=element_text(size=14),axis.text.y=element_text(size=14), 
        axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), plot.margin=unit(c(0.5,0.75,0.3,0.3),"cm"), legend.position = "none") +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) 
  }
}

AR_dens_plot2<- function(dataset, site=NULL, x.lims=NULL, x.breaks=waiver(), y.lims=NULL, y.breaks=waiver(), with_key=F) {
  
if(with_key==T){

  ggplot(dataset, aes(x=AB_risk_pct, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) + ggtitle(site) +
  scale_x_continuous(" ", limits=x.lims, breaks=x.breaks, expand=c(0.03,0)) +
  scale_y_continuous(" ", limits=y.lims, breaks=y.breaks, expand=c(0,0)) +
  theme(plot.title = element_text(face="plain", size = 26, colour = "black", vjust=3, hjust=0.5), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24), 
        legend.direction="vertical", legend.position=c(.4,.65), legend.title=element_blank(), legend.text=element_text(size=26), 
        legend.key.size=unit(0.8,'cm'), plot.margin=unit(c(0.5,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current"," Double"," Quadruple"," Infinite"))

}

  else {
    
  ggplot(dataset, aes(x=AB_risk_pct, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) + ggtitle(site) +
  scale_x_continuous(" ", limits=x.lims, breaks=x.breaks, expand=c(0.03,0)) +
  scale_y_continuous(" ", limits=y.lims, breaks=y.breaks, expand=c(0,0)) +
  theme(plot.title = element_text(face="plain", size = 26, colour = "black", vjust=3, hjust=0.5), axis.text.x=element_text(size=24), axis.text.y=element_text(size=24), 
        plot.margin=unit(c(0.5,0.75,0.3,0.3),"cm"), legend.position = "none") +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current"," Double"," Quadruple"," Infinite")) 
  }
}

# Creating datasets for plotting ####
load('CRCA_PRS_AR.Rda')
  CRCA_ref1=as.data.frame(CRCA_PRS_AR_CurrentSS)
      colnames(CRCA_ref1)='risk'
      CRCA_ref1$AB_risk_pct=CRCA_ref1$risk*100
      CRCA_ref1$size='current'
  CRCA_ref2=as.data.frame(CRCA_PRS_AR_2xSS)
      colnames(CRCA_ref2)='risk'
      CRCA_ref2$AB_risk_pct=CRCA_ref2$risk*100
      CRCA_ref2$size='2X'
  CRCA_ref3=as.data.frame(CRCA_PRS_AR_4xSS)
      colnames(CRCA_ref3)='risk'
      CRCA_ref3$AB_risk_pct=CRCA_ref3$risk*100
      CRCA_ref3$size='4X'
  CRCA_ref4=as.data.frame(CRCA_PRS_AR_InfSS)
      colnames(CRCA_ref4)='risk'
      CRCA_ref4$AB_risk_pct=CRCA_ref4$risk*100
      CRCA_ref4$size='Inf'
CRCA_ref=rbind(CRCA_ref1, CRCA_ref2, CRCA_ref3, CRCA_ref4)
CRCA_ref$size <- factor(CRCA_ref$size, levels=c('current','2X','4X','Inf'))

load('PACA_PRS_AR.Rda')
  PACA_ref1=as.data.frame(PACA_PRS_AR_CurrentSS)
      colnames(PACA_ref1)='risk'
      PACA_ref1$AB_risk_pct=PACA_ref1$risk*100
      PACA_ref1$size='current'
  PACA_ref2=as.data.frame(PACA_PRS_AR_2xSS)
      colnames(PACA_ref2)='risk'
      PACA_ref2$AB_risk_pct=PACA_ref2$risk*100
      PACA_ref2$size='2X'
  PACA_ref3=as.data.frame(PACA_PRS_AR_4xSS)
      colnames(PACA_ref3)='risk'
      PACA_ref3$AB_risk_pct=PACA_ref3$risk*100
      PACA_ref3$size='4X'
  PACA_ref4=as.data.frame(PACA_PRS_AR_InfSS)
      colnames(PACA_ref4)='risk'
      PACA_ref4$AB_risk_pct=PACA_ref4$risk*100
      PACA_ref4$size='Inf'
PACA_ref=rbind(PACA_ref1, PACA_ref2, PACA_ref3, PACA_ref4)
PACA_ref$size <- factor(PACA_ref$size, levels=c('current','2X','4X','Inf'))

load('UECA_PRS_AR.Rda')
  UECA_ref1=as.data.frame(UECA_PRS_AR_CurrentSS)
      colnames(UECA_ref1)='risk'
      UECA_ref1$AB_risk_pct=UECA_ref1$risk*100
      UECA_ref1$size='current'
  UECA_ref2=as.data.frame(UECA_PRS_AR_2xSS)
      colnames(UECA_ref2)='risk'
      UECA_ref2$AB_risk_pct=UECA_ref2$risk*100
      UECA_ref2$size='2X'
  UECA_ref3=as.data.frame(UECA_PRS_AR_4xSS)
      colnames(UECA_ref3)='risk'
      UECA_ref3$AB_risk_pct=UECA_ref3$risk*100
      UECA_ref3$size='4X'
  UECA_ref4=as.data.frame(UECA_PRS_AR_InfSS)
      colnames(UECA_ref4)='risk'
      UECA_ref4$AB_risk_pct=UECA_ref4$risk*100
      UECA_ref4$size='Inf'
UECA_ref=rbind(UECA_ref1, UECA_ref2, UECA_ref3, UECA_ref4)
UECA_ref$size <- factor(UECA_ref$size, levels=c('current','2X','4X','Inf'))

load('RCA_PRS_AR.Rda')
  RCA_ref1=as.data.frame(RCA_PRS_AR_CurrentSS)
      colnames(RCA_ref1)='risk'
      RCA_ref1$AB_risk_pct=RCA_ref1$risk*100
      RCA_ref1$size='current'
  RCA_ref2=as.data.frame(RCA_PRS_AR_2xSS)
      colnames(RCA_ref2)='risk'
      RCA_ref2$AB_risk_pct=RCA_ref2$risk*100
      RCA_ref2$size='2X'
  RCA_ref3=as.data.frame(RCA_PRS_AR_4xSS)
      colnames(RCA_ref3)='risk'
      RCA_ref3$AB_risk_pct=RCA_ref3$risk*100
      RCA_ref3$size='4X'
  RCA_ref4=as.data.frame(RCA_PRS_AR_InfSS)
      colnames(RCA_ref4)='risk'
      RCA_ref4$AB_risk_pct=RCA_ref4$risk*100
      RCA_ref4$size='Inf'
RCA_ref=rbind(RCA_ref1, RCA_ref2, RCA_ref3, RCA_ref4)
RCA_ref$size <- factor(RCA_ref$size, levels=c('current','2X','4X','Inf'))

load('CLL_PRS_AR.Rda')
  CLL_ref1=as.data.frame(CLL_PRS_AR_CurrentSS)
      colnames(CLL_ref1)='risk'
      CLL_ref1$AB_risk_pct=CLL_ref1$risk*100
      CLL_ref1$size='current'
  CLL_ref2=as.data.frame(CLL_PRS_AR_2xSS)
      colnames(CLL_ref2)='risk'
      CLL_ref2$AB_risk_pct=CLL_ref2$risk*100
      CLL_ref2$size='2X'
  CLL_ref3=as.data.frame(CLL_PRS_AR_4xSS)
      colnames(CLL_ref3)='risk'
      CLL_ref3$AB_risk_pct=CLL_ref3$risk*100
      CLL_ref3$size='4X'
  CLL_ref4=as.data.frame(CLL_PRS_AR_InfSS)
      colnames(CLL_ref4)='risk'
      CLL_ref4$AB_risk_pct=CLL_ref4$risk*100
      CLL_ref4$size='Inf'
CLL_ref=rbind(CLL_ref1, CLL_ref2, CLL_ref3, CLL_ref4)
CLL_ref$size <- factor(CLL_ref$size, levels=c('current','2X','4X','Inf'))

load('ECA_PRS_AR.Rda')
  ECA_ref1=as.data.frame(ECA_PRS_AR_CurrentSS)
      colnames(ECA_ref1)='risk'
      ECA_ref1$AB_risk_pct=ECA_ref1$risk*100
      ECA_ref1$size='current'
  ECA_ref2=as.data.frame(ECA_PRS_AR_2xSS)
      colnames(ECA_ref2)='risk'
      ECA_ref2$AB_risk_pct=ECA_ref2$risk*100
      ECA_ref2$size='2X'
  ECA_ref3=as.data.frame(ECA_PRS_AR_4xSS)
      colnames(ECA_ref3)='risk'
      ECA_ref3$AB_risk_pct=ECA_ref3$risk*100
      ECA_ref3$size='4X'
  ECA_ref4=as.data.frame(ECA_PRS_AR_InfSS)
      colnames(ECA_ref4)='risk'
      ECA_ref4$AB_risk_pct=ECA_ref4$risk*100
      ECA_ref4$size='Inf'
ECA_ref=rbind(ECA_ref1, ECA_ref2, ECA_ref3, ECA_ref4)
ECA_ref$size <- factor(ECA_ref$size, levels=c('current','2X','4X','Inf'))

load('TCA_PRS_AR.Rda')
  TCA_ref1=as.data.frame(TCA_PRS_AR_CurrentSS)
      colnames(TCA_ref1)='risk'
      TCA_ref1$AB_risk_pct=TCA_ref1$risk*100
      TCA_ref1$size='current'
  TCA_ref2=as.data.frame(TCA_PRS_AR_2xSS)
      colnames(TCA_ref2)='risk'
      TCA_ref2$AB_risk_pct=TCA_ref2$risk*100
      TCA_ref2$size='2X'
  TCA_ref3=as.data.frame(TCA_PRS_AR_4xSS)
      colnames(TCA_ref3)='risk'
      TCA_ref3$AB_risk_pct=TCA_ref3$risk*100
      TCA_ref3$size='4X'
  TCA_ref4=as.data.frame(TCA_PRS_AR_InfSS)
      colnames(TCA_ref4)='risk'
      TCA_ref4$AB_risk_pct=TCA_ref4$risk*100
      TCA_ref4$size='Inf'
TCA_ref=rbind(TCA_ref1, TCA_ref2, TCA_ref3, TCA_ref4)
TCA_ref$size <- factor(TCA_ref$size, levels=c('current','2X','4X','Inf'))

load('GCA_PRS_AR.Rda')
  GCA_ref1=as.data.frame(GCA_PRS_AR_CurrentSS)
      colnames(GCA_ref1)='risk'
      GCA_ref1$AB_risk_pct=GCA_ref1$risk*100
      GCA_ref1$size='current'
  GCA_ref2=as.data.frame(GCA_PRS_AR_2xSS)
      colnames(GCA_ref2)='risk'
      GCA_ref2$AB_risk_pct=GCA_ref2$risk*100
      GCA_ref2$size='2X'
  GCA_ref3=as.data.frame(GCA_PRS_AR_4xSS)
      colnames(GCA_ref3)='risk'
      GCA_ref3$AB_risk_pct=GCA_ref3$risk*100
      GCA_ref3$size='4X'
  GCA_ref4=as.data.frame(GCA_PRS_AR_InfSS)
      colnames(GCA_ref4)='risk'
      GCA_ref4$AB_risk_pct=GCA_ref4$risk*100
      GCA_ref4$size='Inf'
GCA_ref=rbind(GCA_ref1, GCA_ref2, GCA_ref3, GCA_ref4)
GCA_ref$size <- factor(GCA_ref$size, levels=c('current','2X','4X','Inf'))

load('MSCA_PRS_AR.Rda')
  MSCA_ref1=as.data.frame(MSCA_PRS_AR_CurrentSS)
      colnames(MSCA_ref1)='risk'
      MSCA_ref1$AB_risk_pct=MSCA_ref1$risk*100
      MSCA_ref1$size='current'
  MSCA_ref2=as.data.frame(MSCA_PRS_AR_2xSS)
      colnames(MSCA_ref2)='risk'
      MSCA_ref2$AB_risk_pct=MSCA_ref2$risk*100
      MSCA_ref2$size='2X'
  MSCA_ref3=as.data.frame(MSCA_PRS_AR_4xSS)
      colnames(MSCA_ref3)='risk'
      MSCA_ref3$AB_risk_pct=MSCA_ref3$risk*100
      MSCA_ref3$size='4X'
  MSCA_ref4=as.data.frame(MSCA_PRS_AR_InfSS)
      colnames(MSCA_ref4)='risk'
      MSCA_ref4$AB_risk_pct=MSCA_ref4$risk*100
      MSCA_ref4$size='Inf'
MSCA_ref=rbind(MSCA_ref1, MSCA_ref2, MSCA_ref3, MSCA_ref4)
MSCA_ref$size <- factor(MSCA_ref$size, levels=c('current','2X','4X','Inf'))

load('OCA_PRS_AR.Rda')
  OCA_ref1=as.data.frame(OCA_PRS_AR_CurrentSS)
      colnames(OCA_ref1)='risk'
      OCA_ref1$AB_risk_pct=OCA_ref1$risk*100
      OCA_ref1$size='current'
  OCA_ref2=as.data.frame(OCA_PRS_AR_2xSS)
      colnames(OCA_ref2)='risk'
      OCA_ref2$AB_risk_pct=OCA_ref2$risk*100
      OCA_ref2$size='2X'
  OCA_ref3=as.data.frame(OCA_PRS_AR_4xSS)
      colnames(OCA_ref3)='risk'
      OCA_ref3$AB_risk_pct=OCA_ref3$risk*100
      OCA_ref3$size='4X'
  OCA_ref4=as.data.frame(OCA_PRS_AR_InfSS)
      colnames(OCA_ref4)='risk'
      OCA_ref4$AB_risk_pct=OCA_ref4$risk*100
      OCA_ref4$size='Inf'
OCA_ref=rbind(OCA_ref1, OCA_ref2, OCA_ref3, OCA_ref4)
OCA_ref$size <- factor(OCA_ref$size, levels=c('current','2X','4X','Inf'))

load('LCA_PRS_AR.Rda')
  LCA_ref1=as.data.frame(LCA_PRS_AR_CurrentSS)
      colnames(LCA_ref1)='risk'
      LCA_ref1$AB_risk_pct=LCA_ref1$risk*100
      LCA_ref1$size='current'
  LCA_ref2=as.data.frame(LCA_PRS_AR_2xSS)
      colnames(LCA_ref2)='risk'
      LCA_ref2$AB_risk_pct=LCA_ref2$risk*100
      LCA_ref2$size='2X'
  LCA_ref3=as.data.frame(LCA_PRS_AR_4xSS)
      colnames(LCA_ref3)='risk'
      LCA_ref3$AB_risk_pct=LCA_ref3$risk*100
      LCA_ref3$size='4X'
  LCA_ref4=as.data.frame(LCA_PRS_AR_InfSS)
      colnames(LCA_ref4)='risk'
      LCA_ref4$AB_risk_pct=LCA_ref4$risk*100
      LCA_ref4$size='Inf'
LCA_ref=rbind(LCA_ref1, LCA_ref2, LCA_ref3, LCA_ref4)
LCA_ref$size <- factor(LCA_ref$size, levels=c('current','2X','4X','Inf'))

load('PCA_PRS_AR.Rda')
  PCA_ref1=as.data.frame(PCA_PRS_AR_CurrentSS)
      colnames(PCA_ref1)='risk'
      PCA_ref1$AB_risk_pct=PCA_ref1$risk*100
      PCA_ref1$size='current'
  PCA_ref2=as.data.frame(PCA_PRS_AR_2xSS)
      colnames(PCA_ref2)='risk'
      PCA_ref2$AB_risk_pct=PCA_ref2$risk*100
      PCA_ref2$size='2X'
  PCA_ref3=as.data.frame(PCA_PRS_AR_4xSS)
      colnames(PCA_ref3)='risk'
      PCA_ref3$AB_risk_pct=PCA_ref3$risk*100
      PCA_ref3$size='4X'
  PCA_ref4=as.data.frame(PCA_PRS_AR_InfSS)
      colnames(PCA_ref4)='risk'
      PCA_ref4$AB_risk_pct=PCA_ref4$risk*100
      PCA_ref4$size='Inf'
PCA_ref=rbind(PCA_ref1, PCA_ref2, PCA_ref3, PCA_ref4)
PCA_ref$size <- factor(PCA_ref$size, levels=c('current','2X','4X','Inf'))

load('BCA_PRS_AR.Rda')
  BCA_ref1=as.data.frame(BCA_PRS_AR_CurrentSS)
      colnames(BCA_ref1)='risk'
      BCA_ref1$AB_risk_pct=BCA_ref1$risk*100
      BCA_ref1$size='current'
  BCA_ref2=as.data.frame(BCA_PRS_AR_2xSS)
      colnames(BCA_ref2)='risk'
      BCA_ref2$AB_risk_pct=BCA_ref2$risk*100
      BCA_ref2$size='2X'
  BCA_ref3=as.data.frame(BCA_PRS_AR_4xSS)
      colnames(BCA_ref3)='risk'
      BCA_ref3$AB_risk_pct=BCA_ref3$risk*100
      BCA_ref3$size='4X'
  BCA_ref4=as.data.frame(BCA_PRS_AR_InfSS)
      colnames(BCA_ref4)='risk'
      BCA_ref4$AB_risk_pct=BCA_ref4$risk*100
      BCA_ref4$size='Inf'
BCA_ref=rbind(BCA_ref1, BCA_ref2, BCA_ref3, BCA_ref4)
BCA_ref$size <- factor(BCA_ref$size, levels=c('current','2X','4X','Inf'))

load('HNC_PRS_AR.Rda')
  HNC_ref1=as.data.frame(HNC_PRS_AR_50K)
      colnames(HNC_ref1)='risk'
      HNC_ref1$AB_risk_pct=HNC_ref1$risk*100
      HNC_ref1$size='current'
  HNC_ref2=as.data.frame(HNC_PRS_AR_100K)
      colnames(HNC_ref2)='risk'
      HNC_ref2$AB_risk_pct=HNC_ref2$risk*100
      HNC_ref2$size='2X'
  HNC_ref3=as.data.frame(HNC_PRS_AR_200K)
      colnames(HNC_ref3)='risk'
      HNC_ref3$AB_risk_pct=HNC_ref3$risk*100
      HNC_ref3$size='4X'
  HNC_ref4=as.data.frame(HNC_PRS_AR_InfSS)
      colnames(HNC_ref4)='risk'
      HNC_ref4$AB_risk_pct=HNC_ref4$risk*100
      HNC_ref4$size='Inf'
HNC_ref=rbind(HNC_ref1, HNC_ref2, HNC_ref3, HNC_ref4)
HNC_ref$size <- factor(HNC_ref$size, levels=c('current','2X','4X','Inf'))

save(CLL_ref, ECA_ref, TCA_ref, HNC_ref, PACA_ref, RCA_ref, GCA_ref, MSCA_ref, CRCA_ref, UECA_ref, OCA_ref, LCA_ref, PCA_ref, BCA_ref,
     file='cross_cancer_ref_data.Rda')

# Figure 5 ####
load('cross_cancer_ref_data.Rda')
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/Projection plots')

p1<-AR_dens_plot2(CLL_ref, site='CLL', x.lims=c(-0.1,2), y.lims=c(0,4.8), with_key=T)
p2<-AR_dens_plot2(ECA_ref, site='Esophageal', x.lims=c(-0.03,1.3), x.breaks=seq(0,1.2,0.3))
p3<-AR_dens_plot2(TCA_ref, site='Testicular', x.lims=c(-0.02,0.45), x.breaks=seq(0,0.4,0.1), y.breaks=seq(0,30,5))
p4<-AR_dens_plot2(HNC_ref, site='Oropharyngeal', x.lims=c(-0.03,2.5), x.breaks=seq(0,2.4,0.8), y.breaks=seq(0,8,2))

p5<-AR_dens_plot2(PACA_ref, site='Pancreas', x.lims=c(-0.1,3.5))
p6<-AR_dens_plot2(RCA_ref, site='Renal', x.lims=c(-0.1,4.5))
p7<-AR_dens_plot2(GCA_ref, site='Glioma',x.lims=c(-0.05,1.3), x.breaks=seq(0,1.2,0.3))
p8<-AR_dens_plot2(MSCA_ref, site='Melanoma', x.lims=c(-0.3,13), x.breaks=seq(0,12,3))

p9<-AR_dens_plot2(CRCA_ref, site='Colorectal', x.lims=c(-0.1,7), y.lims=c(0,1.35), y.breaks=seq(0,1.2,0.3))
p10<-AR_dens_plot2(UECA_ref, site='Endometrial', x.lims=c(0,6), x.breaks=seq(0,6,2), y.lims=c(0,1.35), y.breaks=seq(0,1.2,0.3))
p11<-AR_dens_plot2(OCA_ref, site='Ovarian', x.lims=c(0,2.1), y.lims=c(0,3))
p12<-AR_dens_plot2(LCA_ref, site='Lung', x.lims=c(-0.1,12))
p13<-AR_dens_plot2(PCA_ref, site='Prostate', x.lims=c(-1,32))
p14<-AR_dens_plot2(BCA_ref, site='Breast', x.lims=c(-1,40), y.lims=c(0,0.126), y.breaks=seq(0,0.12,0.03))

all_plot=plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, align = "hv", nrow=4, ncol=4)

final=ggdraw(all_plot) + theme(plot.margin=margin(1, 1, 1.5, 1.5, "cm")) +
  draw_label("x: Average residual lifetime risk (%)", fontface='plain', x=0.5, y=-0.005, hjust=0.5, vjust=0.5, size=26) +
  draw_label("y: Density", fontface='plain', angle=90, x=-0.005, y=0.5, hjust=0.5, vjust=0.5, size=26)
  
g = ggplotGrob(final)
g$layout$clip[g$layout$name == "panel"] = "off"

#png(file='Figure5_030420.png', width=20.2, height=20, units="in", res=500)
#grid.draw(g)
#dev.off()

pdf(file='Figure5_030420.pdf', width=20.2, height=20)
grid.draw(g)
dev.off()

# Supplementary Figure 5: GROUP 1 ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
#     CLL ####
load('CLL_PRS_AR_AgeStrat.Rda')

CLL_PRS_AR_30to34_v2=melt(CLL_PRS_AR_30to34)
colnames(CLL_PRS_AR_30to34_v2)=c('size','risk')
CLL_PRS_AR_30to34_v2$size <- factor(CLL_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CLL_PRS_AR_40to44_v2=melt(CLL_PRS_AR_40to44)
colnames(CLL_PRS_AR_40to44_v2)=c('size','risk')
CLL_PRS_AR_40to44_v2$size <- factor(CLL_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CLL_PRS_AR_50to54_v2=melt(CLL_PRS_AR_50to54)
colnames(CLL_PRS_AR_50to54_v2)=c('size','risk')
CLL_PRS_AR_50to54_v2$size <- factor(CLL_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CLL_PRS_AR_60to64_v2=melt(CLL_PRS_AR_60to64)
colnames(CLL_PRS_AR_60to64_v2)=c('size','risk')
CLL_PRS_AR_60to64_v2$size <- factor(CLL_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CLL30=ggplot(CLL_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,2), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,5.1), breaks=seq(0,5,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.direction="vertical", 
        legend.position=c(.4,.8), legend.title=element_blank(), legend.text=element_text(size=26), 
        legend.key.size=unit(0.8,'cm'), plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CLL40=ggplot(CLL_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,2), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,5.1), breaks=seq(0,5,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CLL50=ggplot(CLL_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,2), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,5.1), breaks=seq(0,5,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CLL60=ggplot(CLL_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
   scale_x_continuous(" ", limits=c(-0.1,2), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,5.1), breaks=seq(0,5,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     ESOPHAGEAL ####
load('ECA_PRS_AR_AgeStrat.Rda')

ECA_PRS_AR_30to34_v2=melt(ECA_PRS_AR_30to34)
colnames(ECA_PRS_AR_30to34_v2)=c('size','risk')
ECA_PRS_AR_30to34_v2$size <- factor(ECA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

ECA_PRS_AR_40to44_v2=melt(ECA_PRS_AR_40to44)
colnames(ECA_PRS_AR_40to44_v2)=c('size','risk')
ECA_PRS_AR_40to44_v2$size <- factor(ECA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

ECA_PRS_AR_50to54_v2=melt(ECA_PRS_AR_50to54)
colnames(ECA_PRS_AR_50to54_v2)=c('size','risk')
ECA_PRS_AR_50to54_v2$size <- factor(ECA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

ECA_PRS_AR_60to64_v2=melt(ECA_PRS_AR_60to64)
colnames(ECA_PRS_AR_60to64_v2)=c('size','risk')
ECA_PRS_AR_60to64_v2$size <- factor(ECA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

ECA30=ggplot(ECA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.55), breaks=seq(0,1.5,0.5), expand = c(0.03,0)) +
   scale_y_continuous(" ", limits=c(0,12), breaks=seq(0,12,3), expand = c(0, 0)) +
   theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

ECA40=ggplot(ECA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.55), breaks=seq(0,1.5,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,12), breaks=seq(0,12,3), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

ECA50=ggplot(ECA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.55), breaks=seq(0,1.5,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,12), breaks=seq(0,12,3), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

ECA60=ggplot(ECA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.55), breaks=seq(0,1.5,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,12), breaks=seq(0,12,3), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     TESTICULAR ####
load('TCA_PRS_AR_AgeStrat.Rda')

TCA_PRS_AR_30to34_v2=melt(TCA_PRS_AR_30to34)
colnames(TCA_PRS_AR_30to34_v2)=c('size','risk')
TCA_PRS_AR_30to34_v2$size <- factor(TCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

TCA_PRS_AR_40to44_v2=melt(TCA_PRS_AR_40to44)
colnames(TCA_PRS_AR_40to44_v2)=c('size','risk')
TCA_PRS_AR_40to44_v2$size <- factor(TCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

TCA_PRS_AR_50to54_v2=melt(TCA_PRS_AR_50to54)
colnames(TCA_PRS_AR_50to54_v2)=c('size','risk')
TCA_PRS_AR_50to54_v2$size <- factor(TCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

TCA_PRS_AR_60to64_v2=melt(TCA_PRS_AR_60to64)
colnames(TCA_PRS_AR_60to64_v2)=c('size','risk')
TCA_PRS_AR_60to64_v2$size <- factor(TCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

TCA30=ggplot(TCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.2), breaks=seq(0,1.2,0.4), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,10), breaks=seq(0,10,2), expand = c(0, 0)) +
   theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

TCA40=ggplot(TCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.2), breaks=seq(0,1.2,0.4), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,18.5), breaks=seq(0,15,5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

TCA50=ggplot(TCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.2), breaks=seq(0,1.2,0.4), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,40), breaks=seq(0,40,10), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

TCA60=ggplot(TCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.2), breaks=seq(0,1.2,0.4), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,100), breaks=seq(0,100,25), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     OROPHARYNGEAL ####
load('HNC_PRS_AR_AgeStrat.Rda')

HNC_PRS_AR_30to34_v2=melt(HNC_PRS_AR_30to34)
colnames(HNC_PRS_AR_30to34_v2)=c('size','risk')
HNC_PRS_AR_30to34_v2$size <- factor(HNC_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

HNC_PRS_AR_40to44_v2=melt(HNC_PRS_AR_40to44)
colnames(HNC_PRS_AR_40to44_v2)=c('size','risk')
HNC_PRS_AR_40to44_v2$size <- factor(HNC_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

HNC_PRS_AR_50to54_v2=melt(HNC_PRS_AR_50to54)
colnames(HNC_PRS_AR_50to54_v2)=c('size','risk')
HNC_PRS_AR_50to54_v2$size <- factor(HNC_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

HNC_PRS_AR_60to64_v2=melt(HNC_PRS_AR_60to64)
colnames(HNC_PRS_AR_60to64_v2)=c('size','risk')
HNC_PRS_AR_60to64_v2$size <- factor(HNC_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

HNC30=ggplot(HNC_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,3), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6.5), breaks=seq(0,8,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" 50K cases, 50K controls"," 100K cases, 100K controls"," 200K cases, 200K controls"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" 50K cases, 50K controls"," 100K cases, 100K controls"," 200K cases, 200K controls"," Infinite sample size"))

HNC40=ggplot(HNC_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,3), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6.5), breaks=seq(0,8,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size"))

HNC50=ggplot(HNC_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,3), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6.5), breaks=seq(0,8,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size"))

HNC60=ggplot(HNC_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,3), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6.5), breaks=seq(0,8,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                    labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), 
                     labels=c(" Current sample size"," Current sample size x2"," Current sample size x4"," Infinite sample size"))

#     PANCREATIC ####
load('PACA_PRS_AR_AgeStrat.Rda')

PACA_PRS_AR_30to34_v2=melt(PACA_PRS_AR_30to34)
colnames(PACA_PRS_AR_30to34_v2)=c('size','risk')
PACA_PRS_AR_30to34_v2$size <- factor(PACA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PACA_PRS_AR_40to44_v2=melt(PACA_PRS_AR_40to44)
colnames(PACA_PRS_AR_40to44_v2)=c('size','risk')
PACA_PRS_AR_40to44_v2$size <- factor(PACA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PACA_PRS_AR_50to54_v2=melt(PACA_PRS_AR_50to54)
colnames(PACA_PRS_AR_50to54_v2)=c('size','risk')
PACA_PRS_AR_50to54_v2$size <- factor(PACA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PACA_PRS_AR_60to64_v2=melt(PACA_PRS_AR_60to64)
colnames(PACA_PRS_AR_60to64_v2)=c('size','risk')
PACA_PRS_AR_60to64_v2$size <- factor(PACA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PACA30=ggplot(PACA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,3.4), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PACA40=ggplot(PACA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,3.4), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PACA50=ggplot(PACA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,3.4), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PACA60=ggplot(PACA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,3.4), breaks=seq(0,3,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#           GROUP 1 PLOT GRID ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/Projection plots')

all_plot=plot_grid(CLL30,CLL40,CLL50,CLL60,ECA30,ECA40,ECA50,ECA60,TCA30,TCA40,TCA50,TCA60,
                   HNC30,HNC40,HNC50,HNC60,PACA30,PACA40,PACA50,PACA60, align = "hv", nrow=5, ncol=4)

final=ggdraw(all_plot) + theme(plot.margin=margin(3, 2, 2, 9, "cm")) +
  draw_label("CLL", fontface='plain', x=0.015, y=0.991, hjust=1, vjust=0.5, size=30) +
  draw_label("Esophageal", fontface='plain', x=0.015, y=0.793, hjust=1, vjust=0.5, size=30) +
  draw_label("Testicular", fontface='plain', x=0.015, y=0.592, hjust=1, vjust=0.5, size=30) +
  draw_label("Oropharyngeal", fontface='plain', x=0.015, y=0.394, hjust=1, vjust=0.5, size=30) +
  draw_label("Pancreas", fontface='plain', x=0.015, y=0.192, hjust=1, vjust=0.5, size=30) +
  draw_label("Ages 30 to 34 years", fontface='plain', x=0.145, y=1.03, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 40 to 44 years", fontface='plain', x=0.395, y=1.03, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 50 to 54 years", fontface='plain', x=0.646, y=1.03, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 60 to 64 years", fontface='plain', x=0.897, y=1.03, hjust=0.5, vjust=1, size=30) +
  draw_label("x: Lifetime absolute risk (%)", fontface='plain', x=0.5, y=-0.005, hjust=0.5, vjust=0.5, size=30) +
  draw_label("y: Density", fontface='plain', angle=90, x=-0.15, y=0.5, hjust=0.5, vjust=0.5, size=30)

g = ggplotGrob(final)
g$layout$clip[g$layout$name == "panel"] = "off"
 
png(file='SuppFig5_Group1_72319.png', width=25, height=28, units="in", res=500)
grid.draw(g)
dev.off()

# Supplementary Figure 6: GROUP 2 ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
#     RENAL ####
load('RCA_PRS_AR_AgeStrat.Rda')

RCA_PRS_AR_30to34_v2=melt(RCA_PRS_AR_30to34)
colnames(RCA_PRS_AR_30to34_v2)=c('size','risk')
RCA_PRS_AR_30to34_v2$size <- factor(RCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

RCA_PRS_AR_40to44_v2=melt(RCA_PRS_AR_40to44)
colnames(RCA_PRS_AR_40to44_v2)=c('size','risk')
RCA_PRS_AR_40to44_v2$size <- factor(RCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

RCA_PRS_AR_50to54_v2=melt(RCA_PRS_AR_50to54)
colnames(RCA_PRS_AR_50to54_v2)=c('size','risk')
RCA_PRS_AR_50to54_v2$size <- factor(RCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

RCA_PRS_AR_60to64_v2=melt(RCA_PRS_AR_60to64)
colnames(RCA_PRS_AR_60to64_v2)=c('size','risk')
RCA_PRS_AR_60to64_v2$size <- factor(RCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

RCA30=ggplot(RCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,4.3), breaks=seq(0,4,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.direction="vertical", 
        legend.position=c(.4,.8), legend.title=element_blank(), legend.text=element_text(size=26), 
        legend.key.size=unit(0.7,'cm'), plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

RCA40=ggplot(RCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,4.3), breaks=seq(0,4,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

RCA50=ggplot(RCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,4.3), breaks=seq(0,4,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

RCA60=ggplot(RCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,4.3), breaks=seq(0,4,1), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,2.15), breaks=seq(0,2,0.5), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     GLIOMA ####
load('GCA_PRS_AR_AgeStrat.Rda')

GCA_PRS_AR_30to34_v2=melt(GCA_PRS_AR_30to34)
colnames(GCA_PRS_AR_30to34_v2)=c('size','risk')
GCA_PRS_AR_30to34_v2$size <- factor(GCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

GCA_PRS_AR_40to44_v2=melt(GCA_PRS_AR_40to44)
colnames(GCA_PRS_AR_40to44_v2)=c('size','risk')
GCA_PRS_AR_40to44_v2$size <- factor(GCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

GCA_PRS_AR_50to54_v2=melt(GCA_PRS_AR_50to54)
colnames(GCA_PRS_AR_50to54_v2)=c('size','risk')
GCA_PRS_AR_50to54_v2$size <- factor(GCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

GCA_PRS_AR_60to64_v2=melt(GCA_PRS_AR_60to64)
colnames(GCA_PRS_AR_60to64_v2)=c('size','risk')
GCA_PRS_AR_60to64_v2$size <- factor(GCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

GCA30=ggplot(GCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.1), breaks=seq(0,0.9,0.3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6), breaks=seq(0,6,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

GCA40=ggplot(GCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.1), breaks=seq(0,0.9,0.3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6), breaks=seq(0,6,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

GCA50=ggplot(GCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.1), breaks=seq(0,0.9,0.3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6), breaks=seq(0,6,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

GCA60=ggplot(GCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,1.1), breaks=seq(0,0.9,0.3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,6), breaks=seq(0,6,2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     MELANOMA ####
load('MSCA_PRS_AR_AgeStrat.Rda')

MSCA_PRS_AR_30to34_v2=melt(MSCA_PRS_AR_30to34)
colnames(MSCA_PRS_AR_30to34_v2)=c('size','risk')
MSCA_PRS_AR_30to34_v2$size <- factor(MSCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

MSCA_PRS_AR_40to44_v2=melt(MSCA_PRS_AR_40to44)
colnames(MSCA_PRS_AR_40to44_v2)=c('size','risk')
MSCA_PRS_AR_40to44_v2$size <- factor(MSCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

MSCA_PRS_AR_50to54_v2=melt(MSCA_PRS_AR_50to54)
colnames(MSCA_PRS_AR_50to54_v2)=c('size','risk')
MSCA_PRS_AR_50to54_v2$size <- factor(MSCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

MSCA_PRS_AR_60to64_v2=melt(MSCA_PRS_AR_60to64)
colnames(MSCA_PRS_AR_60to64_v2)=c('size','risk')
MSCA_PRS_AR_60to64_v2$size <- factor(MSCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

MSCA30=ggplot(MSCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.3,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.5), breaks=seq(0,0.5,0.1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

MSCA40=ggplot(MSCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.3,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.5), breaks=seq(0,0.5,0.1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

MSCA50=ggplot(MSCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.3,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.5), breaks=seq(0,0.5,0.1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

MSCA60=ggplot(MSCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.3,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.5), breaks=seq(0,0.5,0.1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     COLORECTAL ####
load('CRCA_PRS_AR_AgeStrat.Rda')

CRCA_PRS_AR_30to34_v2=melt(CRCA_PRS_AR_30to34)
colnames(CRCA_PRS_AR_30to34_v2)=c('size','risk')
CRCA_PRS_AR_30to34_v2$size <- factor(CRCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CRCA_PRS_AR_40to44_v2=melt(CRCA_PRS_AR_40to44)
colnames(CRCA_PRS_AR_40to44_v2)=c('size','risk')
CRCA_PRS_AR_40to44_v2$size <- factor(CRCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CRCA_PRS_AR_50to54_v2=melt(CRCA_PRS_AR_50to54)
colnames(CRCA_PRS_AR_50to54_v2)=c('size','risk')
CRCA_PRS_AR_50to54_v2$size <- factor(CRCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CRCA_PRS_AR_60to64_v2=melt(CRCA_PRS_AR_60to64)
colnames(CRCA_PRS_AR_60to64_v2)=c('size','risk')
CRCA_PRS_AR_60to64_v2$size <- factor(CRCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

CRCA30=ggplot(CRCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,8), expand = c(0, 0)) +
  scale_y_continuous(" ", limits=c(0,1.6), breaks=seq(0,1.6,0.4), expand = c(0.03,0)) +
   theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CRCA40=ggplot(CRCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,8), expand = c(0, 0)) +
  scale_y_continuous(" ", limits=c(0,1.6), breaks=seq(0,1.6,0.4), expand = c(0.03,0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CRCA50=ggplot(CRCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,8), expand = c(0, 0)) +
  scale_y_continuous(" ", limits=c(0,1.6), breaks=seq(0,1.6,0.4), expand = c(0.03,0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

CRCA60=ggplot(CRCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.05,8), expand = c(0, 0)) +
  scale_y_continuous(" ", limits=c(0,1.6), breaks=seq(0,1.6,0.4), expand = c(0.03,0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     ENDOMETRIAL ####
load('UECA_PRS_AR_AgeStrat.Rda')

UECA_PRS_AR_30to34_v2=melt(UECA_PRS_AR_30to34)
colnames(UECA_PRS_AR_30to34_v2)=c('size','risk')
UECA_PRS_AR_30to34_v2$size <- factor(UECA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

UECA_PRS_AR_40to44_v2=melt(UECA_PRS_AR_40to44)
colnames(UECA_PRS_AR_40to44_v2)=c('size','risk')
UECA_PRS_AR_40to44_v2$size <- factor(UECA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

UECA_PRS_AR_50to54_v2=melt(UECA_PRS_AR_50to54)
colnames(UECA_PRS_AR_50to54_v2)=c('size','risk')
UECA_PRS_AR_50to54_v2$size <- factor(UECA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

UECA_PRS_AR_60to64_v2=melt(UECA_PRS_AR_60to64)
colnames(UECA_PRS_AR_60to64_v2)=c('size','risk')
UECA_PRS_AR_60to64_v2$size <- factor(UECA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

UECA30=ggplot(UECA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,7), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,1.6), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

UECA40=ggplot(UECA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,7), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,1.6), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

UECA50=ggplot(UECA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,7), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,1.6), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

UECA60=ggplot(UECA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,7), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,1.6), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     OVARIAN ####
load('OCA_PRS_AR_AgeStrat.Rda')

OCA_PRS_AR_30to34_v2=melt(OCA_PRS_AR_30to34)
colnames(OCA_PRS_AR_30to34_v2)=c('size','risk')
OCA_PRS_AR_30to34_v2$size <- factor(OCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

OCA_PRS_AR_40to44_v2=melt(OCA_PRS_AR_40to44)
colnames(OCA_PRS_AR_40to44_v2)=c('size','risk')
OCA_PRS_AR_40to44_v2$size <- factor(OCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

OCA_PRS_AR_50to54_v2=melt(OCA_PRS_AR_50to54)
colnames(OCA_PRS_AR_50to54_v2)=c('size','risk')
OCA_PRS_AR_50to54_v2$size <- factor(OCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

OCA_PRS_AR_60to64_v2=melt(OCA_PRS_AR_60to64)
colnames(OCA_PRS_AR_60to64_v2)=c('size','risk')
OCA_PRS_AR_60to64_v2$size <- factor(OCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

OCA30=ggplot(OCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,2.3), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,3.6), breaks=seq(0,3,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

OCA40=ggplot(OCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,2.3), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,3.6), breaks=seq(0,3,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

OCA50=ggplot(OCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,2.3), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,3.6), breaks=seq(0,3,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

OCA60=ggplot(OCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.01,2.3), breaks=seq(0,2,0.5), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,3.6), breaks=seq(0,3,1), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#           GROUP 2 PLOT GRID ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/Projection plots')

all_plot=plot_grid(RCA30,RCA40,RCA50,RCA60,GCA30,GCA40,GCA50,GCA60,MSCA30,MSCA40,MSCA50,MSCA60,
                   CRCA30,CRCA40,CRCA50,CRCA60,UECA30,UECA40,UECA50,UECA60,OCA30,OCA40,OCA50,OCA60,
                   align = "hv", nrow=6, ncol=4)

final=ggdraw(all_plot) + theme(plot.margin=margin(3, 2, 2, 8, "cm")) +
  draw_label("Renal", fontface='plain', x=0.015, y=0.992, hjust=1, vjust=0.5, size=30) +
  draw_label("Glioma", fontface='plain', x=0.015, y=0.826, hjust=1, vjust=0.5, size=30) +
  draw_label("Melanoma", fontface='plain', x=0.015, y=0.659, hjust=1, vjust=0.5, size=30) +
  draw_label("Colorectal", fontface='plain', x=0.015, y=0.492, hjust=1, vjust=0.5, size=30) +
  draw_label("Endometrial", fontface='plain', x=0.015, y=0.325, hjust=1, vjust=0.5, size=30) +
  draw_label("Ovarian", fontface='plain', x=0.015, y=0.16, hjust=1, vjust=0.5, size=30) +
  draw_label("Ages 30 to 34 years", fontface='plain', x=0.138, y=1.025, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 40 to 44 years", fontface='plain', x=0.388, y=1.025, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 50 to 54 years", fontface='plain', x=0.639, y=1.025, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 60 to 64 years", fontface='plain', x=0.890, y=1.025, hjust=0.5, vjust=1, size=30) +
  draw_label("x: Lifetime absolute risk (%)", fontface='plain', x=0.5, y=-0.005, hjust=0.5, vjust=0.5, size=30) +
  draw_label("y: Density", fontface='plain', angle=90, x=-0.13, y=0.5, hjust=0.5, vjust=0.5, size=30)

g = ggplotGrob(final)
g$layout$clip[g$layout$name == "panel"] = "off"
 
png(file='SuppFig6_Group2_52019.png', width=25, height=33, units="in", res=500)
grid.draw(g)
dev.off()

# Supplementary Figure 7: GROUP 3 ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/R files')
#     LUNG ####
load('LCA_PRS_AR_AgeStrat.Rda')

LCA_PRS_AR_30to34_v2=melt(LCA_PRS_AR_30to34)
colnames(LCA_PRS_AR_30to34_v2)=c('size','risk')
LCA_PRS_AR_30to34_v2$size <- factor(LCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

LCA_PRS_AR_40to44_v2=melt(LCA_PRS_AR_40to44)
colnames(LCA_PRS_AR_40to44_v2)=c('size','risk')
LCA_PRS_AR_40to44_v2$size <- factor(LCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

LCA_PRS_AR_50to54_v2=melt(LCA_PRS_AR_50to54)
colnames(LCA_PRS_AR_50to54_v2)=c('size','risk')
LCA_PRS_AR_50to54_v2$size <- factor(LCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

LCA_PRS_AR_60to64_v2=melt(LCA_PRS_AR_60to64)
colnames(LCA_PRS_AR_60to64_v2)=c('size','risk')
LCA_PRS_AR_60to64_v2$size <- factor(LCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

LCA30=ggplot(LCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.67), breaks=seq(0,0.6,0.2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.direction="vertical", 
        legend.position=c(.4,.8), legend.title=element_blank(), legend.text=element_text(size=26), 
        legend.key.size=unit(0.7,'cm'), plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

LCA40=ggplot(LCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.67), breaks=seq(0,0.6,0.2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

LCA50=ggplot(LCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.67), breaks=seq(0,0.6,0.2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

LCA60=ggplot(LCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-0.1,14), breaks=seq(0,12,3), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.67), breaks=seq(0,0.6,0.2), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#     PROSTATE ####
load('PCA_PRS_AR_AgeStrat.Rda')

PCA_PRS_AR_30to34_v2=melt(PCA_PRS_AR_30to34)
colnames(PCA_PRS_AR_30to34_v2)=c('size','risk')
PCA_PRS_AR_30to34_v2$size <- factor(PCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PCA_PRS_AR_40to44_v2=melt(PCA_PRS_AR_40to44)
colnames(PCA_PRS_AR_40to44_v2)=c('size','risk')
PCA_PRS_AR_40to44_v2$size <- factor(PCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PCA_PRS_AR_50to54_v2=melt(PCA_PRS_AR_50to54)
colnames(PCA_PRS_AR_50to54_v2)=c('size','risk')
PCA_PRS_AR_50to54_v2$size <- factor(PCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PCA_PRS_AR_60to64_v2=melt(PCA_PRS_AR_60to64)
colnames(PCA_PRS_AR_60to64_v2)=c('size','risk')
PCA_PRS_AR_60to64_v2$size <- factor(PCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

PCA30=ggplot(PCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,35), breaks=seq(0,30,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.17), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PCA40=ggplot(PCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,35), breaks=seq(0,30,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.17), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PCA50=ggplot(PCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,35), breaks=seq(0,30,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.17), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

PCA60=ggplot(PCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,35), breaks=seq(0,30,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.17), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))


#     BREAST ####
load('BCA_PRS_AR_AgeStrat.Rda')

BCA_PRS_AR_30to34_v2=melt(BCA_PRS_AR_30to34)
colnames(BCA_PRS_AR_30to34_v2)=c('size','risk')
BCA_PRS_AR_30to34_v2$size <- factor(BCA_PRS_AR_30to34_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

BCA_PRS_AR_40to44_v2=melt(BCA_PRS_AR_40to44)
colnames(BCA_PRS_AR_40to44_v2)=c('size','risk')
BCA_PRS_AR_40to44_v2$size <- factor(BCA_PRS_AR_40to44_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

BCA_PRS_AR_50to54_v2=melt(BCA_PRS_AR_50to54)
colnames(BCA_PRS_AR_50to54_v2)=c('size','risk')
BCA_PRS_AR_50to54_v2$size <- factor(BCA_PRS_AR_50to54_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

BCA_PRS_AR_60to64_v2=melt(BCA_PRS_AR_60to64)
colnames(BCA_PRS_AR_60to64_v2)=c('size','risk')
BCA_PRS_AR_60to64_v2$size <- factor(BCA_PRS_AR_60to64_v2$size,
                                        levels=c('AR_CurrentSS','AR_2xSS','AR_4xSS','AR_InfSS'))

BCA30=ggplot(BCA_PRS_AR_30to34_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,42), breaks=seq(0,40,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.16), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

BCA40=ggplot(BCA_PRS_AR_40to44_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,42), breaks=seq(0,40,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.16), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

BCA50=ggplot(BCA_PRS_AR_50to54_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,42), breaks=seq(0,40,10), expand = c(0.03,0)) +
  scale_y_continuous(" ", limits=c(0,0.16), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

BCA60=ggplot(BCA_PRS_AR_60to64_v2, aes(x=risk, colour=size, fill=size, group=size)) + 
  geom_density(aes(group=size, colour=size, fill=size), alpha=0.4) +
  scale_x_continuous(" ", limits=c(-1,42), breaks=seq(0,40,10), expand = c(0, 0)) +
  scale_y_continuous(" ", limits=c(0,0.16), breaks=seq(0,0.16,0.04), expand = c(0, 0)) +
  theme(plot.title = element_text(size=26, colour = "black", vjust=4, hjust=0), axis.text.x=element_text(size=26),axis.text.y=element_text(size=26), 
        axis.title.x = element_text(size=26), axis.title.y = element_text(size=26), legend.position = "none", plot.margin=unit(c(1,0.75,0.3,0.3),"cm")) +
  scale_fill_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite")) +
  scale_color_manual(values=c("royalblue","firebrick","forestgreen","lightsalmon"), labels=c(" Current"," Double"," Quadruple"," Infinite"))

#           GROUP 3 PLOT GRID ####
setwd('/Users/wilcoxan/Google Drive/Crosscancer Effect Size/Risk Projections/Projection plots')

all_plot=plot_grid(LCA30,LCA40,LCA50,LCA60,PCA30,PCA40,PCA50,PCA60,BCA30,BCA40,BCA50,BCA60, align = "hv", nrow=3, ncol=4)

final=ggdraw(all_plot) + theme(plot.margin=margin(3, 2, 2, 7, "cm")) +
  draw_label("Lung", fontface='plain', x=0.015, y=0.99, hjust=1, vjust=0.5, size=30) +
  draw_label("Prostate", fontface='plain', x=0.015, y=0.66, hjust=1, vjust=0.5, size=30) +
  draw_label("Breast", fontface='plain', x=0.015, y=0.33, hjust=1, vjust=0.5, size=30) +
  draw_label("Ages 30 to 34 years", fontface='plain', x=0.145, y=1.05, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 40 to 44 years", fontface='plain', x=0.395, y=1.05, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 50 to 54 years", fontface='plain', x=0.646, y=1.05, hjust=0.5, vjust=1, size=30) +
  draw_label("Ages 60 to 64 years", fontface='plain', x=0.897, y=1.05, hjust=0.5, vjust=1, size=30) +
  draw_label("x: Lifetime absolute risk (%)", fontface='plain', x=0.5, y=-0.005, hjust=0.5, vjust=0.5, size=30) +
  draw_label("y: Density", fontface='plain', angle=90, x=-0.1, y=0.5, hjust=0.5, vjust=0.5, size=30)


g = ggplotGrob(final)
g$layout$clip[g$layout$name == "panel"] = "off"
 
png(file='SuppFig7_Group3_52019.png', width=25, height=18, units="in", res=500)
grid.draw(g)
dev.off()
