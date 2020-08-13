README.md
====
  

1. `code`: contains all code files to get the figures/tables in the paper. 
   * `a_read_summarydata`:  read the raw summary GWAS data, then convert it to some format that can be load to GENESIS R package. Take prostate cancer as an example, the raw data can be downloaded from http://practical.icr.ac.uk/blog/?page_id=8164. The summary GWAS data is summarized into `summarydata.RData`. 
   * `b_clump_outlier`: doing LD clumping for outlier SNPs with effect sizes chi-square larger than 80, i.e., take out those SNPs with chi-square larger than 80, as well as their correlated SNPs. Denote  Also calculate the heritability of those outlier SNPs.   The summary data after removing the outlier SNPs are stored into `summarydata_RemoveOutlier.RData`.  The heritability that are explained by those removed outlier SNPs are stored into `herit_OutlierIndep.RData`. 
   * `b1_LDclump_gw_signif`: doing LD clumping to get genome-wide significant SNPs information. 
   * `c_fit_GENESIS`: `fit_GENESIS` fit GENESIS model on `summarydata_RemoveOutlier.RData`, and stored object into `fit2_RemoveOutlier_pic*.RData`, `bestfit2_RemoveOutlier.RData`, and `fit3_RemoveOutlier.RData`. If there is no outlier for a particular GWAS data, then `fit_GENESIS 2` can be fitted on `summarydata.RData`, and the fitted object are stored in `fit2_pic*.RData`, `bestfit2.RData`, and `fit3.RData`. 
   * `d_summary_figure`: contains the R code to get the figures/tables in the paper. 
   * `e_LifeTimeRisk_Fig5_SF567`: contains the R code to get Figure 5 and Supplementary Figure 5, 6 & 7.  
   
2.   `data_new` contains the `.RData` files running from `code/a_read_summarydata` and `code/b_clump_outlier`. 

3.  `data_files` contains the REFERENCE panel, i.e., 1000 Genome Hapmap3 common SNPs information.

4. `data_samplesize` contains `cancer_sample_size.csv` file which contains the sample size information of the cancer GWAS.  
   
2. `genesis_result_new` saved GENESIS results running from `code/c_fit_GENESIS`

5. `result_png_csv_new_clump` constains the summarized figures/tables running from `code/d_summary_figure`. 
 

## Citation

Zhang, Y.D., Hurson, A.N., Zhang, H. et al. Assessment of polygenic architecture and risk prediction based on common variants across fourteen cancers. Nat Commun 11, 3353 (2020). [https://doi.org/10.1038/s41467-020-16483-3](https://doi.org/10.1038/s41467-020-16483-3)

## Contact the Author

Yan Dora Zhang (doraz@hku.hk or yandorazhang@gmail.com)

(Updated on July 3, 2020)



