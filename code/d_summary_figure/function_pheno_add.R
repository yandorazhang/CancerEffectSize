#-------------------------------------------------------------------
# Update Date: 05/05/2019
# Create Date: 05/05/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
pheno_add <- function(est,n,summarydata_OutlierIndep=NULL,gwas.significance=5e-8){
  
  if(length(est)==3)components=2
  if(length(est)==5)components=3
  
  if(components==2){
    pic = est[1]
    sig = sqrt(est[2])
    den <- function(x){return(dnorm(x/sig)/sig )}
    herit0 <- pic*M*sig^2
  }
  
  if(components==3){
    pic = est[1]
    p1 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return(p1 * dnorm(x/s1)/s1 + (1-p1)*dnorm(x/s2) /s2)}
    herit0 <- pic*M*(p1*est[3] + (1-p1)*est[4])
  }
  
  tem0 <- function(x){return(x^2*den(x))}
  c_gwsignificance = abs(qnorm(gwas.significance/2))
  pow <- function(x){return(1 - pnorm(c_gwsignificance - sqrt(n)*x) + pnorm(-c_gwsignificance - sqrt(n)*x) )}
  
  
  if(is.null(summarydata_OutlierIndep)){
    pheno.var1 = 0
  }
  
  if(!is.null(summarydata_OutlierIndep)){
    beta = as.numeric(as.character(summarydata_OutlierIndep$z))/sqrt(as.numeric(as.character(summarydata_OutlierIndep$n)))
    n_OutlierIndep = length(beta)
    tau = sqrt(1/as.numeric(as.character(summarydata_OutlierIndep$n)))
    pheno.var1 = sum( (beta^2 - tau^2) * pow(beta) )
  }
  
  return(pheno.var1)
}
