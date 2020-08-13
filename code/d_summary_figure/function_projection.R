#-------------------------------------------------------------------
# Update Date: 04/30/2019
# Create Date: 04/21/2019
# Author: Yan (Dora) Zhang
#-------------------------------------------------------------------
projection <- function(est,v=NULL,n,summarydata_OutlierIndep,herit_OutlierIndep,
                       gwas.significance=5e-8,tol=c(1e-12,1e-15),
                       M=1070777,CI=FALSE,nsim=1000,CI.coverage=0.95,seeds=123){
  
  # within function
  pp <- function(est,n,summarydata_OutlierIndep,herit_OutlierIndep,gwas.significance=5e-8,tol=c(1e-12,1e-15),M=1070777){
    
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
    
    heritall <- herit0 + herit_OutlierIndep
    tem0 <- function(x){return(x^2*den(x))}
    c_gwsignificance = abs(qnorm(gwas.significance/2))
    pow <- function(x){return(1 - pnorm(c_gwsignificance - sqrt(n)*x) + pnorm(-c_gwsignificance - sqrt(n)*x) )}
    tem <- function(x){return(pow(x)*den(x))}
    tem1 <- function(x){return(pow(x)*den(x)*x^2)}
    
    Numdiscoveries = M*pic * integrate(tem, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]
    pheno.var =  M*pic * integrate(tem1, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]
    
    
    if(!is.null(summarydata_OutlierIndep)){
      beta = as.numeric(as.character(summarydata_OutlierIndep$z))/sqrt(as.numeric(as.character(summarydata_OutlierIndep$n)))
      n_OutlierIndep = length(beta)
      Numdiscoveries1 = sum(pow(beta))
      tau = sqrt(1/as.numeric(as.character(summarydata_OutlierIndep$n)))
      
      pheno.var1 = sum( (beta^2 - tau^2) * pow(beta) )
    
      Numdiscoveries = Numdiscoveries + Numdiscoveries1
      pheno.var = pheno.var + pheno.var1
    }
   
    GVpercentage = pheno.var/heritall*100
    return(list(Numdicoveries=Numdiscoveries, GVpercentage=GVpercentage, pheno.variance=pheno.var,herit=heritall))
  }
  
  
  library(MASS)
  logest = log(est)
  
  if(CI==TRUE){
    set.seed((seeds))
    alpha = (1-CI.coverage)/2
    logv = diag(1/est)%*%v%*%diag(1/est)
    estmat = exp(mvrnorm(nsim,mu=logest,Sigma=logv))
    
    tem = pp(est,n,summarydata_OutlierIndep,herit_OutlierIndep,gwas.significance,tol,M)
    tem1 = apply(estmat, 1, function(t) {pp(t,n,summarydata_OutlierIndep,herit_OutlierIndep,gwas.significance,tol,M)})
    
    pest = tem$Numdicoveries;
    gvest = tem$GVpercentage;
    pheno.variance = tem$pheno.variance;
    herit = tem$herit
    
    re = unlist(lapply(tem1,function(t) t[1]))
    rere = apply(matrix(re,ncol=1),1,function(t) rbinom(1,size=M,prob=t/M))
    regv = unlist(lapply(tem1,function(t)t[2]))
    
    return(list(heritability = herit, 
                Numdiscoveries = c(pest,quantile(rere,alpha),quantile(rere,1-alpha)),
                pheno.variance = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))*herit/100, 
                GVpercentage = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))
    ))
  }
  
  if(CI==FALSE){
    tem = pp(est,n,summarydata_OutlierIndep,herit_OutlierIndep,gwas.significance,tol,M)
    pest = tem$Numdicoveries;
    gvest = tem$GVpercentage;
    pheno.variance = tem$pheno.variance; 
    herit = tem$herit;
    
    return(list(heritability = herit, Numdiscoveries = pest, 
                pheno.variance=pheno.variance, GVpercentage = gvest))
  }
  
}