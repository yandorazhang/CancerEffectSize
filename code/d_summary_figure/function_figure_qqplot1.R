#-------------------------------------------------------------------
# Update Date: 04/21/2019
# Create Date: 04/21/2019
# Author: Yan (Dora) Zhang
# Goal: summary table results
#-------------------------------------------------------------------

figure_qq2com <- function(trait_name,trait_name_plot){
  LDwindow = 1; LDcutoff = 0.1
  
  output_path = paste0("../../genesis_result_new/",trait_name)
  if(!file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
    load(paste0(output_path,"/bestfit2.RData"))
    print(c(trait_name,"no"))
  }
  
  if(file.exists(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))){
    load(paste0(output_path,"/bestfit2_RemoveOutlier.RData"))
  }
  
  
  df = fit$qqplotdata$QQdata
  obs_lambda = fit$qqplotdata$observedlambda
  m.lambda = fit$qqplotdata$meanEXPlambda
  l.lambda = fit$qqplotdata$lowEXPlambda
  h.lambda = fit$qqplotdata$highEXPlambda
  
  # inx = seq(1,nrow(df),10)
  # df = df[inx,]
  plot(df$mean_log_exp_pvalues, df$log_obs_pvalues, xlab="",ylab="", type="l", xlim=c(0,10),ylim=c(0,10), main=trait_name_plot)
  polygon(c(df$lower,rev(df$upper)),c(df$log_obs_pvalues,rev(df$log_obs_pvalues)),col = "grey75", border = FALSE)
  # points(df$mean_log_exp_pvalues, df$log_obs_pvalues)
  abline(a=0,b=1 ,col = "gray50",lty=2)
  # #add red lines on borders of polygon
  # lines(y=df$log_obs_pvalues,col="red", df$upper,lty=2)
  # lines(y=df$log_obs_pvalues,col="red", df$lower,lty=2)
  
  mylabel =  bquote(italic(lambda[obs])== .(formatC(obs_lambda,format="f", digits = 3)))
  text(x = 8.5, y = 2.5, labels = mylabel)
  
  mylabel =  bquote(italic(lambda[fit])== .(formatC(m.lambda,format="f", digits = 3)))
  text(x = 8.5, y = 1.5, labels = mylabel)
}



figure_qq3com <- function(trait_name,trait_name_plot){
  LDwindow = 1; LDcutoff = 0.1
  
  output_path = paste0("../../genesis_result_new/",trait_name)
  if(!file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
    load(paste0(output_path,"/fit3.RData"))
    print(c(trait_name,"no"))
  }
  
  if(file.exists(paste0(output_path,"/fit3_RemoveOutlier.RData"))){
    load(paste0(output_path,"/fit3_RemoveOutlier.RData"))
    load(paste0("../../data_new/",trait_name,"/herit_OutlierIndep.RData"))
  }
  
  df = fit$qqplotdata$QQdata
  obs_lambda = fit$qqplotdata$observedlambda
  m.lambda = fit$qqplotdata$meanEXPlambda
  l.lambda = fit$qqplotdata$lowEXPlambda
  h.lambda = fit$qqplotdata$highEXPlambda
  
  
  # inx = seq(1,nrow(df),10)
  # df = df[inx,]
  plot(df$mean_log_exp_pvalues, df$log_obs_pvalues, type="l",xlab="",ylab="",  xlim=c(0,10),ylim=c(0,10),main=trait_name_plot)
  polygon(c(df$lower,rev(df$upper)),c(df$log_obs_pvalues,rev(df$log_obs_pvalues)),col = "grey75", border = FALSE)
  # points(df$mean_log_exp_pvalues, df$log_obs_pvalues)
  abline(a=0,b=1 ,col = "gray50",lty=2)
  # #add red lines on borders of polygon
  # lines(y=df$log_obs_pvalues,col="red", df$upper,lty=2)
  # lines(y=df$log_obs_pvalues,col="red", df$lower,lty=2)
  
  mylabel =  bquote(italic(lambda[obs])== .(formatC(obs_lambda,format="f", digits = 3)))
  text(x = 8.5, y = 2.5, labels = mylabel)
  
  mylabel =  bquote(italic(lambda[fit])== .(formatC(m.lambda,format="f", digits = 3)))
  text(x = 8.5, y = 1.5, labels = mylabel)
}


