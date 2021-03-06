#Wrapper function for the two factor ANOVA function
run_twofactor_cpp <- function(data,gpData,red_df){
  #Create design matrix for reduced model, i.e., no interaction effect
  colnames(gpData)[-c(1,2)] <- c("Factor1","Factor2")
  gpData <- cbind(gpData,y=1:nrow(gpData))
  
  #Create design matrix for the full model, i.e., all first order and interaction effects
  Xred <- unname(model.matrix(lm(y~Factor1+Factor2-1,data=gpData)))
  Xfull <- unname(model.matrix(lm(y~Factor1*Factor2-1,data=gpData)))
  
  #Run the two factor ANOVA model 
  #Edit 2/16: add "group_ids" so C++ knows which groups belong to which columns, helps with handling NaNs
  res <- two_factor_anova_cpp(data,Xfull,Xred,red_df, group_ids=as.numeric(factor(gpData$Group,levels=unique(gpData$Group))))
  
  #Get the unique group levels to translate ANOVA parameters into group means
  red_gpData <- dplyr::distinct(gpData,Group,.keep_all=TRUE)
  red_gpData <- dplyr::arrange(red_gpData,y)
  
  means <- res$par_estimates[,1:length(red_gpData$Group)]
  colnames(means) <- red_gpData$Group
  
  res$group_means <- means
  return(res)
}
