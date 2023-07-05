singlemapping <- function(X,y,covariate=NULL){
  
  #This function implements the single locus mapping based on linear regression. 
  
  #Author: Zitong Li, University of Helsinki, Email:zitong.li@helsinki.fi
  
  #Input values: 
  
  #-X: a data matrix containing genotypes, rows represent individuals, and column represent markers.
  
  #-y: a data vector containing phenotypes.
  
  #-covariate: environmental or other confounding factors such as age, sex, locations..., default is null
  
  # num_perm: the number of replicates in the permutation test, default is 10000
  
  # core: the number of cores to be used for parallel computation, to speed up the permutation test 
  
  # Output values:
  
  # pval: the original p-value of each SNP, calculated from the t-test
  
  # pval.corr: the adjusted p-value of each SNP by the permutation test
  
  # Coef: the estimated effect size of each SNP
  
  
  
  
  
  X <- as.matrix(X)
  
  y <- as.matrix(y)
  
  k <- dim(covariate)[2]
  
  p <- dim(X)[2]
  n <- dim(y)[1]
  T <- NULL
  Pval <- NULL
  Coef <- NULL
  for (j in 1:p){
    
    if (!is.null(covariate)){
      covariate <- as.matrix(covariate)
      fit <- summary(lm(y~covariate+X[,j]))
      Pval[j] <-  fit$coefficients[k+2,4]
      T[j] <-  abs(fit$coefficients[k+2,3])
      Coef[j] <- abs(fit$coefficients[k+2,1])
    }else{
      fit <- summary(lm(y~X[,j]))
      Pval[j] <-  fit$coefficients[2,4]
      T[j] <-  abs(fit$coefficients[2,3])
      Coef[j] <- abs(fit$coefficients[2,1])
    }
    
  }
  
 
  
  returnList <- list("pval" = Pval, "Coef"=Coef)
  return(returnList)
  
}


