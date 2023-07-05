rm(list=ls())
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BRR2")
#setwd("~/DATA/thesis/Project_3_Grundtj?rn/analysis/genotype_changes/BLasso/CV_BRR2")

#load the libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(LDheatmap)
library(Hmisc)
library(GGally)
library(synbreed)
library(RColorBrewer)
library(rrBLUP)
library(asreml)
library(asremlPlus)
library(Matrix)
library(BGLR)
library(coda)
library(MCMCpack)
library(gdata)

#load("CV_BRR2_wood.RData")

#Setting path to working directory
path<-"~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BRR2"

#Save temp files in the same working directory
outpath<-"~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BRR2"


##### selecting the trait
head(phenot)
dim(phenot)
# trait for cv
EBV<-phenot[, "EBV_MOEd"]           
tt<-as.data.frame(phenot[, "MOEd"]) 
colnames(tt)<-c("MOEd")           

trait<-names(tt)
y <- tt[, trait]  
#traits <- c("MOEd")


samplesize <- 347
nfolds <- 2 # no of folds

results<-as.matrix(1:6,6, 1)

#***************************************************
# start first Loop 1
#**************************************************
for(j in 1:10){
  
  ## Matrix to store the summary of on each trait
  tabRes  <- data.frame(row.names=c('PredAbi', 'Accuracy', 'RankCor', 'MSE', 'Bias', 'Best10'))
  
  
  # 2nd Loop --> FOLDS          
  for(fold in 1:nfolds)
  {
    # Sample from y starting from observation 1, using samplesize
    whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
    yNa <- y              # make a copy of y. yNa and y are numeric  
    yNa[whichNa] <- NA    # set the values of sampled trees to NA
    
    
    # Run, See http://genomics.cimmyt.org/BGLR-extdoc.pdf
    fmR <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
                ETA = list(MRK = list(X = geno_centered, model = model)),
                nIter = 20000, burnIn = 1000, thin = 100,
                saveAt = "") 
    
    # Data frame, combine predicted (yHat and observed y)    
    df1 <- data.frame(fmR$yHat[whichNa], EBV[whichNa], y[whichNa])
    
    colnames(df1) <- c( 'GEBV', 'EBV', "Phenotype")
    
    # Calculate statistics
    PredAbi = round(cor(df1$Phenotype, df1$GEBV, use = "complete.obs"), 2) # PredAbi
    Accuracy=round(cor(df1$EBV,df1$GEBV, use = "complete.obs"), 2)
    RankCor = round(cor(df1$Phenotype, df1$GEBV, method = 'spearman', use = "complete.obs"), 2)
    Best10 <- round(mean(tail(sort(df1$Phenotype), n = ceiling(nrow(df1)*0.1))), 2) # Best 10%
    Bias   <- round(coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2], 2)
    MSE  <-  round(mean((df1$Phenotype - df1$GEBV)^2, na.rm=TRUE), 2)
    
    pred.sum <- rbind(PredAbi, Accuracy, RankCor,Best10, Bias, MSE)
    pred.sum <-round(pred.sum,2) 
    results<-cbind(results, pred.sum) 
  }
  
}

results

as.data.frame(apply(results[, -1], 1, mean))  # average of 10 times
sum_results<-round(cbind(as.data.frame(apply(results[, -1], 1, function(x){ mean(x,na.rm=TRUE) })),as.data.frame(apply(results[, -1], 1, function(x){ sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x))) })) ),2)
colnames(sum_results)<-c("Estimate", "SE")
sum_results

write.csv(sum_results, file=paste(trait,model,samplesize, ".csv" , sep=""), row.names=TRUE)
save.image('CV_BRR2_MOEd.RData')
