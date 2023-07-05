rm(list=ls())

setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BL")

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
#library(asreml)
library(asremlPlus)
library(Matrix)
library(BGLR)
library(coda)
library(MCMCpack)
library(gdata)

# save.image("CV10_BL_wood.RData")
#load("CV10_BL_wood.RData")

#Setting path to working directory
path<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BL"

#Save temp files in the same working directory
outpath<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BL"

#Loading genotypic data
genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])
dim(genos)

genos2 <- as.matrix(genos)
is.matrix(genos2)

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])

#Loading the wood phenotypes. 
phenot<-read.csv("wood.csv", row.names = 1)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(ad_den,MFA,MOEd,MOEs))
colnames(phenotypes)<-c("ad_den", "MFA", "MOEd", "MOEs")
head(phenotypes)
dim(phenotypes)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

# * Extract first trait 
trait <- colnames(phenotypes)[3]  # select first trait [1]
y <- phenotypes[, trait]

# create a vector of traits to analyze 
traits <- c("MOEd")

###########################################
### Run Bayesian Lasso for 10 fold Cross Validation (CV)
phenot2<-read.csv(file = "growth_pheno.csv", row.names = 1)
head(phenot2)
# trait for density
EBV<-phenot2[, "EBV_h30"]           
tt<-as.data.frame(phenot2[, "ht30"]) 
colnames(tt)<-c("ht30")           

trait<-names(tt)
y <- tt[, trait]  
traits <- c("ht30")


# Choose parameters here
# Select model: from "BRR", "BL", "BayesA", "BayesB", "BayesC"
model <- "BL"
# 
set.seed(12342) #maybe not necessary

# sampling of population for cross validation 
samplesize <- 139 #69, 139, 208, 277, 347
nfolds <- 5

## Matrix to store the summary of on each trait
tabRes  <- data.frame(row.names=c('PredAbi','PredAbi_SE','PA_min', 'PA_max', 'Accuracy','Acc_SE','PAcc','PAcc_SE', 'RankCor', 'MSE', 'Bias', 'Best10'))
head(tabRes)

## data frame for sampled and predicted data  
rawdf <- data.frame() 

# data frame for summary statistics 
tabdf <- data.frame()

results<-as.matrix(1:12, 12, 1)

for (i in 1:10) {
  # 1st Loop --> TRAITS 
  for (trait in traits){
    # Phenotype of a single trait
    y <- phenot2[, trait]	 
    
    # statistics for all folds of a trait are saved in ftab 
    ftab <- data.frame()
    
    # 2nd Loop --> FOLDS  	
    for(fold in 1:nfolds)
    {
      # Sample from vector y starting from observation 1, using samplesize
      whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
      yNa <- y              # make a copy of y. yNa and y are numeric  
      yNa[whichNa] <- NA    # set the values of sampled trees to NA
      
      # Run, 
      fmR <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
                  ETA = list(MRK = list(X = geno_centered, model = model)),
                  nIter = 20000, burnIn = 1000, thin = 100,
                  saveAt = "")  
      
      # Data frame (sample size x 2): combine predicted (yHat and observed y)    
      df1 <- data.frame(rep(fold, samplesize), rep(trait, samplesize),
                        fmR$yHat[whichNa],EBV[whichNa],y[whichNa])
      
      colnames(df1) <- c('fold', 'trait', 'GEBV','EBV', 'Phenotype')
      
      # Calculate 5 statistics and put them in matrix tab, previously created 
      tab = list()
      tab$fold <- fold
      tab$trait <- trait
      tab$PredAbi <- round(cor(df1$Phenotype, df1$GEBV, use = "complete.obs"), 2) # PredAbi
      tab$Accuracy <- round(cor(df1$EBV, df1$GEBV, use = "complete.obs"), 2) #Accuracy
      tab$PAcc <- round((cor(df1$Phenotype, df1$GEBV, use = "complete.obs")/sqrt(0.36)), 2)
      tab$RankCor <- round(cor(df1$Phenotype, df1$GEBV, method = 'spearman', use = "complete.obs"), 2) # RankCor
      tab$MSE <- round(mean((df1$Phenotype - df1$GEBV)^2, na.rm=TRUE), 2) # MSE 
      tab$Bias <- round(coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2], 2) # Bias
      tab$Best10 <- round(mean(tail(sort(df1$Phenotype), n = ceiling(nrow(df1)*0.1))), 2) # Best 10%
      
      # Save df of GEBV and phenotype
      rawdf = rbind(rawdf, df1)  
      # Save statistics
      ftab <- rbind(ftab, tab)
    }
    
    # Save tab results in a df 
    tabdf <- rbind(tabdf, ftab)
    
    # Function to create a mean for the folds
    lambda  <- function(x) paste(round(mean(x), 2))
    lambdaSE<-function(x) paste(round((sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))),2))
    
    # Function to create a max for the folds
    lammin  <- function(x) paste(round( min(x), 2))
    lammax <- function (x) paste(round( max(x), 2))
    
    # # Function to create a summary for the folds
    # lambda  <- function(x) paste(round(mean(x), 2), " (", 
    #                              round( min(x), 2), "-", 
    #                              round( max(x), 2), ")", sep = "")
    
    # Make a summary of tabs for each folds of the trait
    tabs <- data.frame()
    tabs['PredAbi', trait] <- lambda(ftab$PredAbi)
    tabs['PredAbi_SE', trait] <- lambdaSE(ftab$PredAbi)
    tabs['PA_min', trait] <- lammin(ftab$PredAbi)
    tabs['PA_max', trait] <-lammax(ftab$PredAbi)
    tabs['Accuracy', trait] <- lambda(ftab$Accuracy)
    tabs['Acc_SE', trait] <- lambdaSE(ftab$Accuracy)
    tabs['PAcc', trait] <- lambda(ftab$PAcc)
    tabs['PAcc_SE', trait] <- lambdaSE(ftab$PAcc)
    tabs['RankCor', trait] <- lambda(ftab$RankCor)
    tabs['MSE', trait] <- lambda(ftab$MSE)
    tabs['Bias', trait] <- lambda(ftab$Bias)
    tabs['Best10',trait] <- lambda(ftab$Best10)
    
    # See results 
    tabRes[, trait] <- tabs[, trait]
    results <- cbind(results, tabRes)
  }
}

results
str(results)
results[] <- lapply(results, function(x) as.numeric(as.character(x)))

results.average<-as.data.frame(apply(results[,-1],1,mean))
results.average
write.csv(results.average, file="CV8_BL_Ht30.csv", row.names = TRUE)


# Sample ebv and GEBV 
colnames(rawdf) <- c('fold', 'trait', 'GEBV', 'EBV', 'Phenotype')  
head(rawdf) 
tail(rawdf)
str(rawdf)
write.csv(rawdf, file="TS_CV8_BL_Ht30.csv", row.names=TRUE)
save.image("CV_BRR_wood.RData")


# # 1st Loop --> TRAITS 
# for (trait in traits){
#   # Phenotype of a single trait
#   y <- tt[, trait]	 
#   
#   # statistics for all folds of a trait are saved in ftab 
#   ftab <- data.frame()
#   
#   # 2nd Loop --> FOLDS  	
#   for(fold in 1:nfolds)
#   {
#     # Sample from vector y starting from observation 1, using samplesize
#     whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
#     yNa <- y              # make a copy of y. yNa and y are numeric  
#     yNa[whichNa] <- NA    # set the values of sampled trees to NA
#     
#     # Run, See http://genomics.cimmyt.org/BGLR-extdoc.pdf
#     fmR <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
#                 ETA = list(MRK = list(X = geno_centered, model = model)),
#                 nIter = 20000, burnIn = 1000, thin = 100,
#                 saveAt = "")  
#     
#     # Data frame (sample size x 2): combine predicted (yHat and observed y)    
#     df1 <- data.frame(rep(fold, samplesize), rep(trait, samplesize),
#                       fmR$yHat[whichNa],EBV[whichNa], y[whichNa])
#     
#     colnames(df1) <- c('fold', 'trait', 'GEBV', 'EBV', 'Phenotype')
#     
#     # Calculate 5 statistics and put them in matrix tab, previously created 
#     tab = list()
#     tab$fold <- fold
#     tab$trait <- trait
#     tab$PredAbi <- round(cor(df1$Phenotype, df1$GEBV, use = "complete.obs"), 2) # PredAbi
#     tab$Accuracy <- round(cor(df1$EBV, df1$GEBV, use = "complete.obs"), 2) #Accuracy
#     tab$RankCor <- round(cor(df1$Phenotype, df1$GEBV, method = 'spearman', use = "complete.obs"), 2) # RankCor
#     tab$MSE <- round(mean((df1$Phenotype - df1$GEBV)^2, na.rm=TRUE), 2) # MSE 
#     tab$Bias <- round(coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2], 2) # Bias
#     tab$Best10 <- round(mean(tail(sort(df1$Phenotype), n = ceiling(nrow(df1)*0.1))), 2) # Best 10%
#     
#     # Save df of GEBV and phenotype
#     rawdf = rbind(rawdf, df1)  
#     # Save statistics
#     ftab <- rbind(ftab, tab)
#   }
#   
#   # Save tab results in a df 
#   tabdf <- rbind(tabdf, ftab)
#   
#   # Function to create a summary for the folds
#   lambda  <- function(x) paste(round(mean(x), 2), " (", 
#                                round( min(x), 2), "-", 
#                                round( max(x), 2), ")", 
#                                round((sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))),2) ,sep = "")
#   # Make a summary of tabs for each folds of the trait
#   tabs <- data.frame()
#   tabs['PredAbi', trait] <- lambda(ftab$PredAbi)
#   tabs['Accuracy', trait] <- lambda(ftab$Accuracy)
#   tabs['RankCor', trait] <- lambda(ftab$RankCor)
#   tabs['MSE', trait] <- lambda(ftab$MSE)
#   tabs['Bias', trait] <- lambda(ftab$Bias)
#   tabs['Best10',trait] <- lambda(ftab$Best10)
#   
#   # See results 
#   tabRes[, trait] <- tabs[, trait]
# }
# 
# # Sample ebv and GEBV 
# colnames(rawdf) <- c('fold', 'trait', 'GEBV', 'EBV', 'Phenotype')  
# head(rawdf) 
# tail(rawdf)
# str(rawdf)
# 
# 
# # Print final SUMMARY stats      
# tabRes
# write.csv(tabRes, file="tabres_BL_MOEd.csv", row.names = TRUE)
# write.csv(rawdf,file = "10fold_BL_MOEd.csv", row.names = TRUE )

save.image("CV10_BL_wood.RData")





###########################################
##### CV for Different TS AND VS
### Run Bayesian Lasso for Cross Validation (CV)

# Choose parameters here
# Select model: from "BRR", "BL", "BayesA", "BayesB", "BayesC"
model <- "BL"

# sampling of population for cross validation 
samplesize <- 69
nfolds <- 10

## Matrix to store the summary of on each trait
tabRes  <- data.frame(row.names=c('PredAbi', 'PA_min', 'PA_max', 'RankCor', 'MSE', 'Bias', 'Best10'))
head(tabRes)

## data frame for sampled and predicted data  
rawdf <- data.frame() 

# data frame for summary statistics 
tabdf <- data.frame()

results<-as.matrix(1:7, 7, 1)

for (i in 1:10) {
  # 1st Loop --> TRAITS 
  for (trait in traits){
    # Phenotype of a single trait
    y <- phenotypes[, trait]	 
    
    # statistics for all folds of a trait are saved in ftab 
    ftab <- data.frame()
    
    # 2nd Loop --> FOLDS  	
    for(fold in 1:nfolds)
    {
      # Sample from vector y starting from observation 1, using samplesize
      whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
      yNa <- y              # make a copy of y. yNa and y are numeric  
      yNa[whichNa] <- NA    # set the values of sampled trees to NA
      
      # Run, See http://genomics.cimmyt.org/BGLR-extdoc.pdf
      fmR <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
                  ETA = list(MRK = list(X = geno_centered, model = model)),
                  nIter = 20000, burnIn = 1000, thin = 100,
                  saveAt = "")  
      
      # Data frame (sample size x 2): combine predicted (yHat and observed y)    
      df1 <- data.frame(rep(fold, samplesize), rep(trait, samplesize),
                        fmR$yHat[whichNa], y[whichNa])
      
      colnames(df1) <- c('fold', 'trait', 'GEBV', 'Phenotype')
      
      # Calculate 5 statistics and put them in matrix tab, previously created 
      tab = list()
      tab$fold <- fold
      tab$trait <- trait
      tab$PredAbi <- round(cor(df1$Phenotype, df1$GEBV, use = "complete.obs"), 2) # PredAbi
      tab$RankCor <- round(cor(df1$Phenotype, df1$GEBV, method = 'spearman', use = "complete.obs"), 2) # RankCor
      tab$MSE <- round(mean((df1$Phenotype - df1$GEBV)^2, na.rm=TRUE), 2) # MSE 
      tab$Bias <- round(coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2], 2) # Bias
      tab$Best10 <- round(mean(tail(sort(df1$Phenotype), n = ceiling(nrow(df1)*0.1))), 2) # Best 10%
      
      # Save df of GEBV and phenotype
      rawdf = rbind(rawdf, df1)  
      # Save statistics
      ftab <- rbind(ftab, tab)
    }
    
    # Save tab results in a df 
    tabdf <- rbind(tabdf, ftab)
    
    # Function to create a mean for the folds
    lambda  <- function(x) paste(round(mean(x), 2))
    # Function to create a max for the folds
    lammin  <- function(x) paste(round( min(x), 2))
    lammax <- function (x) paste(round( max(x), 2))
    
    # # Function to create a summary for the folds
    # lambda  <- function(x) paste(round(mean(x), 2), " (", 
    #                              round( min(x), 2), "-", 
    #                              round( max(x), 2), ")", sep = "")
    # Make a summary of tabs for each folds of the trait
    tabs <- data.frame()
    tabs['PredAbi', trait] <- lambda(ftab$PredAbi)
    tabs['PA_min', trait] <- lammin(ftab$PredAbi)
    tabs['PA_max', trait] <-lammax(ftab$PredAbi)
    tabs['RankCor', trait] <- lambda(ftab$RankCor)
    tabs['MSE', trait] <- lambda(ftab$MSE)
    tabs['Bias', trait] <- lambda(ftab$Bias)
    tabs['Best10',trait] <- lambda(ftab$Best10)
    
    # See results 
    tabRes[, trait] <- tabs[, trait]
    results <- cbind(results, tabRes)
  }
}

results
str(results)
results[] <- lapply(results, function(x) as.numeric(as.character(x)))

results.average<-as.data.frame(apply(results[,-1],1,mean))

write.csv(results.average, file="CV8_BL_MOEd.csv", row.names = TRUE)


# Sample ebv and GEBV 
colnames(rawdf) <- c('fold', 'trait', 'GEBV', 'Phenotype')  
head(rawdf) 
tail(rawdf)
str(rawdf)
write.csv(rawdf, file="TS_BL_CV8_MOEd.csv", row.names=TRUE)


# Print final SUMMARY stats      
#tabRes  #only shows the last one


save.image("CV_BL_wood.RData")
