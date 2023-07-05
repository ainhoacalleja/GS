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

#load("CV_BRR2.RData")

#Setting path to working directory
path<-"~/DATA/thesis/Project_3_Grundtj?rn/analysis/genotype_changes/BLasso/CV_BRR2"

#Save temp files in the same working directory
outpath<-"~/DATA/thesis/Project_3_Grundtj?rn/analysis/genotype_changes/BLasso/CV_BRR2"

#save.image("CV_BRR2.RData")
#load("CV_BRR2.RData")

genos <- read.csv('genod.csv', row.names=1) #imputation performed by rrBLUP EM algorithm and in format [-1,0,1]
head(genos[,1:6])
dim(genos)

load('pine_codeGeno.rda')


#Loading the phenotypes. I will use the same as in the GBlup cross validation.
phenot<-as.data.frame(gp.num$pheno)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(dbh30.1,dbh36.1,ht10.1,ht30.1))

#change the name of the columns
colnames(phenotypes)<-c("dbh30", "dbh36", "ht10", "ht30")
head(phenotypes)
dim(phenotypes)

genos <- as.matrix(genos)
is.matrix(genos)
identical(sort(rownames(genos)), sort(rownames(phenotypes)))   #FALSE


 # To have same individuals in genotype and phenotype
new.data<-merge.data.frame(genos, phenotypes, by=0)  #merge phenotype and genotype by rows
head(new.data[, 1:6])
dim(new.data)
head(new.data[,8720:8724])  #Check las columns of the dataframe. Need to delete trait columns, i.e col from 8721:8724
new.data$ID=new.data$Row.names
dim(new.data)
head(new.data[,8720:8725])

new.data<-subset(new.data, select = -c(8721:8724)) #to delete trait columns
new.data[1:7, 8719:8721]
head(new.data[,1:6])
dim(new.data)

#Make first column of new.data as rownames for the data
new.data2 <- new.data[,-1]
rownames(new.data2) <- new.data[,1]
head(new.data2[,1:6])
str(new.data2)
is.data.frame(new.data2)
dim(new.data2)  #I still need to detelte the ID column

#Last step is just delete the ID column
new.data2$ID<-NULL
dim(new.data2)
head(new.data2[,1:6])

#Now I need the dataframe as matrix
#Faster and easiest way is to save as csv and read it again
write.csv(new.data2, file='small_genos.csv')

genotypes <- read.csv('small_genos.csv', row.names=1)
head(genotypes[,1:6])
genos2 <- as.matrix(genotypes)
is.matrix(genos2)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))   #Now we have the same individuals in both genotypes and phenotypes
save.image("CV_BRR2.RData")
phenotypes<-drop.levels(phenotypes)

geno_centered=genos2
dim(geno_centered)
str(geno_centered)
head(geno_centered[,1:6])
save.image("CV_BRR2.RData")

head(phenotypes)
# # * Extract first trait 
# trait <- colnames(phenotypes)[1]  # select first trait [1]
# y <- phenotypes[, trait]
# 
# # create a vector of traits to analyze 
# traits <- c("ht30")


## reading the datafile 
phenot2<-read.csv("growth_pheno.csv", row.names = 1)
head(phenot2)
dim(phenot2)

# * Extract first trait 
trait <- colnames(phenot2)[4]  # select first trait [1]
y <- phenot2[, trait]

# create a vector of traits to analyze 
traits <- c("ht30.1")

# trait for cv
EBV<-phenot2[, "EBV_h30"]           
#tt<-as.data.frame(phenot2[, "dbh30.1"]) 
#colnames(tt)<-c("dbh30.1")           

#trait<-names(tt)
#y <- tt[, trait]  
#traits <- c("dbh30.1")

###########################################
###########################################
### Run Bayesian Ridge Regression for Cross Validation (CV)

# Choose parameters here
# Select model: from "BRR", "BL", "BayesA", "BayesB", "BayesC"
model <- "BRR"

# sampling of population for cross validation 
samplesize <- 347
nfolds <- 2

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
      tab$PAcc <- round((cor(df1$Phenotype, df1$GEBV, use = "complete.obs")/sqrt(0.37)), 2)
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
write.csv(results.average, file="CV2_BRR_Ht30_EM.csv", row.names = TRUE)


# Sample ebv and GEBV 
colnames(rawdf) <- c('fold', 'trait', 'GEBV', 'EBV', 'Phenotype')  
head(rawdf) 
tail(rawdf)
str(rawdf)
write.csv(rawdf, file="TS_EM_CV2_BRR_Ht30.csv", row.names=TRUE)
save.image("CV_BRR2.RData")


