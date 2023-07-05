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
library(asreml)
library(asremlPlus)
library(Matrix)
library(BGLR)
library(coda)
library(MCMCpack)
library(gdata)

# save.image("CV10rep_BL.RData")
#load("CV10rep_BL.RData")

#Setting path to working directory
path<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BL"

#Save temp files in the same working directory
outpath<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/CV_BL"

#Load the imputed genotypes from synbreed that I saved as pine_codeGeno.rda

load('pine_codeGeno.rda')


#Loading the phenotypes. I will use the same as in the GBlup cross validation.
phenot<-as.data.frame(gp.num$pheno)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(dbh30.1,dbh36.1,ht10.1,ht30.1))
colnames(phenotypes)<-c("dbh30", "dbh36", "ht10", "ht30")
head(phenotypes)
dim(phenotypes)
# phenotypes<-drop.levels(phenotypes)
# phenotypes<-phenotypes[!is.na(phenotypes$dbh30),]
# phenotypes<-phenotypes[!is.na(phenotypes$dbh36),]
# phenotypes<-phenotypes[!is.na(phenotypes$ht10),]
# phenotypes<-phenotypes[!is.na(phenotypes$ht30),]

# names(phenotypes)<-c('trait1', 'trait2', 'trait3', 'trait4')

# * Extract first trait 
trait <- colnames(phenotypes)[4]  # select first trait [1]
y <- phenotypes[, trait]

# create a vector of traits to analyze 
traits <- c("ht30")


genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])


genos2 <- as.matrix(genos)
is.matrix(genos2)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])

#### CROSS VALIDATION FOR BAYESIAN ##########
set.seed(110716)

### Run Bayesian Lasso for Cross Validation (CV)

# Choose parameters here
# Select model: from "BRR", "BL", "BayesA", "BayesB", "BayesC"
model <- "BL"
# 
# set.seed(12342) maybe not necessary

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
    y <- phenotypes$ht30	 
    
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
    
    # Make a summary of tabs for each folds of the trait
    tabs <- data.frame()
    tabs['PredAbi',trait] <- lambda(ftab$PredAbi)
    tabs['PA_min', trait] <- lammin(ftab$PredAbi)
    tabs['PA_max', trait] <-lammax(ftab$PredAbi)
    tabs['RankCor',trait] <- lambda(ftab$RankCor)
    tabs['MSE',trait] <- lambda(ftab$MSE)
    tabs['Bias',trait] <- lambda(ftab$Bias)
    tabs['Best10', trait] <- lambda(ftab$Best10)
    
    # See results 
    tabRes <- tabs
    results <- cbind(results, tabRes)
  }

}

results
str(results)
results[] <- lapply(results, function(x) as.numeric(as.character(x)))

results.average<-as.data.frame(apply(results[,-1],1,mean))

write.csv(results.average, file="CV10_BL_h30.csv", row.names = TRUE)


# Sample ebv and GEBV 
colnames(rawdf) <- c('fold', 'trait', 'GEBV', 'Phenotype')  
head(rawdf) 
tail(rawdf)
str(rawdf)
write.csv(rawdf, file="rawdf_CV10_BL_ht30.csv", row.names=TRUE)


# Print final SUMMARY stats      
tabRes  #only shows the last one

save.image("CV10rep_BL.RData")
