# Remove everything in the working environment.
rm(list=ls())

setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/Bayesian_rrblup")

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

#Setting path to working directory
path<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/bayesian_rrblup"

#Save temp files in the same working directory
outpath<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/bayesian_rrblup"

# save.image("BL_rrblup.RData")
#load("BL_rrblup.RData")

#Load the imputed genotypes from synbreed that I saved as pine_codeGeno.rda

load('pine_codeGeno.rda')


#Loading the phenotypes. I will use the same as in the GBlup cross validation.
phenot<-as.data.frame(gp.num$pheno)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(dbh30.1,dbh36.1,ht10.1,ht30.1))
head(phenotypes)
dim(phenotypes)

names(phenotypes)<-c('trait1', 'trait2', 'trait3', 'trait4')

#Loading the genotypes
genos <- read.csv('small_genos.csv', row.names=1)
head(genos[,1:6])
genos2 <- as.matrix(genos)
is.matrix(genos2)


#Select first trait (dbh30)
y <- phenotypes$trait1
head(y)
str(y)


# * sampling of population for cross validation 
samplesize <- 100
fold <- 2 

# * sample vector y from observation 1, using samplesize 
whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
yNa <- y              # make a copy of y. yNa and y are numeric  
yNa[whichNa] <- NA    # set the values of sampled trees to NA

#############################################
###     Bayesian LASSO 
#############################################
set.seed(123)
d30_BL <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
               ETA = list(MRK = list(X = genos2, model = "BL")),
               nIter = 20000, burnIn = 1000, thin = 100,
               saveAt = outpath)

summary(d30_BL)
str(d30_BL)

####################################################
####  Calculate statistics  for the LASSO results

df1 <- data.frame(d30_BL$yHat[whichNa], y[whichNa])

PredAbi = round( cor(df1[, 2], df1[, 1], use = "complete.obs"), 2)                              
Rank = round( cor(df1[, 2], df1[, 1], method = 'spearman', use = "complete.obs"), 2)         
MSE =  round( mean((df1[, 2] - df1[,1])^2, na.rm=TRUE), 2)                          
Bias  = round( coef(lm(y[whichNa] ~ d30_BL$yHat[whichNa]))[2], 2)        
Best10 =   round( mean(tail(sort(df1[, 2]), n = ceiling(nrow(df1)*0.1))), 2)   
# The MSE of training set 
# MSE.trn<-mean((fmR$yHat[-whichNa]-y[-whichNa])^2,na.rm=T)

# organize in a table 
table <- cbind(PredAbi, Rank, MSE, Bias, Best10)    
rownames(table) <- c('d30')
table 

save.image("BL_rrblup.RData")
