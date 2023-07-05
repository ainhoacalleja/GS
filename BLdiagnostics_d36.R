rm(list=ls())

setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso")

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

# save.image("BLdiag_d36.RData")
#load("BLdiag_d36.RData")

#Setting path to working directory
path<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso"

#Save temp files in the same working directory
outpath<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/outBLdiag_d36"

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


genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])
genos2 <- as.matrix(genos)
is.matrix(genos2)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])

#Select first trait (dbh36)
y <- phenotypes$trait2
head(y)
str(y)

str(geno_centered)

# * sampling of population for cross validation 
# samplesize <- 100
# fold <- 2 
# 
# # * sample vector y from observation 1, using samplesize 
# whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
# yNa <- y              # make a copy of y. yNa and y are numeric  
# yNa[whichNa] <- NA    # set the values of sampled trees to NA

###############################################################################
############################ MODEL DIAGNOSTICS ################################
set.seed(123)

# Bayesian LASSO Regression
# Number of samples
nIter = 20000
# Burn-in period for the Gibbs sampler
burnIn = 0
# Number of chains
n.chains=10

# r, proportion of phenotypic variance attributed to model residuals
r=0.5

# Prior hyperparameter values
# sigmaE2 (residual variance)
mode.sigE=r*var(y, na.rm = TRUE)
dfe=3
Se=mode.sigE*(dfe + 2)

# lambda
mode.sigL=(1-r)*var(y, na.rm=TRUE) # proportion attributed to genetics
lambda.hat=sqrt(2*mode.sigE/mode.sigL*sum(colMeans(geno_centered)^2, na.rm = TRUE))
delta.lambda=0.05                # rate
r.lambda=lambda.hat*delta.lambda # shape
shape=1.1
rate=2*Se/mode.sigL*sum(colMeans(geno_centered)^2, na.rm = TRUE)


# Set priors 
prior=list( varE=list(S0=Se,df0=dfe,value=runif(1,min=0,max=100)),
            lambda=list(type='random',
                        value=runif(1,min=0,max=100),
                        shape=shape, 
                        rate=rate) )

ETA<-list(MRK=list(X=geno_centered, model="BL", prior))

for (k in 1:n.chains){
  # Fit Bayesian LASSO Regression
  fmR<-BGLR(y=y,ETA=ETA,
            nIter=nIter,
            burnIn=burnIn,
            thin=1,
            rmExistingFiles = TRUE ,
            saveAt=paste(outpath,'/Gelman-Rubin_Chain_L_',k,sep=''))
}

save.image("BLdiag_d36.RData")

# Load MCMC draws of varE in matrix varE_L
varE_L=matrix(NA,ncol=n.chains,nrow=nIter)
for (k in 1:n.chains){
  varE_L[,k]=read.table( file=paste(outpath, '/Gelman-Rubin_Chain_L_',k,'varE.dat',sep='') )$V1
}

# Load MCMC draws of lambda in matrix Lambda
Lambda=matrix(NA,ncol=n.chains,nrow=nIter)
for (k in 1:n.chains){
  Lambda[,k]=read.table( file=paste(outpath, '/Gelman-Rubin_Chain_L_',k,'ETA_MRK_lambda.dat',sep='') )$V1
}




# Select parameter
draws=varE_L
#draws=Lambda

# Statistics for all chains
summary(mcmc(draws))

# MCMC object with all chains
idx=500:20000
THETA=mcmc.list(mcmc(draws[idx,1]),
                mcmc(draws[idx,2]),
                mcmc(draws[idx,3]),
                mcmc(draws[idx,4]),
                mcmc(draws[idx,5]))



# Gelman-Plot
# Create output file 
jpeg(file = "GelmanPlot_BL_dbh36.jpg", height = 10, width = 12, res = 1000, units = 'cm')
pdf(paste(outpath,"/GelmanPlot_dbh36_BL.pdf",sep=""), width=4, height=3)

par(mfrow=c(1,1))
par(mar = c(4.2, 4.2, 1, 1), mgp=c(3,1,0) )
gelman.plot(THETA,
            main='',  # Bayesian LASSO
            xlab='Iteration',
            cex.lab=1.2,
            cex=1.5,
            lwd=2)
dev.off()

# Gelman-Rubin statistic
gelman.rubin=gelman.diag(THETA)
gelman.rubin

# Geweke-Plot of one chain
geweke.diag(draws)
# geweke.plot(mcmc(as.matrix(draws[,5]))) 

save.image("BLdiag_d36.RData")

