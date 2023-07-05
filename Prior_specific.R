rm(list=ls())

setwdsetwd("~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso")

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

# save.image("Lasso_prior.RData")
#load("Lasso_prior.RData")

#Setting path to working directory
path<-"~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso"

#Save temp files in the same working directory
outpath<-"~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/out_bayes_synbreed"

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
save.image("Lasso_prior.RData")


genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])
genos2 <- as.matrix(genos)
is.matrix(genos2)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])
save.image("Lasso_prior.RData")

#Select first trait (dbh30)
y <- phenotypes$trait1
head(y)
str(y)

# # * sampling of population for cross validation 
# samplesize <- 100
# fold <- 2 

# # * sample vector y from observation 1, using samplesize 
# whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
# yNa <- y              # make a copy of y. yNa and y are numeric  
# yNa[whichNa] <- NA    # set the values of sampled trees to NA


###############################################################################
####################### PRIOR SPECIFICATION ##########################
###############################################################################
nIter = 20000
# Burn-in period for the Gibbs sampler
burnIn = 5000

# r, proportion of phenotypic variance attributed to model residuals
R2=0.5

# Prior hyperparameter values
# sigmaE2 (residual variance)
mode.sigE=R2*var(y, na.rm=TRUE)
dfe=5
Se=mode.sigE*(dfe + 2)

# lambda
mode.sigL=(1-R2)*var(y, na.rm=TRUE) # variance explained by markers
rate=2*Se/mode.sigL*sum(colMeans(geno_centered)^2)
shape=1.1
# Set priors
prior=list( varE=list(S0=Se,df0=dfe),
            lambda=list(type='random',
                        shape=shape,
                        rate=rate) )

ETA<-list(MRK=list(X=geno_centered, model="BL", prior))

fmBL <- BGLR(y = y,  ETA = ETA,  nIter=nIter,  burnIn=burnIn, saveAt=outpath )

str(fmBL)
summary(fmBL)
save.image("Lasso_prior.RData")
# ###############################################################################
# ############################ PRIOR INFLUENCE ################################
# ###############################################################################
# 
# Bayesian LASSO
# propE, proportion of phenotypic variance attributed to model residuals
# Number of samples
nIter = 20000
# Burn-in period for the Gibbs sampler
burnIn = 5000
# propE, proportion of phenotypic variance attributed to model residuals
propE=c(0.8,0.6,0.4)
for (r in propE ){
  # Prior hyperparameter values
  # sigmaE2
  mode.sigE=r*var(y, na.rm=TRUE)
  dfe=3
  Se=mode.sigE*(dfe + 2)

  # lambda
  mode.sigL=(1-r)*var(y, na.rm=TRUE)
  lambda.hat=sqrt(2*mode.sigE/mode.sigL*sum(colMeans(geno_centered)^2))
  delta.lambda=0.05                # rate
  r.lambda=lambda.hat*delta.lambda # shape

  # Set priors
  prior=list( varE=list(S=Se,df=dfe),
              lambda=list(type='random',
                          value=lambda.hat,
                          shape=r.lambda,
                          rate=delta.lambda) )

  ETA<-list(MRK=list(X=geno_centered, model="BL", prior))
  # Fit Bayesian LASSO Regression
  fmR<-BGLR(y=y,ETA=ETA,
            nIter=nIter,
            burnIn=burnIn,
            thin=1,
            saveAt=paste(outpath,'/LASSO_r_',r,'_',sep=''))
}

# Load MCMC draws of varE in matrix varE_L
k=0
varE_L=matrix(NA,ncol=length(propE),nrow=nIter)
for (r in propE){
  k=k+1
  varE_L[,k]=read.table(file=paste(outpath,'/LASSO_r_',r,'_varE.dat',sep=''),sep=' ')$V1
}

# BL posterior means and SEs
round( colMeans(varE_L[(burnIn+1):nIter,]) , 2 )
round( apply(varE_L[(burnIn+1):nIter,],2,sd) , 2 )

colMeans(Lambda[(burnIn+1):nIter,])
apply(Lambda[(burnIn+1):nIter,],2,sd)


# Prior densities for varE
k=0
N=2^9
priorE=matrix(NA,ncol=length(propE),nrow=N)
for(r in propE){
  k=k+1
  mode.sigE=r*var(y, na.rm=TRUE)
  dfe=3
  Se=mode.sigE*(dfe + 2)
  xE=seq(from=0.01,to=0.8,length.out=N)
  priorE[,k]=dinvgamma(x=xE,shape=dfe/2, scale=Se/2)
}

# Posterior densities for varE
for (k in 1:length(propE)){
  namx=paste('dLx_',k,sep='')
  namy=paste('dLy_',k,sep='')
  dL=density(varE_L[(burnIn+1):nIter,k],n=N)
  assign(namx,dL$x)
  assign(namy,dL$y)
}

##### FIGURE 12.6_Priors Effects
# Put the output to a file rather than the standard output
# jpeg(file = "BL_SWEEP.jpg", height = 8, width = 10, res = 1000, units = 'in')
#pdf(paste(image.path,"/Fig12-6_PriorsInfluence.pdf",sep=""), width = 6, height = 4)

# Plot
par(mfrow=c(1,1))
par(mar = c(4.3, 4.3, 0.9, 0.9), mgp=c(3,1,0) )

plot(xE,priorE[,1]/max(priorE[,1]),
     xlab=expression(sigma[e]^2),ylab='Scaled density',
     main=' ', xlim = c(0,1),
     type='l',lty=2,lwd=2, cex.lab=1.2,col=2)

for (k in 2:length(propE)){
  lines(xE,priorE[,k]/max(priorE[,k]),type='l',lty=2,lwd=2,col=k+1) }

for (k in 1:length(propE)){
  xx=get(paste('dLx_',k,sep=''))
  yy=get(paste('dLy_',k,sep=''))
  lines(xx,yy/max(yy),type='l',lwd=2,col=k+1)
}

legend('bottomright',lty=c(1,1,1), lwd=2, col=c(2,3,4),
       legend=c('Posterior (Prior: 20% Genetics)',
                'Posterior (Prior: 40% Genetics)',
                'Posterior (Prior: 60% Genetics)'),
       bty='n',cex=0.9)

#dev.off()
save.image("Lasso_prior.RData")
