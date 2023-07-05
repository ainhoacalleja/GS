# Remove everything in the working environment.
rm(list=ls())

setwd("~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso")

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


# save.image("bayes_compare.RData")
# load("bayes_compare.RData")

#Setting path to working directory
#path<-"~/Documents/PhD_UPSC/Project_3_GrundtjÃ¤rn/analysis/genotype_changes/BLasso"
path <- "~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso"

#Save temp files in the same working directory
outpath<-"~/DATA/thesis/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/out_bayes_compare"

#Load the imputed genotypes from synbreed that I saved as pine_codeGeno.rda

load('pine_codeGeno.rda')

# #BLRG requires matrix of markers
# genotype<-as.data.frame(gp.num$geno)
# dim(genotype)


#Loading the phenotypes. I will use the same as in the GBlup cross validation.
phenot<-as.data.frame(gp.num$pheno)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(dbh30.1,dbh36.1,ht10.1,ht30.1))
head(phenotypes)
dim(phenotypes)

names(phenotypes)<-c('trait1', 'trait2', 'trait3', 'trait4')

#Need same rows in genotypes and phenotypes

genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])
genos2 <- as.matrix(genos)
is.matrix(genos2)

identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])
save.image("bayes_compare.RData")

#Select first trait (dbh30)
y <- phenotypes$trait1
#is.data.frame((phenotypes))
head(y)
str(y)

str(geno_centered)

# * sampling of population for cross validation 
samplesize <- 100
fold <- 2 

# * sample vector y from observation 1, using samplesize 
whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
yNa <- y              # make a copy of y. yNa and y are numeric  
yNa[whichNa] <- NA    # set the values of sampled trees to NA

set.seed(12342)


###############################################################################
############################# COMPARING MODELS ################################
nIter = 20000; burnIn = 5000

#### Bayesian LASSO ####
ETA_BL<- list(MRK=list(X=geno_centered, model="BL",saveEffects=TRUE))
fmBL<-BGLR(y=yNa, ETA=ETA_BL, nIter=nIter, burnIn = burnIn, saveAt = "BL_")
# B_BL=readBinMat("BL_ETA_MRK_b.bin")
# BL_var <- getVariances(B=B_BL, X=geno_centered, sets=sample(1:20,size=1279,replace=T))
# head(BL_var)
# plot(BL_var[,"total"])
summary(fmBL)
str(fmBL)
predict(fmBL)
fmBL$varE
str(ETA_BL)

#### Bayesian Ridge Regression ####
ETA_BRR<-list(MRK=list(X=geno_centered, model="BRR",saveEffects=TRUE))
fmBRR <- BGLR (y=yNa, ETA=ETA_BRR, nIter = nIter, burnIn = burnIn, saveAt = "BRR_")
summary(fmBRR)
str(fmBRR)
varU=scan('BRR_ETA_MRK_varB.dat')
varE=scan('varE.dat')
mean(varU)
varE
h2=varU/(varU+varE)

#### Bayes A (Scaled-t prior) ####
ETA$MRK$model<-"BayesA" 
fmBA<-BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,saveAt="BA_") 

#### Bayes B (point of mass at zero + scaled-t slab)  ####
ETA$MRK$model<-"BayesB" 
fmBB<-BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,saveAt="BB_") 

#### Bayes C (point of mass at zero + scaled-t slab)  ####
ETA$MRK$model<-"BayesC" 
fmBC<-BGLR(y=yNa,ETA=ETA,nIter=nIter,burnIn=burnIn,saveAt="BC_") 

####   Calculate correlations  for each model
r_BL  <- cor(y[whichNa], fmBL$yHat[whichNa] )
r_BRR <- cor(y[whichNa],fmBRR$yHat[whichNa] )
r_BA  <- cor(y[whichNa], fmBA$yHat[whichNa] ) 
r_BB  <- cor(y[whichNa], fmBB$yHat[whichNa] ) 
r_BC  <- cor(y[whichNa], fmBC$yHat[whichNa] ) 


# organize in a table 
table <- data.frame(rbind( r_BL,          r_BRR,         r_BA,         r_BB,  r_BC) , 
                    rbind( fmBL$varE,     fmBRR$varE,    fmBA$varE,    fmBB$varE, fmBC$varE) ,
                    rbind( fmBL$fit$pD,   fmBRR$fit$pD,  fmBA$fit$pD,  fmBB$fit$pD, fmBC$fit$pD)  ,
                    rbind( fmBL$fit$DIC,  fmBRR$fit$DIC, fmBA$fit$DIC, fmBB$fit$DIC, fmBC$fit$DIC)  ) 

colnames(table) <- c( 'PredAbi', "varE", "pD", "DIC")
rownames(table) <- c( 'BLasso', "BRidge", "BayesA", "BayesB", "BayesC")
round(table, 3)
save.image("bayes_compare.RData")

# Estimated effects   #this does not work
pdf("~/Documents/PhD_UPSC/Project_3_GrundtjÃ¤rn/analysis/genotype_changes/BLasso/out_bayes_compare/Pred_BRR_BB_graphs.pdf", width = 4, height = 4)

par(mar=c(4,4.5,0,0)+0.2)#sets margins of plotting area

p<-ncol(geno_centered)
#   tmp<-range(abs(b0))
plot(numeric()~numeric(),ylim=c(-0.05,0.05),xlim=c(1,p),
     ylab=expression(paste("|",beta[j],"|")),
     xlab="Marker Possition (order)",yaxt='n')
axis(2,at=seq(-0.05,0.05,len=5),lab=c(".05",".025","0",".025",".05"))
# abline(v=whichQTL,lty=2,col=4)
# points(x=whichQTL,y=abs(b0[whichQTL]),pch=19,col=4)

points(x=1:p,y=abs(fmBRR$ETA$MRK$b),col=1,cex=0.6, type="o")
points(x=1:p,y=-abs(fmBB$ETA$MRK$b),col=2,cex=0.6,type="o")
text(x=300,y=-0.05,label="Bayes-B", col ="red")
text(x=300,y=0.05,label="Bayes-RR")

# dev.off()
save.image("bayes_compare.RData")
