# Remove everything in the working environment.
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


# save.image("bayes_synbreed.RData")
#load("bayes_synbreed.RData")

#Setting path to working directory
path<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso"

#Save temp files in the same working directory
outpath<-"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/BLasso/out_bayes_synbreed"

#Load the imputed genotypes from synbreed that I saved as pine_codeGeno.rda

load('pine_codeGeno.rda')

#BLRG requires matrix of markers
genotype<-as.data.frame(gp.num$geno)
dim(genotype)


#Loading the phenotypes. I will use the same as in the GBlup cross validation.
phenot<-as.data.frame(gp.num$pheno)
head(phenot)
dim(phenot)
phenotypes<-subset(phenot,select = c(dbh30.1,dbh36.1,ht10.1,ht30.1))
head(phenotypes)
dim(phenotypes)

names(phenotypes)<-c('trait1', 'trait2', 'trait3', 'trait4')
save.image("bayes_synbreed.RData")

#Need same rows in genotypes and phenotypes
identical(sort(rownames(genotype)), sort(rownames(phenotypes)))
new.data<-merge.data.frame(genotype, phenotypes, by=0)
new.data$ID=new.data$Row.names

head(genotype[,1:6])

dim(new.data)
new.data[1:7, 8706:8712]
new.data[1:7, 1:6]

new.data<-subset(new.data, select = -c(8708:8711)) #to delete trait columns
new.data[1:7, 8706:8708]
dim(new.data)

#Make first column of new.data as rownames for the data
new.data2 <- new.data[,-1]
rownames(new.data2) <- new.data[,1]
head(new.data2[,1:6])
str(new.data2)
is.data.frame(new.data2)
dim(new.data2)

#Last step is just delete the ID column
new.data2$ID<-NULL
dim(new.data2)
head(new.data2[,1:6])

#Now I need the dataframe as matrix
#Faster and easiest way is to save as csv and read it again
write.csv(new.data2, file='new.genotype.csv')

genos <- read.csv('new.genotype.csv', row.names=1)
head(genos[,1:6])
genos2 <- as.matrix(genos)
is.matrix(genos2)

#identical(sort(rownames(genos2)), sort(rownames(phenotypes)))  #TRUE

#markers should be centered on zero for easier calculation, i.e., coded as {-1, 0, 1}
geno_centered<- genos2-1
head(geno_centered[,1:6])
save.image("bayes_synbreed.RData")

#Select first trait (ht30)
y <- phenotypes$trait4
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

#############################################
###     Bayesian LASSO 
#############################################
set.seed(123)
h30_BL <- BGLR(y  = yNa, response_type = "gaussian", a = NULL, b = NULL,
            ETA = list(MRK = list(X = geno_centered, model = "BL")),
            nIter = 20000, burnIn = 1000, thin = 100,
            saveAt = outpath)

summary(h30_BL)
str(h30_BL)
save.image("bayes_synbreed.RData")

####################################################
####  Calculate statistics  for the LASSO results

df1 <- data.frame(h30_BL$yHat[whichNa], y[whichNa])

PredAbi = round( cor(df1[, 2], df1[, 1], use = "complete.obs"), 2)                              
Rank = round( cor(df1[, 2], df1[, 1], method = 'spearman', use = "complete.obs"), 2)         
MSE =  round( mean((df1[, 2] - df1[,1])^2, na.rm=TRUE), 2)                          
Bias  = round( coef(lm(y[whichNa] ~ h30_BL$yHat[whichNa]))[2], 2)        
Best10 =   round( mean(tail(sort(df1[, 2]), n = ceiling(nrow(df1)*0.1))), 2)   
# The MSE of training set 
# MSE.trn<-mean((fmR$yHat[-whichNa]-y[-whichNa])^2,na.rm=T)

# organize in a table 
table <- cbind(PredAbi, Rank, MSE, Bias, Best10)    
rownames(table) <- c('ht30')
table 
save.image("bayes_synbreed.RData")

###############################################################################
####################### GRAPHS
###############################################################################
### * 1- phenotype (y=phenotype) and direct genetic values (yHat)
plot( h30_BL$yHat ~ y,   xlab = "Phenotype", ylab = "GEBV", main='Height 30',
      col = densCols(y), pch = 20, cex = 1)

# add GEBV 
points(x = y[whichNa], y = h30_BL$yHat[whichNa], 
       col = 2, bg = 'red', cex = .9, pch = 21)
# add correlation to the legend 
txt <- paste("r_vset = ", round(cor(h30_BL$yHat[whichNa], y[whichNa], 
                                    use = 'complete.obs'), 2), sep = "" ) 
legend("topleft", txt, col = 'red3', bty = 'n')


### 2 Predictions
yHat<-h30_BL$yHat
tmp<-range(c(y,yHat), na.rm=TRUE, finite=FALSE)
plot(yHat~y, xlab='Observed', ylab='Predicted',
     col=2, xlim=tmp,ylim=tmp, main= 'Height 30')
abline(a=0,b=1,col=4,lwd=2) 

#############################################
### * 2 Estimated Marker Effects & posterior SDs
bHat<- h30_BL$ETA$MRK$b
SD.bHat<- h30_BL$ETA$MRK$SD.b
plot(bHat^1, ylab='Estimated Marker Effect', 
     ylim=c(0,0.06),
     type='o',cex=.5,col=4,main='Marker Effects')

points(SD.bHat^1, ylab='Estimated SD',
       type='o',cex=.5,col=6,main='Marker Effects')


########## Visualize marker effects
# Create data frame 
meff <- data.frame(bHat)
names(meff) <- c('Effect')
# scatter plots of marker effects 
plot(abs(meff$Effect), type='o', col = densCols(meff), pch = 20, cex = 1.6, 
     ylab = 'Marker effect',  main = 'Height 30 years')
save.image("bayes_synbreed.RData")
###############################################################################
###############################################################################

#############################################
###  Godness of fit and related statistics
h30_BL$fit
h30_BL$varE # compare to var(y)


### plot of residuals convergence/ TRACE PLOTS FOR LASSO
#pdf("~/Google Drive/Book/Images/ch12_gs/Fig12-4_TracePlots.pdf", width = 7, height = 3)

varE <-scan(paste(outpath,'varE.dat',sep=''))

par(mfrow=c(1,2))
#par(mar=c(4,4.5,0,0)+0.2)#sets margins of plotting area

plot(varE, type="o", col = "gray60", pch = 20, cex = 1,
     ylab = expression(paste(sigma[epsilon]^2)))

abline(h=h30_BL$varE, col=1, lwd=2, lty=1) 
abline(v=h30_BL$burnIn/h30_BL$thin,col=1, lwd=2, lty=1)


### lambda (regularization parameter of BL) 
lambda <- scan(paste(outpath,'ETA_MRK_lambda.dat',sep='')) 

plot(lambda, type='o', col="gray60",cex=.5, 
     ylab=expression(lambda))

abline(v=h30_BL$burnIn/h30_BL$thin, col=1, lwd=2, lty=1) # vertical line
abline(h=h30_BL$ETA$MRK$lambda, col=1, lwd=2, lty=1) #horizontal line

#dev.off()
save.image("bayes_synbreed.RData")



