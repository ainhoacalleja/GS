
#***********************************************************************************************************
#~~~~~~~~~~~~~~~~RKHS~~~~~~~~~
#***********************************************************************************************************

#remove everything in the working environment
rm(list=ls())

# install BGLR package 
require("BGLR")
library(gdata)
library(coda) # provide a set of functions to help the user to monitor the status of Markov Chain convergence (Plummer et al)
library(MCMCpack)

# 1370 individuals
mdat<-read.table("Adjusted_EBV_mdat_F.csv",head=T,sep=",", na.strings=c("",".", "NA"))

# 1382 individuals after we delete 2 samples without phenotypic value
mdat<-mdat[!is.na(mdat$adj_Hjd_17),]
#mdat<-mdat[!is.na(mdat$adj_Pilodyn),]
#mdat<-mdat[!is.na(mdat$adj_Velocity),]
#mdat<-mdat[!is.na(mdat$adj_MOE),]

# noSNP<-c(10, 25, 50, 100, 250, 500, 1000, 2000, 4000, 6000, 8000, 10000, 50000, 100000)
NosubsetSNP=10

# here I need to use drop level to reduce factor.
# set up as factor, it is important in ASReml, ASREml_R and R for estimating factor 
mdat<-drop.levels(mdat)
ped_progeny <- mdat[, c(1,2,3)]


######################################################
# Change folder according your computer!
outpath<-getwd()
load("Norway spruce codeGeno data.rda")
new_ped<-ped_progeny
IND<-match(row.names(gp.num$geno), new_ped$Genotype_id)

#gp.num$geno[1:10, 1:10] ; dim(gp.num$geno)
gp.num$genotype<-gp.num$geno[!is.na(IND), ]; gp.num$genotype<-drop.levels(gp.num$genotype); dim(gp.num$genotype)  # gp.num$genotype 

new_ped$Genotype_id<-as.character(new_ped$Genotype_id) # Genotype_id should be character

gp.num$genotype<-gp.num$genotype[new_ped$Genotype_id,]  # right

#******************************************************************
gp.num$geno<-gp.num$genotype # right

# to see if pedigree has the same order as genotype
identical(as.character(new_ped[,1]), row.names(gp.num$geno))


############################################
### *  GENOTYPES 
genotype <- as.matrix(gp.num$geno)
genotype_b <- as.matrix(gp.num$geno)
# trait for height
EBV<-mdat[, "EBVHeight"]           #^^^^^^^^^3
tt<-as.data.frame(mdat[, "adj_Hjd_17"]) #^^^^^^^^^4
colnames(tt)<-c("adj_Hjd_17")           #^^^^^^^^^5

trait<-names(tt)
y <- tt[, trait]  

samplesize <- ceiling(dim(mdat)[1]/10)
nfolds <- 10 # no of folds

nIter = 150000;  burnIn = 50000
thin=1000

#nIter = 50;  burnIn =10 
#thin=10

results<-as.matrix(1:6,6, 1) 
model <- "RKHS"
set.seed(12345)

#**************************************************
 # start first Loop 1
  # 1st Loop --> FOLDS          
  for(fold in 1:nfolds){
    # Sample from y starting from observation 1, using samplesize
    whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
    yNa <- y              # make a copy of y. yNa and y are numeric  
    yNa[whichNa] <- NA    # set the values of sampled trees to NA
   
    ###########################################################################
    
    #estimate the sample covariance matrix between individuals########
    # ^^^^^here we change it 
    X <- genotype[-whichNa,]  #**** 
    p <- dim(X)[2] #total number of SNPs
    Xs <- scale(X,center=TRUE,scale=TRUE)
    
    A <- tcrossprod(Xs)/p
  
    #load source code for Single SNP analysis
    source('Singlemapping_withoutimpuation.R')
    
    # adjusted phenotypic MOE 
    
    Y <- y[-whichNa] 
    
    #Alternative to LASSO and Stability selection, we can also use the single SNP approach
    
    #Basically, this means that we use linear regression to test each SNP at a time. For each SNP we get a orignal p-value, plus an adjusted p-value by using permutation multiple testing method.
    
    #Permutation test is very time consuming for large data sets. So I make it running in parallel, you need to install package 'doParallel'. Chosse number of cores, based on number of processors of your computer.
    
    #This will take quite long computational time, so please be patient. 
    
    #As discussed, we include the first two PCs as the covariates in the model
    
    # result2 <- singlemapping(X=X,y=Y ,covariate=PC[,c(1,2, 3, 4, 5)])
    result2 <- singlemapping(X=X,y=Y)
    
    
    length(result2$Coef)
    #For adjusted p-values, we may use the classical threshold 0.05. And we will detect the significant SNPs
    
    new<-cbind(as.data.frame(result2$Coef), 1)
    row.names(new)<-colnames(X)
    new<-new[order(new$`result2$Coef`, decreasing = TRUE),]
    
    #*** here we change ^^^^^^^^
    new.X<-genotype_b[,row.names(new)] 
    
    ### *  GENOTYPES 
    genotype <- as.matrix(new.X)
    
    #***************************************************
    ### * BGLR requires matrix of markers. 
    # marker matrix is centered on zero for ease of calculation
    subgenotype<-genotype[, 1:NosubsetSNP]
    geno_centered <- subgenotype- 1
    
    
    X=geno_centered
    p<-ncol(X)
    D<-(as.matrix(dist(X, method="euclidean"))^2)/p
    h<-0.25
    K<-exp(-h*D)
    ETA<-list(K1=list(K=K, model="RKHS"))
    
    
    #_______________________________________________________________
    # end the cut
    ############################################################  
    
    # RKHS   
    # computing the distance matrix and then the kernel
    fmR <- BGLR(y = yNa, response_type = "gaussian", a = NULL, b = NULL, ETA = ETA,
                nIter = nIter, burnIn = burnIn, thin = thin, saveAt = "RKHS_h=0.25_")
    
    # Data frame, combine predicted (yHat and observed y)    
    df1 <- data.frame(fmR$yHat[whichNa], EBV[whichNa], y[whichNa])
    
    colnames(df1) <- c( 'GEBV', 'Phenotype', "adjustedPheno")
    
    # Calculate statistics
    PredAbi = round(cor(df1$adjustedPheno, df1$GEBV), 2) # PredAbi
    Accuracy=round(cor(df1$Phenotype,df1$GEBV), 2)
    RankCor = round(cor(df1$Phenotype,df1$GEBV,method='spearman'), 2)
    Best10 <- round(mean(tail(sort(df1$Phenotype), n = ceiling(nrow(df1)*0.1))), 2) # Best 10%
    Bias   <- round(coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2], 2)
    MSE  <-  round(mean((df1$Phenotype - df1$GEBV)^2), 2) # MSE   # MSE
    
    pred.sum <- rbind(PredAbi, Accuracy, RankCor,Best10, Bias, MSE)
    pred.sum <-round(pred.sum,2) 
    results<-cbind(results, pred.sum) 
  }
  

results  
as.data.frame(apply(results[, -1], 1, mean))  # average of 10 times
sum_results<-round(cbind(as.data.frame(apply(results[, -1], 1, function(x){ mean(x,na.rm=TRUE) })),as.data.frame(apply(results[, -1], 1, function(x){ sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x))) })) ),2)
colnames(sum_results)<-c("Estimate", "SE")
sum_results

write.csv(sum_results, file=paste(trait,model, NosubsetSNP, "largest.csv" , sep=""), row.names=FALSE)
