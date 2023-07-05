# Remove everything in the working environment.
rm(list=ls())


#################################################
### Using a function to check if packages are installed. 
# If they are not, the function download them from CRAN and 
# load them into the R session.
## Source: https://gist.github.com/stevenworthington/3178163
# function


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("plyr", "dplyr",  "ggplot2", "RColorBrewer", "LDheatmap", "Hmisc", "GGally", "synbreed", "rrBLUP")
ipak(packages)


#However, we still need to load the libraries
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

#Setting the working directory
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")
save.image(".RData")
load(".RData")


### Load GENOTYPES and PEDIGREE file
# The first 4 columns are ID, Parent1, Parent2, Generation
genoped <- read.csv("genoped.csv", header= T,stringsAsFactors = F, 
                    colClasses = c(Par1 = "character", Par2 = "character"), row.names =1,check.names=F)
# Show first six entries and six columns of genoped
head(genoped[, 1:6])

#rm(pheno)
### Load PEHNOTYPE data
pheno <- read.csv("pheno1.csv", header=TRUE,
                  strip.white = TRUE, na.strings = c("NA",""),
                  stringsAsFactors = F, row.names = 1, check.names=F)
# print the first 6 lines
head(pheno)
tail(pheno)

#check that rownames for geno and pheno are the same
identical(sort(rownames(genoped)),sort(rownames(pheno)))
save.image(".RData")

#create a pedigree file
ped<-create.pedigree(rownames(genoped),genoped$Par1,genoped$Par2,genoped$gener)
summary(ped)
head(ped)
tail(ped)
save.image(".RData")

write.csv(ped, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/ped.csv", row.names=FALSE)

#Create matrix of markers and match with ped$ID
gmat <- as.matrix(subset(genoped, select = -c(Par1, Par2, gener)))
# Check if row names identical?
identical(sort(rownames(gmat)), sort(ped$ID))
# Show first six entries and 5 columns
head(gmat[,1:5])
tail(gmat[,1:5])
save.image(".RData")

gp<-create.gpData(geno=gmat, pheno=pheno, ped=ped, map=NULL)
str(gp)
summary(gp)
save.image(".RData")

#code genotypic data and imputatio when nmissing=10%
gp.num<-codeGeno(gp, label.heter = "alleleCoding", maf=0.01, nmiss = 0.1, impute=TRUE,
                    impute.type = "random", verbose = TRUE)

summary(gp.num)
save(gp.num, file="pine_codeGeno.rda")  #saving the imputation
save.image(".RData")

###############  Relatedness matrices  ##################
#### Pedigree based matrix
Additive<-kin(gp.num, ret = "add") 
summary(Additive)

# Realized (based on markers) relationship matrix.

Realized<-kin(gp.num, ret="realized") #based on the GOF method
summary(Realized)

Realized2 <- kin(gp.num, ret="realizedAB") #based on the GD method
summary(Realized2)

### ... Pedigree based dominance relationship matrix (D).... ###
Dominance <- kin(gp.num, ret="dom")
summary(Dominance)

save.image(".RData")

### Notice that the order of rows/columns are not same for realized and pedigree based matrices
identical(rownames(Realized), rownames(Additive))
identical(sort(rownames(Realized)), sort(rownames(Additive)))

#sort both matrices to match order of pheno object
#this is required later in the cross-validation step
Additive = Additive[rownames(gp.num$pheno), rownames(gp.num$pheno)]
Realized = Realized[rownames(gp.num$pheno), rownames(gp.num$pheno)]

save.image(".RData")


#Moreover, perform spectral decomposition of these matrices.
#We wish to know the eigen-values in order to determine whether the 
#matrices are positive definite and therefore invertable or if there are
#singularities

evalAdd <- eigen(Additive)   #Spectral decomposition
evalDom <- eigen(Dominance)
evalReal <- eigen(Realized)

#inspect the summaries of the eigenvalues
summary(evalAdd$values)
summary(evalDom$values)
summary(evalReal$values)

hist(evalAdd$values)
hist(evalReal$values)
hist(evalG$values)

#~~~~~Code to remove singularities from a matrix ~~~~~~~
#use nearPD function from Matrix package to get a positive definite matrix
RealizedPD = nearPD(Realized, keepDiag = T)
#The output from nearPD is a very specialized object and needs to be 
#reformatted to original matrix format again.
G = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
G = G + diag(0.01, nrow(G)) #A small value is added to the diagonal for safety's sake
attr(G, "dimnames") = RealizedPD[[1]]@Dimnames
class(G) = "relationshipMatrix"
str(G)
summary(G)
save.image(".RData")



#Now perform spectral decomposition on the resultant G-matrix and check
#the summary of its eigenvalues. 
evalG <- eigen(G)
summary(evalG$values)
head(evalG)
save.image(".RData")

#Now also make and inspect heatmaps of your relationship matrices
#~~~~~~~~  plot A-relationship matrix with equal colorkeys ~~~~~~~
pdf(paste("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/AmatrixHeat.pdf",sep=""), width=7, height=7)
par(mar=c(4,4,4,1)+0.1)#sets margins of plotting area
plot(Additive, levelbreaks=seq(-0.2,1.2,length=15))  
title(main='A matrix')
dev.off()  
save.image(".RData")

#~~~~~~~~  plot A-relationship matrix with equal colorkeys ~~~~~~~
pdf(paste("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/GmatrixHeat.pdf",sep=""), width=7, height=7)
par(mar=c(4,4,4,1)+0.1)#sets margins of plotting area
plot(G, levelbreaks=seq(-0.2,1.2,length=15))  
title(main='G matrix')
dev.off()   
save.image(".RData")

#Lets now also compare G coefficients to A coefficients in another way
#by boxplotting
png(paste("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/GvsAboxplot.png",sep=""), width=5, height=5, units ="in", res = 200)
#pdf(paste(gpath,"/Fig11-5_GvsAboxplot.pdf",sep=""), width=5, height=5)
boxplot(as.vector(G) ~ as.vector(Additive), col="gray", xlab = "A matrix", ylab = "G matrix")
dev.off()
save.image(".RData")

#Also compile a summary about G coefficients by each category of A coefficient
by(as.vector(G), as.factor(as.vector(Additive)), FUN = summary)
save.image(".RData")
load(".RData")


#convert from 'dense' matrix representation to 'table' format
Atab = write.relationshipMatrix(Additive, sorting = "ASReml", type = "none")
colnames(Atab)[3] = "A"
head(Atab)

Gtab = write.relationshipMatrix(G, sorting = "ASReml", type = "none")
colnames(Gtab)[3] = "G"
head(Gtab)

#combine A and G coefficients and compare them
AG = merge(Atab, Gtab, all = T)
AG[is.na(AG)] = 0

#show the most extreme differences between A and G
AG$diff = AG$A - AG$G
AG = AG[order(AG$diff),]

round(head(AG),3)

write.csv(AG, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/AG.csv", row.names = FALSE)
save.image(".RData")

#show the pedigrees of some pairs with extreme differences
gp.num$pedigree[gp.num$pedigree$ID %in% rownames(Additive)[c(696, 26)],]

round(tail(AG),3)
gp.num$pedigree[gp.num$pedigree$ID %in% rownames(Additive)[c(606, 607)],]
save.image(".RData")

#Histograms of the kinship coefficients at each matrix
hist(G, main="Realized",breaks = 20,
     xlab="kinship coefficient",ylab="Freq", col="gray", axes=F)
axis(2)
axis(1, at=seq(0,0.75, by=0.25))

hist(Additive,  main="Additive",breaks = 20,
     xlab="kinship coefficient",ylab="Freq", col="gray", axes=F)
axis(2)
axis(1, at=seq(0,0.75, by=0.25))

save.image(".RData")

#------------------------------------------------

#Final preparations for run with ASReml
#rm(pheno)
### Load PEHNOTYPE data
new.ph <- read.csv("pheno2.csv", header=TRUE,
                  strip.white = TRUE, na.strings = c("NA",""),
                  stringsAsFactors = F, row.names = 1, check.names=F)
# print the first 6 lines
head(new.ph)
tail(new.ph)
save.image(".RData")

#check that rownames for geno and pheno are the same
identical(sort(rownames(genoped)),sort(rownames(new.ph)))
save.image(".RData")

#Add row names to the dataframe proper as factors
new.ph$ID = factor(rownames(new.ph))
head(new.ph)  #check the result
tail(new.ph)
save.image(".RData")



#Also create inverses of the A- and G-matrices and reformat them
#to ASReml table format (write.relationshipMatrix)
A.giv <- write.relationshipMatrix(Additive,file=NULL,sorting=c("ASReml"),type=c("ginv"),digits=10)
G.giv <- write.relationshipMatrix(G,file=NULL,sorting=c("ASReml"),type=c("ginv"),digits=10)
names(G.giv) <- c("row", "column", "coefficient")  
head(G.giv)
save.image(".RData")

names(A.giv)<- c("row", "column", "coefficient")  
head(A.giv)
save.image(".RData")

#Now run ASReml-R, both for the A-matrix and the G-matrix. First univariate.
#Factors need to be assign in the phenotypic data to block, plot, Family, Mother (Par1) and father(Par2). 








