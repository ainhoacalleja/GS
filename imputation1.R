# Remove everything in the working environment.
rm(list=ls())

#Setting the working directory
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")

load("imputation1.RData")


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


########################################### IMPUTATION WITH SYNBREED #################################################
####merging genotypes and pedigree
#First. Reading the data files gendata.csv and indi_ped2.csv
gendata<-read.csv('gendata.csv',header=T,row.names=1, sep=",",
                          stringsAsFactors = FALSE)
head(gendata[, 1:6])
save.image("imputation1.RData")

#reading the pedigree file
ped1 <- read.csv('indi_ped2.csv',header=T,sep=",", colClasses=c(Par1="character",
                Par2="character"), row.names = 1, stringsAsFactors = FALSE)
head(ped1[, 1:3])
save.image("imputation1.RData")


identical(sort(rownames(gendata)),sort(rownames(ped1)))



#Check that the row names are identical for pedigree and genotypes
save.image("imputation1.RData")
ped <-create.pedigree(rownames(ped1),ped1$Par1, ped1$Par2, ped1$gener)
summary(ped)
head(ped)
tail(ped)
save.image("imputation1.RData")

gmat<- as.matrix(gendata)
head(gmat[, 1:6])
save.image("imputation1.RData")

identical(sort(rownames(gmat)),sort(rownames(ped)))

gp1 <- create.gpData(geno=gmat, pheno=NULL, ped=ped, map=NULL)
str(gp1)
summary(gp1)

#code genotypic data and imputatio when nmissing=50%
gp1.coded<-codeGeno(gp1, label.heter = "alleleCoding", maf=0.01, nmiss = 0.5, impute=TRUE,
                    impute.type = "random", verbose = TRUE)

summary(gp1.coded)
save.image("imputation1.RData")

#code genotypic data and imputatio when nmissing=30%
gp2.coded<-codeGeno(gp1, label.heter = "alleleCoding", maf=0.01, nmiss = 0.3, impute=TRUE,
                    impute.type = "random", verbose = TRUE)

summary(gp2.coded)
save.image("imputation1.RData")

#code genotypic data and imputatio when nmissing=20%
gp3.coded<-codeGeno(gp1, label.heter = "alleleCoding", maf=0.01, nmiss = 0.2, impute=TRUE,
                    impute.type = "random", verbose = TRUE)

summary(gp3.coded)
save.image("imputation1.RData")

#code genotypic data and imputatio when nmissing=10%
gp4.coded<-codeGeno(gp1, label.heter = "alleleCoding", maf=0.01, nmiss = 0.1, impute=TRUE,
                    impute.type = "random", verbose = TRUE)

summary(gp4.coded)

save.image("imputation1.RData")

#saving the objects of imputation for future use.
save(gp1.coded, file="imputation1.rda")
save(gp2.coded, file="imputation2.rda")
save(gp3.coded, file="imputation3.rda")
save(gp4.coded, file="imputation4.rda")

save.image("imputation1.RData")

####################################################################################################################

######## To perform imputation with rrBLUP the genotypes should be code as {-1,0,1}, so I must recode my data
###### reading the txt datafile obtained from TASSEl where genotypes are coded as {0,1,2}

genotdat <- read.delim("reference_probability.txt", header=F, sep="\t", stringsAsFactors = FALSE, row.names = 1)
head(genotdat[, 1:6])  #the column name are not correct
colnames(genotdat) = genotdat[1, ] # the first row will be the header
genotdat = genotdat[-1, ]  #remove the row with incorrect column names
head(genotdat[, 1:6]) 

new_genot <- cbind(Taxa = rownames(genotdat), genotdat)  #assign title to rownames column = Taxa
rownames(new_genot)<-NULL #delete the original rownames column that did not contained a title
head(new_genot[, 1:6])
save.image("imputation1.RData")

#read datafile with only new ID, which contains the same Taxa column as the previous file 
dataid <-read.csv('indiv.csv',header=T,row.names=,  sep=",")
head(dataid)

genodata <- merge(dataid,new_genot,by="Taxa") #merging the two dataframes (new ID and old ID) by common column "Taxa" to add the new ID column
head(genodata[, 1:6])  #to check the first 6 columns of the datafile
tail(genodata[, 1:6]) #to check the last 6 columns of the data file

save.image("imputation1.RData")

#Now delete the "Taxa" column to keep only the ID column
genodata$Taxa <- NULL  #that will assign no column Taxa to the dataframe
head(genodata[, 1:6])  #to check the first 6 columns of the datafile
tail(genodata[, 1:6]) #to check the last 6 columns of the data file

#Saving the datafile  
write.csv(genodata, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/genotypes.csv",col.names = TRUE)
str(genodata)


#To not modify my dataframe genodata I will create a new dataframe
genodata1<-genodata
#Transforming from numeric from TASSEL {0, 0.5, 1} to the numeric format for rrBLUP package {-1, 0, 1}.
genodata1[genodata1=="0"]<- "-1"
genodata1[genodata1=="0.5"]<- "0"
head(genodata1[, 1:6])  #to check the first 6 columns of the datafile
tail(genodata1[, 1:6]) #to check the last 6 columns of the data file
save.image("imputation1.RData")

#Saving the datafile  
write.csv(genodata1, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/genotypes_rrBLUP.csv", row.names = FALSE)
save.image("imputation1.RData")

#To start the session a different date without running all scripts again
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")
save.image("imputation1.RData")

########################### imputation with rrBLUP ################################
genodata2 <- read.csv("genotypes_rrBLUP.csv", header = T, 
                    stringsAsFactors = F, row.names = 1, check.names = F)

head(genodata2[, 1:6])
str(genodata2)
save.image("imputation1.RData")

markerimput<-A.mat(genodata2,min.MAF=0.01,max.missing=0.1,impute.method="EM",tol=0.02,shrink=FALSE,
                    return.imputed=TRUE)

markerimput2<-A.mat(genodata2,min.MAF=0.01,max.missing=0.05,impute.method="EM",tol=0.02,shrink=FALSE,
                   return.imputed=TRUE)

save.image("imputation1.RData")

# markers_imput <- markerimput$imputed
# head(markers_imput[,1:6])

markers_imput2 <-markerimput2$imputed

# Gm<-markerimput$A
# head(Gm[,1:6])

Gm2<-markerimput2$A
head(Gm2[,1:6])

save.image("imputation1.RData")
load("imputation1.RData")

write.csv(Gm,"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/Gmat_rrBLUP.csv", row.names = TRUE )
write.csv(Gm2,"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/Gmat2_rrBLUP.csv", row.names = TRUE )


dim(genodata2)
# dim(markers_imput)
dim(markers_imput2)
#Create a csv file with the imputed SNPs (8719 SNPs), without an ID column, just the numbers of each ID
#write.csv(markers_imput, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/markers_imput.csv", row.names = TRUE)

write.csv(markers_imput2, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/markers_imput2.csv", row.names = TRUE)
save.image("imputation1.RData")



#Remove NA values from the imputed matrix of SNPs (it takes long time)
# myImputedData<-na.omit(markers_imput)
# dim(myImputedData)
# write.csv(myImputedData, "Imputed_file.csv")
# save.image("imputation1.RData")

myImputedData2<-na.omit(markers_imput2)
dim(myImputedData2)
write.csv(myImputedData2, "Imputed_file2.csv")
save.image("imputation1.RData")

# for (i in 1: nrow(markers_imput)) {
#   for (j in 1: ncol(markers_imput)){
#     myImputedData<- markers_imput [which(rowSums(markers_imput[,],
#                                                  na.rm = FALSE, dims=1) !="NA") ,]
#   }
# }

##########################################################################


head(gp4.coded$geno[,1:6])
tail(gp4.coded$geno[,1:6])

#rm(codedgeno)
write.csv(gp4.coded$geno, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/codedgeno.csv", row.names = TRUE)
save.image("imputation1.RData")
###########################################################################


head(myImputedData2[,1:6])
genod2<-as.data.frame(round(myImputedData2,digits=0)) #rounding of imputed values
head(genod2[,1:6])
write.csv(genod2, "genod2.csv")  #imputed genotypes rounded as saved coded as -1, 0, 1
genod3<-genod2+1 #add 1 to get genotype codes equat to 0, 1 and 2
head(genod3[,1:6])
tail(genod3[,1:6])
save.image("imputation1.RData")

write.csv(genod3, "imputed_round2.csv")
save.image("imputation1.RData")

     