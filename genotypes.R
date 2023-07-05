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
packages <- c("plyr", "dplyr",  "ggplot2", "RColorBrewer", "LDheatmap", "Hmisc", "GGally", "synbreed")
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

#Setting the working directory
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/Imputation/after_imputation_tassel")

#reading the SNPs data created after imputation with TASSEL. If The file is on .txt format
data <- read.table(file='knn.txt',header=TRUE,row.names=1,stringsAsFactors=F)
save.image(".RData")



#The data in TASSEL has only one nucleotide but a list of the correspondant base pair calls is provided by the software.
#Transforming from one nucleotide to two with the information provided by TASSEL.
data[data=="N"]<-NA
data[data=="A"]<-"AA"
data[data=="C"]<-"CC"
data[data=="G"]<-"GG"
data[data=="T"]<-"TT"
data[data=="R"]<-"AG"
data[data=="Y"]<-"CT"
data[data=="S"]<-"CG"
data[data=="W"]<-"AT"
data[data=="K"]<-"GT"
data[data=="M"]<-"AC"
save.image(".RData")

#Writing the marker data in a .csv file that contains the new information in a new directory
write.csv(data, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/data2.csv",col.names = TRUE)
save.image(".RData")

#Reading the new marker data file
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")
data2 <-read.csv('data2.csv',header=T,row.names=1,  sep=",")
#nacheck <- is.na(data2)
head(data2[, 1:6])  #to check the first 6 columns of the datafile
tail(data2[, 1:6]) #to check the last 6 columns of the data file
save.image(".RData")

#ord.data3<-cbind(rownames(data3)[order(rownames(data3))], data3[order(rownames(data3)),])
#head(ord.data3[, 1:6])
#tail(ord.data3[, 1:6])
data3 <-read.csv('data2.csv',header=T,row.names=1,  sep=",")  #read again the file to keep original untouched data2
head(data3[, 1:6])
class(data3)
str(data3)
save.image(".RData")


data4 <- cbind(Taxa = rownames(data3), data3)  #assign title to rownames column = Taxa
rownames(data4)<-NULL #delete the original rownames column that did not contained a title
head(data4[, 1:6])
save.image(".RData")

#read new datafile with only new ID, which contains the same Taxa column as the previous file 
dataid <-read.csv('indiv.csv',header=T,row.names=,  sep=",")
head(dataid)

datanew <- merge(dataid,data4,by="Taxa") #merging the two dataframes (new ID and old ID) by common column "Taxa" to add the new ID column
head(datanew[, 1:6])  #to check the first 6 columns of the datafile
tail(datanew[, 1:6]) #to check the last 6 columns of the data file

save.image(".RData")

#Now delete the "Taxa" column to keep only the ID column
datanew$Taxa <- NULL  #that will assign no column Taxa to the dataframe
head(datanew[, 1:6])  #to check the first 6 columns of the datafile
tail(datanew[, 1:6]) #to check the last 6 columns of the data file

#Saving the datafile with the new ID column
write.csv(datanew, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/indi_genot.csv",col.names = TRUE)

#Convert the ID column into the rownames column
#row.names(datanew) <- datanew$ID
#datanew[1] <- NULL

#Saving the datafile 
write.csv(datanew, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/indi_genot1.csv",col.names = TRUE)

#Adding pedigree information for individuals with marker information
dataped <-read.csv('dataplan_ped.csv',header=T,row.names=,  sep=",")
head(dataped)

dataped1 <- merge(dataid,dataped,by="ID") #merging the two dataframes by common column "ID" 
dataped1$Taxa <- NULL 
head(dataped1)  #to check the first 6 columns of the datafile
tail(dataped1) #to check the last 6 columns of the data file

#Saving the datafile with the pedigree as "indi_ped1.csv"
write.csv(dataped1, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/indi_ped1.csv",col.names = TRUE)


#data4 <-read.csv('data2.csv',header=T,  sep=",")
#head(data4[, 1:6])  #to check the first 6 columns of the datafile
#tail(data2[, 1:6]) 


