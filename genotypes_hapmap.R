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
library(rrBLUP)

#Setting the working directory
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")

################################### MODIFY MARKER DATA FILE FROM TASSEL###########################
library(readr)
knn_hmp <- read_delim("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/knn.hmp.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#reading the SNPs data created after imputation with TASSEL. If The file is on .txt format
#datahmp <- read.table(file='knn.hmp.txt',header=TRUE,row.names=1,stringsAsFactors=F, sep = " ")
View(knn_hmp)
save.image(".RData")

#Saving the datafile in .csv format
write.csv(knn_hmp, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/datahmp.csv",col.names = TRUE)

#reading the .csv datafile
datahmp <-read.csv('datahmp.csv',header=T,row.names=1,  sep=",")
head(datahmp[, 1:6])

#Extracting snp, chromosome and position from the datafile under the new dataframe called map
map<- select(datahmp, rs., chrom, pos)
head(map[, 1:3])
colnames(map)[colnames(map)=="rs."] <- "locus_synb"  #change name of the first column
head(map)

write_csv(map,"~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/map.csv") 

map2 <- map[,-1]   #deleting SNP column from datafile and creating new datafile
rownames(map2) <- map[,1] # creating rownames with the SNP column from before

head(map2)

#saving the map datafile created in .csv format
write.csv(map2, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/map2.csv", col.names = TRUE)


#create a datafile with snps and individuals, using the transpose of the previous datahmp file
gendata1<- t(datahmp)
head(gendata1[, 1:6])

#Droping the columns that I do not need
gendata2 <- gendata1[-c(1,3,4,5,6,7,8,9,10,11,12), ,drop=FALSE]
head(gendata2[, 1:6])

#Converting first row as the column names
colnames(gendata2) = gendata2[1, ] # the first row will be the header
gendata2 = gendata2[-1, ]          # removing the first row.

head(gendata2[, 1:6]) 


#saving the gendata2 file created in .csv format
write.csv(gendata2, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/gendata2.csv", col.names = TRUE)

#reading the datafile
gendata3 <-read.csv('gendata2.csv',header=T,row.names= NULL,  sep=",", stringsAsFactors = FALSE)
head(gendata3[, 1:6])
save.image(".RData")

#The data in TASSEL has only one nucleotide but a list of the correspondant base pair calls is provided by the software.
#Transforming from one nucleotide to two with the information provided by TASSEL.
gendata3[gendata3=="N"]<-NA
gendata3[gendata3=="A"]<-"AA"
gendata3[gendata3=="C"]<-"CC"
gendata3[gendata3=="G"]<-"GG"
gendata3[gendata3=="T"]<-"TT"
gendata3[gendata3=="R"]<-"AG"
gendata3[gendata3=="Y"]<-"CT"
gendata3[gendata3=="S"]<-"CG"
gendata3[gendata3=="W"]<-"AT"
gendata3[gendata3=="K"]<-"GT"
gendata3[gendata3=="M"]<-"AC"
head(gendata3[, 1:6])
save.image(".RData")

#In this dataframe the rownames are assigned as individuals without tittle for the column so I create a new column names
gendata3 <- cbind(names = rownames(gendata3), gendata3)

#Now I can assign a tittle for the column name Taxa.
colnames(gendata3)[2] <- "Taxa"
head(gendata3[, 1:6])

#Now I need to delete the comun names created before
gendata3$names <- NULL
head(gendata3[, 1:6])


#reading the datafile that contains  ID column and Taxa column 
dataid <-read.csv('indiv.csv',header=T,row.names=,  sep=",")
head(dataid)

#Merge dataid and gendata3 by Taxa 
gendata4 <- merge(dataid,gendata3,by="Taxa") #merging the two dataframes (new ID and old ID) by common column "Taxa" to add the new ID column
head(gendata4[, 1:6])  #to check the first 6 columns of the datafile
tail(gendata4[, 1:6]) #to check the last 6 columns of the data file

#Deleting the Taxa column because and keep only the ID column
gendata4$Taxa <- NULL  #that will assign no column Taxa to the dataframe
head(gendata4[, 1:6])  #to check the first 6 columns of the datafile
tail(gendata4[, 1:6]) #to check the last 6 columns of the data file
save.image(".RData")

#Saving the definitive marker file to continue with more anlaysis
write.csv(gendata4, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/gendata.csv",col.names = TRUE)
save.image(".RData")

####merging genotypes and pedigree
#First. Reading the data files gendata.csv and indi_ped2.csv
gendata5 <-read.csv('gendata.csv',header=T,row.names=1,sep=",", stringsAsFactors = FALSE)
head(gendata5[, 1:6])



#reading the pedigree file
#library(readr)
indi_ped2 <- read.csv('indi_ped2.csv',header=T,sep=",",stringsAsFactors = FALSE)
head(indi_ped2[, 1:4])

#Check that the row names are identical for pedigree and genotypes
identical(sort(rownames(gendata5)),sort(rownames(indi_ped2)))

#Merge both files by id 
genoped <- merge(indi_ped2,gendata5,by="ID") #merging the two dataframes by common column "ID" 
head(genoped[, 1:6])  #to check the first 6 columns of the datafile
tail(genoped[, 1:6]) #to check the last 6 columns of the data file

genoped2 <- genoped[,-1]
head(genoped2[, 1:6])

rownames(genoped2) <- genoped[,1]
head(genoped2[ ,1:6])
tail(genoped2[, 1:6])
save.image(".RData")

#Saving the definitive marker+pedigree file to continue with more anlaysis
#write.csv(genoped, "~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/genoped.csv",col.names = TRUE, row.names = TRUE)
write.csv(genoped2, "genoped.csv")
save.image(".RData")
load(".RData")


