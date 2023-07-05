setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes")
#data <- read.table(file='manhattan.txt',header=TRUE,row.names=1,stringsAsFactors=F)

#library(vcfR)
#data <-read.vcfR("/Users/ainhoa/Documents/PhD_UPSC/Project_3_Grundtjärn/Imputation/Imputed_K-near_NoHeader_mafo.01_SNPs_only.vcf")
#write.vcf(data, "data.vcf.gz")


library(vcfR)
vcf <- read.vcfR("/Users/ainhoa/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/genotype_changes/Imputed_KNN.vcf")
#write.vcf(data, "data.vcf.gz")1.vcf")

summary(vcf)

library(VariantAnnotation)


library(synbreed)

data[data=="N"]<-NA
data[data=="A"]<-"AA"

nacheck <- is.na(data)
sum(data)


