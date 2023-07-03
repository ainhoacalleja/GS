###################################################
library(synbreed)
library(asreml)
library(asremlPlus)
library(Matrix)
library(gdata) # for drop.levles
library("matrixcal")

rm(list=ls())   

load("C:/Users/aica/OneDrive - Skogforsk/Documents/DATA/thesis/Project_3_Grundtjärn/analysis/non_additive/NA_analysis_2020/Chen_scprits/IMputedGdata.RData")

# ####################************************************
#... genotype based...
### realized G matrix************************************

ped<-read.csv("C:/Users/aica/OneDrive - Skogforsk/Documents/DATA/thesis/Project_3_Grundtjärn/analysis/non_additive/NA_analysis_2020/Chen_scprits/CV_non_additive/pineped.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)

# we need to delete those bad samples wich are not included in the new G matrix 
delete_ped<-drop.levels(ped[is.na(match(ped$ID, row.names(G))),])

# to get parents
# 46 parents, 731-46= 685 individuals 
ped_parent<-ped[1:46,]
ped_parent<-droplevels(ped_parent)

# to get progeny
# mdat<-read.csv("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP/pine_S1.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)
mdat<-read.csv("C:/Users/aica/OneDrive - Skogforsk/Documents/DATA/thesis/Project_3_Grundtjärn/analysis/non_additive/NA_analysis_2020/Chen_scprits/CV_non_additive/pine.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)
head(mdat)

library(dplyr)
phenodat <- select(mdat,ID,dbh30, dbh36,ht10,ht30,ad_den,MFA,MOEs,MOEd)
head(phenodat)
str(mdat)
mdat[,'ID'] <- as.factor(mdat[,'ID'])

peddata<- select(mdat,ID,Par1,Par2)
ped2<-create.pedigree(ID, Par1, Par2, gener=NULL,sex=NULL,add.ancestors=FALSE,unknown=0)
gp<-create.gpData(geno=NULL, pheno=phenodat, ped=ped_PP, map=NULL)
str(gp)
summary(gp)


sub_mdat<-subset(mdat, !is.element(mdat$ID, delete_ped$ID))

ped_progeny<-sub_mdat[, 1:3]
ped_progeny$ID<-as.character(ped_progeny$ID)

ped_PP<-rbind(ped_parent, ped_progeny[,1:3]) #  new pedigree
names(ped_PP)<-c("ID", "Par1", "Par2")
head(ped_PP)
dim(ped_PP)
identical(sort(as.character(ped_PP$ID)), sort(as.character(row.names(G))))



Realized<-G[as.character(ped_PP$Genotype_id), as.character(ped_PP$Genotype_id)]
summary(Realized)
Realized[1:10, 1:10]
dim(Realized)
Gorder<-row.names(Realized)

# use it for GBLUP in ASReml Standard alone
write.table(Gorder, file = "C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/gorder.txt", sep="\t",col.names = FALSE, row.names = FALSE)

#**************************************** 
# we need all of below scripts in order to run correct in ASReml R 
RealizedPD = nearPD(Realized, keepDiag = T)
G = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
G = G + diag(0.01, nrow(G)) 
attr(G, "dimnames") = RealizedPD[[1]]@Dimnames
class(G) = "relationshipMatrix"


#**************************************** 
# for  ASReml-R
G.giv<-write.relationshipMatrix(G, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)

# For ASReml Stand alone
write.relationshipMatrix(G, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/G-Pedgree.grm", sorting="ASReml",  type="none")



# for extract inbreeding coefficent +1   '
# diag_G<-diag(G)
saveRDS(G, "C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/G_Height_matrix.rds")    

#******** end of additve #********************



#*********************************************
#******************************************
# wrote by chen
# to produce Genomic Dominance matrix 
#******************************************
# to produce genomic Dominance (Dg)
# gp.num$geno[1:5, 1:5]
# X=gp.num$geno
 
 X<-Gdata[as.character(ped_PP$Genotype_id),]

# function for producing dominance matrix, it takes 7 minutes
D <- function(M) {
  M2<-M # the genotype is 0, 1, 2
  
  #formulate the w matrix
  wj<-function(g) { 
    p<-sum(g)/(length(g)*2) 
    q<-1-p
    wij <- ifelse(g == 0,-2*p^2 , ifelse(g == 1, 2*p*q, -2*q^2))
    return(wij)
  }
  w<-apply(M2,2,wj) # w matrix 
  p<-colSums(M2)/(nrow(M2)*2)
  D<-(w%*%t(w))/sum((2*p*(1-p))^2) # sum((2*p*(1-p))^2) it works for two vectors 
  D<-round(D,3)
  return(D)
}

GDom<-D(X)
GDom[1:10, 1:10]

identical(as.character(ped_PP$Genotype_id), row.names(GDom))

RealizedPD = nearPD(GDom, keepDiag = T)
Gd = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
Gd = Gd + diag(0.001, nrow(Gd))  
attr(Gd, "dimnames") = RealizedPD[[1]]@Dimnames
class(Gd) = "relationshipMatrix"

Gd.giv<-write.relationshipMatrix(Gd, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)

# write.relationshipMatrix(Gd, file="C:/Users/zhen0001/Desktop/AINHOA/Gd.giv", sorting="ASReml", type=c("ginv"), digits=10)

# For ASReml Stand alone
write.relationshipMatrix(Gd, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/Gd-Pedgree.grm", sorting="ASReml", type="none")


#**************************************************************
# Epistasis
#**********************************
# install.packages("matrixcalc")
# library("matrixcalc")

Geaa<- hadamard.prod(Realized, Realized )

RealizedPD = nearPD(Geaa, keepDiag = T)
Geaa = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
Geaa = Geaa + diag(0.01, nrow(Geaa))  
attr(Geaa, "dimnames") = RealizedPD[[1]]@Dimnames
class(Geaa) = "relationshipMatrix"

Gead<- hadamard.prod(Realized , GDom )
RealizedPD = nearPD(Gead, keepDiag = T)
Gead = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
Gead = Gead + diag(0.01, nrow(Gead))  
attr(Gead, "dimnames") = RealizedPD[[1]]@Dimnames
class(Gead) = "relationshipMatrix"

Gedd<- hadamard.prod( GDom, GDom )
RealizedPD = nearPD(Gedd, keepDiag = T)
Gedd = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
Gedd = Gedd + diag(0.01, nrow(Gedd))  
attr(Gedd, "dimnames") = RealizedPD[[1]]@Dimnames
class(Gedd) = "relationshipMatrix"

#For ASReml-R
Geaa.giv<- write.relationshipMatrix(Geaa, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)
Gead.giv<-write.relationshipMatrix(Gead, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)
Gedd.giv<-write.relationshipMatrix(Gedd, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)


# For ASReml Stand alone
write.relationshipMatrix(Geaa, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F//Geaa-Pedgree.grm", sorting="ASReml", type="none")
write.relationshipMatrix(Gead, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F//Gead-Pedgree.grm", sorting="ASReml", type="none")
write.relationshipMatrix(Gedd, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F//Gedd-Pedgree.grm", sorting="ASReml", type="none")

# to run additive, dominance, epistasis GBLUP models 
# GBLUP -ADE

mdat<-read.csv("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/pine.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)
mdat<-subset(mdat, !is.element(mdat$ID, delete_ped$ID))
write.csv(mdat, file="C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/pine_d9_ind.csv", row.names = FALSE)

str(mdat)
dim(mdat)
dim(G)


mdat$Genotype_id <-as.factor(mdat$ID)
mdat$Genotype_idd<-as.factor(mdat$ID)
mdat$Genotype_EAA<-as.factor(mdat$ID)
mdat$Genotype_EAD<-as.factor(mdat$ID)
mdat$Genotype_EDD<-as.factor(mdat$ID)

# GBLUP-A
GBLUPA.asr <- asreml(fixed = ht10 ~ 1, 
              random = ~ giv(Genotype_id), 
              ginverse=list(Genotype_id=G.giv),
              rcov= ~ units, 
              data=mdat, 
              na.method.X='include',
              workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them


summary(GBLUPA.asr)

# GBLUP-AD
GBLUPAD.asr <- asreml(fixed = ht10 ~ 1 , 
              random = ~ giv(Genotype_id) + giv(Genotype_idd) , 
              ginverse=list(Genotype_id=G.giv,  Genotype_idd=Gd.giv),
              rcov= ~ units, 
              data=mdat, 
              na.method.X='include',
              workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them

summary(GBLUPAD.asr)

# GBLUP-ADE
GBLUPADE.asr <- asreml(fixed = ht10 ~ 1, 
              random = ~ giv(Genotype_id) + giv(Genotype_idd) + giv(Genotype_EAA) + giv(Genotype_EAD)+giv(Genotype_EDD), 
              ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv, Genotype_EAA=Geaa.giv, Genotype_EAD=Gead.giv, Genotype_EDD=Gedd.giv),
              rcov= ~ units, 
              data=mdat, 
              na.method.X='include',
              workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them


summary(GBLUPADE.asr)

# -2894.6001   1458.6794   684  17:53:22     0.5 GBLUP-A
# -2894.6001   1458.6789   684  17:54:29     2.8 GBLUP-AD
# -2894.5355   1359.8738   684  18:07:46    63.6 GBLUP