library(asreml)
library(lattice)
library(gdata) # for drop.levels and read xls file 
library(lattice)
require(plyr)

rm(list=ls())
marker<-c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10")


func_PreAbi_add <- function(xx)
{
  return(data.frame(COR = cor(xx[,2], xx$pred_add, use="complete.obs")))
}

func_PreAbi_dom<- function(xx)
{
  return(data.frame(COR = cor(xx[,2], xx$solution_add+xx$solution_dom, use="complete.obs")))
}

func_PreAbi_epi<- function(xx)
{
  return(data.frame(COR = cor(xx[,2], xx$pred_epi, use="complete.obs")))
}

load("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/before_fitering_chen/Ainhoa/IMputedGdata.RData")

ped<-read.csv("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AINHOA_GS_Nonadditive/AINHOA-1/pineped.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)

# we need to delete those bad samples wich are not included in the new G matrix 
delete_ped<-drop.levels(ped[is.na(match(ped$ID, row.names(G))),])

# to get parents
# 46 parents, 731-46= 685 individuals 
ped_parent<-ped[1:46,]
ped_parent<-droplevels(ped_parent)

# to get progeny
# mdat<-read.csv("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP/pine_S1.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)
mdat<-read.csv("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/pine.csv",  header=TRUE,sep=",", na.strings = c("", ".", "NA"),  strip.white=TRUE)

mdat<-subset(mdat, !is.element(mdat$ID, delete_ped$ID))
mdat<-droplevels(mdat)
# sort by Site, this is required by ASReml-R when we tried to fit two residuals for each of two sites 
# mdat<-mdat[order(mdat$Site),];dim(mdat)

# in order to accurately estimate each parameters, 
# we only need those individuals which have phenotypic values
# of course, Asreml can estimate any values without phenotypic values
# mdat<-mdat[!is.na(mdat$ht10),]

# here I need to use drop level to reduce factor.
# set up as factor, it is important in ASReml, ASREml_R and R for estimating factor 
# mdat<-drop.levels(mdat) # need package "gdata"
head(mdat); dim(mdat)

# to prepare the inverse G matrix (G.ginv) for asreml-r, here in mac, i can not use synbreed package, when i use R 3.0 version
# also save G.ginv file as .txt could not be uses for ASReml-R and ASReml, bur G.grm works in ASReml-R
# so this step will do in my pc and save as. RDS

G.giv<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/G_Height.rds")
head(G.giv)

# matrix format, in order to esitmate inbreeding coefficient 
G<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/G_Height_matrix.rds")
diag_G<-diag(G)

str(diag_G)

Gd.giv<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/Gd_Height.rds")
head(Gd.giv)

# epistatic matrices
Geaa.giv<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/Geaa.rds")
Gead.giv<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/Gead.rds")
Gedd.giv<-readRDS("C:/Users/zhen0001/Documents/Chen Zhi-qiang/chen/11_4 AINHOA_GS/AInhoa_GS_Nonadditive/GBLUP-F/Gedd.rds")



################################################## 
### ~~~~ Create a 10-Fold Cross-Validation set 
# Create a copy of the phenotype data 
# assign fold to each individual
# we can repeat 10 times;
############################************************************

set.seed(110716)
results<-as.matrix(1:9,9, 1) 
# traits will be used in this trials
# dbh30	dbh36	ht10	ht30	ad_den	MFA	MOEs	MOEd

i=5 # which means dbh30, if run anothr trait, just change i=6 and so .....


#********************************************************************************************

# for (i in 1:1) {

nfold=10  # new of fold


# to get a data frame for phenotype value (ph)
new.ph = mdat[, c(1:4, i) ]
new.ph[, 5]<-as.numeric(new.ph[,5])

# in order to do produce NA for each fold,  
# count how many columns in new.ph data frame 
ncolumn<-dim(new.ph)[2]

# to produce fold column
new.ph$fold = sample(1:nfold,size = dim(new.ph)[1], replace = T)

# here is just a multple for matrices to produce the same number of columns (nfold) as column one
new.y = as.matrix(new.ph[,5, drop = F])%*%t(rep(1, nfold))   
colnames(new.y) = paste0("Y", 1:nfold)
new.ph = cbind(new.ph, new.y); head(new.ph)

# set one of the fold observations to NA for each line corresponding to its fold number
# jump first n columns of datafile + 1 column for fold
jump_Nocolumn=ncolumn+1
 missingInFold = function(r){
  fold = as.numeric(r["fold"])
  
  index = fold + jump_Nocolumn
  r[index] = NA
  return(r)
 }
 
 
new.ph2 = t(apply(new.ph, 1, missingInFold))
new.ph.df = data.frame(new.ph2)
head(new.ph.df); dim(new.ph.df); str(new.ph.df)
new.ph.df[, 5:16]<-apply(new.ph.df[,5:16],2, as.numeric)

new.ph.df$Genotype_id<-as.factor(mdat$ID)
new.ph.df$Genotype_idd<-as.factor(mdat$ID)
new.ph.df$Genotype_idEaa<-as.factor(mdat$ID)
new.ph.df$Genotype_idEad<-as.factor(mdat$ID)
new.ph.df$Genotype_idEdd<-as.factor(mdat$ID)
head(new.ph.df)
str(new.ph.df)

#************************************************************************
# GBLUP-ADE model for scots Pine
#***********************************************************
# execute ASReml for each fold data set and save results in a list

result.list = list()
pred.list = list()


# from the column which is repeated
for (trait in names(new.ph.df)[(jump_Nocolumn+1):(jump_Nocolumn+nfold)]) {
  asr <- asreml(fixed = get(trait) ~ 1, 
                random = ~giv(Genotype_id) + giv(Genotype_idd) + giv(Genotype_idEaa) + giv(Genotype_idEad)+ giv(Genotype_idEdd),
                ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv, Genotype_idEaa=Geaa.giv, Genotype_idEad=Gead.giv, Genotype_idEdd=Gedd.giv),
                rcov= ~ units, 
                data=new.ph.df, 
                na.method.X='include',
                workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them
  result.list[[trait]] = asr
  pred.list[[trait]] = summary(asr, all = T)$coef.rand
  rm(list = "asr")
}

#****************************************
# example, to check first cross-validation
temp.asr <- asreml(fixed = Y1 ~ 1, 
                   random = ~giv(Genotype_id) + giv(Genotype_idd) + giv(Genotype_idEaa) + giv(Genotype_idEad)+ giv(Genotype_idEdd),
                   ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv, Genotype_idEaa=Geaa.giv, Genotype_idEad=Gead.giv, Genotype_idEdd=Gedd.giv),
                   rcov= ~ units, 
                   data=new.ph.df, 
                   na.method.X='include',
                   workspace=20e8  ) # 32G rum, maximum space is 40e8, here i only use half of them
# summary(temp.asr)

# names(vcs) =  c( "Vadd",  "Vd", "Vep", "Vres")

get.varcomps = function(c){
  #c is a component of result.list
  vcs = summary(c)$varcomp$component
  names(vcs) =  c( "Vadd",  "Vd", "Veaa", "Vead", "Vedd", "Vres")
  return(vcs)
}
varcomps = sapply(result.list, FUN = get.varcomps)
varcomps = as.data.frame(t(varcomps))
varcomps$fold = as.numeric(gsub("Y", "", rownames(varcomps)))

#get the mean for each estimation (training) set for each fold to use as mean of predictions
mu = colMeans(new.ph.df[, marker],na.rm = T)
mu = as.data.frame(mu)
mu$fold = as.numeric(gsub("Y", "", rownames(mu)))

vc.mu = merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
#************* we should change it if we want to get at(site,handanberg)
#******************* giv(Genotype_idd)_
new.ph.df$ID_add<-paste("giv(Genotype_id)_", new.ph.df$Genotype_id,sep="")    # main
new.ph.df$ID_Dom<-paste("giv(Genotype_idd)_", new.ph.df$Genotype_id,sep="")   # dominance
new.ph.df$ID_Eaa<-paste("giv(Genotype_idEaa)_", new.ph.df$Genotype_id,sep="") # epistasis additive X additive 
new.ph.df$ID_Ead<-paste("giv(Genotype_idEad)_", new.ph.df$Genotype_id,sep="") # dominance additive x dominance 
new.ph.df$ID_Edd<-paste("giv(Genotype_idEdd)_", new.ph.df$Genotype_id,sep="") # dominance dominance x dominance 

dim(new.ph.df)
tempt<-data.frame()
tempt1<-data.frame()
tempt2<-data.frame()
tempt3<-data.frame()
tempt4<-data.frame()

# for additive  
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_add"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  a<- pred.list[[trait]][ID2,]
  tempt<-rbind(tempt,a) 
}


# for dominance
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_Dom"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt1<-rbind(tempt1,d) 
}  

# for epistasis additive by additive
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_Eaa"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt2<-rbind(tempt2,d) 
}  


# for epistasis additive by dominance
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_Ead"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt3<-rbind(tempt3,d) 
}  

# for epistasis dominance by dominance
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_Edd"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt4<-rbind(tempt4,d) 
}  

tempt<-cbind(tempt, tempt1, tempt2, tempt3, tempt4)
colnames(tempt)<-c("solution_add","std error_add", "z ratio_add",  "solution_dom","std error_dom", "z ratio_dom", "solution_Eaa","std error_Eaa", "z ratio_Eaa", "solution_Ead","std error_Ead", "z ratio_Ead", "solution_Edd","std error_Edd", "z ratio_Edd") 
# n1<-length(diag_G)

preds<-cbind(tempt, new.ph.df[, "Genotype_id"], as.matrix(diag_G[47:n1]))
colnames(preds)<-c("solution_add","std error_add", "z ratio_add",  "solution_dom","std error_dom", "z ratio_dom", "solution_Eaa","std error_Eaa", "z ratio_Eaa", "solution_Ead","std error_Ead", "z ratio_Ead", "solution_Edd","std error_Edd", "z ratio_Edd", "Genotype_id","inbreeding" ) 

# this script is not right. 
# preds.obs = merge(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds), by  = 0); dim(preds.obs) # here is noly 1381 
preds.obs <- cbind(new.ph.df[, 5:6], as.data.frame(preds))

# merge pred.obs with varcomps by fold to get the variance components for each fold
preds.obs <- merge(preds.obs, vc.mu, by = "fold")
head(preds.obs)
 
preds.obs$pred_add = preds.obs$mu + preds.obs$solution_add  # additive model 
preds.obs$pred_dom<-preds.obs$mu+ preds.obs$solution_add+preds.obs$solution_dom # dominiance
preds.obs$pred_epi<-preds.obs$mu+ preds.obs$solution_add+preds.obs$solution_dom+preds.obs$solution_Eaa+preds.obs$solution_Ead+preds.obs$solution_Edd # epistasis



# use round function to round numeric columns 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

preds.obs <- round_df(preds.obs, digits=3) 


pred.obs_backup1<-preds.obs

temp.name<-colnames(pred.obs_backup1)[2]
str(pred.obs_backup1)

#****************** Predictive ability *******************
# with additive model 
PA_A= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_add, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_add)
SE_Accuracy_A<-sd(temp$COR)/sqrt(10)
a1<-round(mean(temp$COR),3);   b1<-round(SE_Accuracy_A,3)

# with dominance model 
PA_AD= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_dom, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_dom)
SE_Accuracy_AD<-sd(temp$COR)/sqrt(10)
a2<-round(mean(temp$COR),3);   b2<-round(SE_Accuracy_AD,3)


# with epistasis model , we only use this in the table as biyue did in her paper
PA_ADE= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_dom, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_epi)
SE_Accuracy_ADE<-sd(temp$COR)/sqrt(10)
a3<-round(mean(temp$COR),3);   b3<-round(SE_Accuracy_ADE,3)



print(c(a1, b1, a2, b2, a3, b3))


#********************************************************************************
# GBLUP-AD model in scots Pine 
# execute ASReml for each fold data set and save results in a list
result.list = list()
pred.list = list()


# from the column which is repeated
for (trait in names(new.ph.df)[(jump_Nocolumn+1):(jump_Nocolumn+nfold)]) {
  asr <- asreml(fixed = get(trait) ~ 1, 
                random = ~giv(Genotype_id) + giv(Genotype_idd) ,
                ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv),
                rcov= ~ units, 
                data=new.ph.df, 
                na.method.X='include',
                workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them
  result.list[[trait]] = asr
  pred.list[[trait]] = summary(asr, all = T)$coef.rand
  rm(list = "asr")
}

# temp.asr <- asreml(fixed = Y1 ~ 1, 
#               random = ~giv(Genotype_id) + giv(Genotype_idd) ,
#               ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv),
#               rcov= ~ units, 
#               data=new.ph.df, 
#               na.method.X='include',
#               workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them
# summary(temp.asr)

#   names(vcs) =  c( "Vadd",  "Vd", "Vep", "Vres")

get.varcomps = function(c){
  #c is a component of result.list
  vcs = summary(c)$varcomp$component
  names(vcs) =  c( "Vadd",  "Vd",  "Vres")
  return(vcs)
}
varcomps = sapply(result.list, FUN = get.varcomps)
varcomps = as.data.frame(t(varcomps))
varcomps$fold = as.numeric(gsub("Y", "", rownames(varcomps)))

#get the mean for each estimation (training) set for each fold to use as mean of predictions
mu = colMeans(new.ph.df[, marker],na.rm = T)
mu = as.data.frame(mu)
mu$fold = as.numeric(gsub("Y", "", rownames(mu)))

vc.mu = merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
#************* we should change it if we want to get at(site,handanberg)
#******************* giv(Genotype_idd)_
new.ph.df$ID_add<-paste("giv(Genotype_id)_", new.ph.df$Genotype_id,sep="") # main
new.ph.df$ID_Dom<-paste("giv(Genotype_idd)_", new.ph.df$Genotype_id,sep="") # dominance



dim(new.ph.df)
tempt<-data.frame()
tempt1<-data.frame()

# for additive  
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_add"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  a<- pred.list[[trait]][ID2,]
  tempt<-rbind(tempt,a) 
}


# for dominance
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_Dom"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt1<-rbind(tempt1,d) 
}  

tempt<-cbind(tempt, tempt1)
colnames(tempt)<-c("solution_add","std error_add", "z ratio_add",  "solution_dom","std error_dom", "z ratio_dom") 
n1<-length(diag_G)

preds<-cbind(tempt, new.ph.df[, "Genotype_id"], as.matrix(diag_G[47:n1]))
colnames(preds)<-c("solution_add","std error_add", "z ratio_add", "solution_dom","std error_dom", "z ratio_dom", "Genotype_id","inbreeding" ) 

# this script is not right. 
# preds.obs = merge(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds), by  = 0); dim(preds.obs) # here is noly 1381 
preds.obs <- cbind(new.ph.df[, 5:6], as.data.frame(preds))

# merge pred.obs with varcomps by fold to get the variance components for each fold
preds.obs <- merge(preds.obs, vc.mu, by = "fold")
head(preds.obs)

preds.obs$pred_add = preds.obs$mu + preds.obs$solution_add
preds.obs$pred_dom<-preds.obs$mu+ preds.obs$solution_add+preds.obs$solution_dom



# use round function to round numeric columns 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

preds.obs <- round_df(preds.obs, digits=3) 


pred.obs_backup1<-preds.obs

temp.name<-colnames(pred.obs_backup1)[2]
str(pred.obs_backup1)
#****************** Predictive ability *******************
# with additive model 
PA_A= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_add, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_add)
SE_Accuracy_A<-sd(temp$COR)/sqrt(10)
a1<-round(mean(temp$COR),3);   b1<-round(SE_Accuracy_A,3)

# with dominance model 
PA_AD= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_dom, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_dom)
SE_Accuracy_AD<-sd(temp$COR)/sqrt(10)
a2<-round(mean(temp$COR),3);   b2<-round(SE_Accuracy_AD,3)

print(c(a1, b1, a2, b2))

#*******************************************_________________________________________________
# ##############################################                  GBLUP+A
# execute ASReml for each fold data set and save results in a list
result.list = list()
pred.list = list()

# from the column which is repeated
for (trait in names(new.ph.df)[(jump_Nocolumn+1):(jump_Nocolumn+nfold)]) {
  asr <- asreml(fixed = get(trait) ~ 1, 
                random = ~giv(Genotype_id) ,
                ginverse=list(Genotype_id=G.giv),
                rcov= ~ units, 
                data=new.ph.df, 
                na.method.X='include',
                workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them
  result.list[[trait]] = asr
  pred.list[[trait]] = summary(asr, all = T)$coef.rand
  rm(list = "asr")
}

# temp.asr <- asreml(fixed = Y1 ~ 1, 
#               random = ~giv(Genotype_id) + giv(Genotype_idd) ,
#               ginverse=list(Genotype_id=G.giv, Genotype_idd=Gd.giv),
#               rcov= ~ units, 
#               data=new.ph.df, 
#               na.method.X='include',
#               workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them
# summary(temp.asr)

#   names(vcs) =  c( "Vadd",  "Vd", "Vep", "Vres")

get.varcomps = function(c){
  #c is a component of result.list
  vcs = summary(c)$varcomp$component
  names(vcs) =  c( "Vadd",   "Vres")
  return(vcs)
}
varcomps = sapply(result.list, FUN = get.varcomps)
varcomps = as.data.frame(t(varcomps))
varcomps$fold = as.numeric(gsub("Y", "", rownames(varcomps)))

#get the mean for each estimation (training) set for each fold to use as mean of predictions
mu = colMeans(new.ph.df[, marker],na.rm = T)
mu = as.data.frame(mu)
mu$fold = as.numeric(gsub("Y", "", rownames(mu)))

vc.mu = merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
#************* we should change it if we want to get at(site,handanberg)
#******************* giv(Genotype_idd)_
new.ph.df$ID_add<-paste("giv(Genotype_id)_", new.ph.df$Genotype_id,sep="") # main


dim(new.ph.df)
tempt<-data.frame()

# for additive  
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID_add"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  a<- pred.list[[trait]][ID2,]
  tempt<-rbind(tempt,a) 
}



tempt<-cbind(tempt, tempt1)
colnames(tempt)<-c("solution_add","std error_add", "z ratio_add") 
n1<-length(diag_G)

preds<-cbind(tempt, new.ph.df[, "Genotype_id"], as.matrix(diag_G[47:n1]))
colnames(preds)<-c("solution_add","std error_add", "z ratio_add", "Genotype_id","inbreeding" ) 

preds.obs <- cbind(new.ph.df[, 5:6], as.data.frame(preds))

# merge pred.obs with varcomps by fold to get the variance components for each fold
preds.obs <- merge(preds.obs, vc.mu, by = "fold")
head(preds.obs)

preds.obs$pred_add = preds.obs$mu + preds.obs$solution_add

# use round function to round numeric columns 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

preds.obs <- round_df(preds.obs, digits=3) 

pred.obs_backup1<-preds.obs

temp.name<-colnames(pred.obs_backup1)[2]
str(pred.obs_backup1)

#****************** Predictive ability *******************
# with additive model 
PA_A= cor(pred.obs_backup1[,2], pred.obs_backup1$pred_add, use="complete.obs")

temp<-ddply(pred.obs_backup1, .(fold), func_PreAbi_add)
SE_Accuracy_A<-sd(temp$COR)/sqrt(10)
a1<-round(mean(temp$COR),3);   b1<-round(SE_Accuracy_A,3)

print(c(a1, b1))



###################***************************************************___________________ ABLUP+AD

# produce a pedigree matrix
ped<-mdat[,c(1,2,3)]

# produce an inverse of a A matrix
# actually it is based on the factors to produce a .ginv file 
pedinv <- asreml.Ainverse(ped)$ginv 
A_inbreeding<-asreml.Ainverse(ped)$inbreeding

# to produce domiance matrix
pedigree<-asreml.Ainverse(ped)$pedigree
pedigree$gener<-c(rep(0, 55),rep(1,dim(pedigree)[1]-55) )
load("/Users/zhen0001/Box Sync/Genomic_selection_Fikret/data/Norway spruce codeGeno data.rda")
gp.num$pedigree<-pedigree
colnames(gp.num$pedigree)=c("ID"  ,  "Par1",  "Par2",  "gener")
Dominance <- kin(gp.num, ret="dom")
dim(Dominance) 
RealizedPD = nearPD(Dominance, keepDiag = T)
Dom = matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
Dom = Dom + diag(0.01, nrow(Dom)) 
attr(Dom, "dimnames") = RealizedPD[[1]]@Dimnames
class(Dom) = "relationshipMatrix"
#**************************************** end
Dom.giv<-write.relationshipMatrix(Dom, file=NULL, sorting="ASReml", type=c("ginv"), digits=10)
result.list = list()
pred.list = list()

# from the column which is repeated
for (trait in names(new.ph.df)[(jump_Nocolumn+1):(jump_Nocolumn+nfold)])
{
  asr <- asreml(fixed = get(trait) ~ 1, 
                random = ~ ped(Genotype_id) + diag(Site):ped(Genotype_id) + diag(Site):giv(Genotype_idd) ,
                ginverse=list(Genotype_id=pedinv,Genotype_idd= Dom.giv),
                rcov= ~ at(Site):units, 
                data=new.ph.df, 
                na.method.X='include',
                workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them + at(Site,1):giv(Genotype_idd) ,Genotype_idd= Gd.giv
  
  result.list[[trait]] = asr
  pred.list[[trait]] = summary(asr, all = T)$coef.rand
  rm(list = "asr")
}

# asr1 <- asreml(fixed = get(trait) ~ 1, 
#               random = ~ ped(Genotype_id) + diag(Site):ped(Genotype_id) + diag(Site):giv(Genotype_idd) ,
#               ginverse=list(Genotype_id=pedinv,Genotype_idd= Dom.giv),
#               rcov= ~ at(Site):units, 
#               data=new.ph.df, 
#               na.method.X='include',
#               workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them + at(Site,1):giv(Genotype_idd) ,Genotype_idd= Gd.giv
# summary(asr1)
# 
########### Summarize cross-validations ###########################

# extract the variance components estimates from each fold
get.varcomps = function(c){
  #c is a component of result.list
  vcs = summary(c)$varcomp$component
  names(vcs) =  c( "Vadd_main",  "Vadd_GXE_SH","Vadd_GXE_SV", "Vd_SH", "Vd_SV", "Vres_SV", "Vres_SH")
  return(vcs)
}
varcomps = sapply(result.list, FUN = get.varcomps)
varcomps = as.data.frame(t(varcomps))
varcomps$fold = as.numeric(gsub("Y", "", rownames(varcomps)))

#get the mean for each estimation (training) set for each fold to use as mean of predictions
mu = colMeans(new.ph.df[, marker],na.rm = T)
mu = as.data.frame(mu)
mu$fold = as.numeric(gsub("Y", "", rownames(mu)))

vc.mu = merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
#************* we should change it if we want to get at(site,handanberg)
#******************* giv(Genotype_idd)_
new.ph.df$ID2_main<-paste("ped(Genotype_id)_", new.ph.df$Genotype_id,sep="") 
new.ph.df$ID2_add_GXE_SV<-paste("Site_Vinden:ped(Genotype_id)_", new.ph.df$Genotype_id,sep="")  # site Vindeln
new.ph.df$ID2_add_GXE_SH<-paste("Site_Hadanberg:ped(Genotype_id)_", new.ph.df$Genotype_id,sep="")  # Hadanberg
new.ph.df$ID2_Dom_SH<-paste("Site_Hadanberg:giv(Genotype_idd)_", new.ph.df$Genotype_id,sep="")  # Hadanberg
new.ph.df$ID2_Dom_SV<-paste("Site_Vinden:giv(Genotype_idd)_", new.ph.df$Genotype_id,sep="")  # site Vindeln


dim(new.ph.df)
tempt<-data.frame()
tempt1<-data.frame()
tempt2<-data.frame()
tempt3<-data.frame()
tempt4<-data.frame()
# for additive  in main
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_main"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  a<- pred.list[[trait]][ID2,]
  tempt<-rbind(tempt,a) 
}
# for GXE_additive  in Hadanberg
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_add_GXE_SH"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  b<- pred.list[[trait]][ID2,]
  tempt1<-rbind(tempt1,b)   
}

for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_add_GXE_SV"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  c<- pred.list[[trait]][ID2,]
  tempt2<-rbind(tempt2,c)   
}

# for dominance in Site Hadanberg
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_Dom_SH"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  d<- pred.list[[trait]][ID2,]
  tempt3<-rbind(tempt3,d) 
}  

# for dominance in Site Vindeln
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_Dom_SV"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  e<- pred.list[[trait]][ID2,]
  tempt4<-rbind(tempt4,e) 
}  

tempt<-cbind(tempt, tempt1, tempt2, tempt3, tempt4)
# SV: site vindeln, SH: site Hadanberg
# SV: site vindeln, SH: site Hadanberg
colnames(tempt)<-c("solution_add_main","std error_add_main", "z ratio_add_main", "solution_add_GXE_SH","std error_add_GXE_SH", "z ratio_add_GXE_SH","solution_add_GXE_SV","std error_add_GXE_SV", "z ratio_add_GXE_SV", "solution_dom_SH","std error_dom_SH", "z ratio_dom_SH", "solution_dom_SV","std error_dom_SV", "z ratio_dom_SV") 
#n1<-length(diag_G)


preds<-cbind(tempt, new.ph.df[, c("Site","Genotype_id")], as.matrix(A_inbreeding[56:length(A_inbreeding)])+1)
colnames(preds)<-c("solution_add_main","std error_add_main", "z ratio_add_main", "solution_add_GXE_SH","std error_add_GXE_SH", "z ratio_add_GXE_SH","solution_add_GXE_SV","std error_add_GXE_SV", "z ratio_add_GXE_SV", "solution_dom_SH","std error_dom_SH", "z ratio_dom_SH", "solution_dom_SV","std error_dom_SV", "z ratio_dom_SV", "Site","Genotype_id","inbreeding" ) 

# this script is not right. 
# preds.obs = merge(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds), by  = 0); dim(preds.obs) # here is noly 1381 
preds.obs <- cbind(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds)) 

# merge pred.obs with varcomps by fold to get the variance components for each fold
preds.obs <- merge(preds.obs, vc.mu, by = "fold")
head(preds.obs)
preds.obs$pred_add_main = preds.obs$mu + preds.obs$solution_add_main
preds.obs$pred_add_GXE_SH<-preds.obs$mu+ preds.obs$solution_add_main +preds.obs$solution_add_GXE_SH

# Compute reliabilities using Mrode's formula
names(preds.obs)[names(preds.obs) == "std error_add_main"] = "se_main"
names(preds.obs)[names(preds.obs) == "std error_add_GXE_SH"] = "se_GXE_SH"

preds.obs$reliability_main<-with(preds.obs,1-se_main^2/(Vadd_main*inbreeding))
# no meaning to do this????
preds.obs$reliability_GXE_SH<-with(preds.obs,1-se_GXE_SH^2/(Vadd_GXE_SH*inbreeding)) 


#  preds.obs = preds.obs[order(preds.obs$Site),]

# use round function to round numeric columns 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

preds.obs <- round_df(preds.obs, digits=3) 
# head(preds.obs); dim(preds.obs)
#  View(preds.obs)

# accuracy?? , here i think  everyone use Corr(AEBV, EBV), accuracy
# one paper in GS, use reliability in ABLUP??????? 
mean(preds.obs$reliability_main)# 0.1977
sqrt(mean(preds.obs$reliability_main)) # 0.44

mean(preds.obs$reliability_GXE_SH)# 0.1977
sqrt(mean(preds.obs$reliability_GXE_SH)) # 0.44


#*******************************************
#*************** predicitive ability*******

# additive
PredAbi = with(preds.obs, cor(adj_Hjd_17, pred)) 
temp<-ddply(preds.obs, .(fold), func_PreAbi)
SE_PredAbi<-sd(temp$COR)/sqrt(10)

pred.obs_backup3<-preds.obs


#****************** accuracy*******************
Accuracy=with(pred.obs_backup3, cor(EBV_Hjd_17, pred_add_main))
temp<-ddply(preds.obs, .(fold), func)
SE_Accuracy<-sd(temp$COR)/sqrt(10)
a3<-round(mean(temp$COR),3);   b3<-round(SE_Accuracy,3)

###################______________________________________________________ ABLUP+A

result.list = list()
pred.list = list()

# from the column which is repeated
for (trait in names(new.ph.df)[(jump_Nocolumn+1):(jump_Nocolumn+nfold)])
{
  asr <- asreml(fixed = get(trait) ~ 1, 
                random = ~ ped(Genotype_id) + diag(Site):ped(Genotype_id) ,
                ginverse=list(Genotype_id=pedinv),
                rcov= ~ at(Site):units, 
                data=new.ph.df, 
                na.method.X='include',
                workspace=20e8 ) # 32G rum, maximum space is 40e8, here i only use half of them + at(Site,1):giv(Genotype_idd) ,Genotype_idd= Gd.giv
  result.list[[trait]] = asr
  pred.list[[trait]] = summary(asr, all = T)$coef.rand
  rm(list = "asr")
}

# asr1 <- asreml(fixed = get(trait) ~ 1, 
#               random = ~ ped(Genotype_id) + diag(Site):ped(Genotype_id) ,
#               ginverse=list(Genotype_id=pedinv),
#               rcov= ~ at(Site):units, 
#               data=new.ph.df, 
#               na.method.X='include',
#               workspace=20e8 ) # 32G rum, maximum space is 40e8
# 
# summary(asr1)

########### Summarize cross-validations ###########################

# extract the variance components estimates from each fold
get.varcomps = function(c){
  #c is a component of result.list
  vcs = summary(c)$varcomp$component
  names(vcs) =  c( "Vadd_main", "Vadd_GXE_SH","Vadd_GXE_SV",  "Vres_SV", "Vres_SH")
  return(vcs)
}
varcomps = sapply(result.list, FUN = get.varcomps)
varcomps = as.data.frame(t(varcomps))
varcomps$fold = as.numeric(gsub("Y", "", rownames(varcomps)))

#get the mean for each estimation (training) set for each fold to use as mean of predictions
mu = colMeans(new.ph.df[, marker],na.rm = T)
mu = as.data.frame(mu)
mu$fold = as.numeric(gsub("Y", "", rownames(mu)))

vc.mu = merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
#************* we should change it if we want to get at(site,handanberg)
#******************* giv(Genotype_idd)_
new.ph.df$ID2_main<-paste("ped(Genotype_id)_", new.ph.df$Genotype_id,sep="") 
new.ph.df$ID2_add_GXE_SV<-paste("Site_Vinden:ped(Genotype_id)_", new.ph.df$Genotype_id,sep="")  # site Vindeln
new.ph.df$ID2_add_GXE_SH<-paste("Site_Hadanberg:ped(Genotype_id)_", new.ph.df$Genotype_id,sep="")  # Hadanberg


dim(new.ph.df)
tempt<-data.frame()
tempt1<-data.frame()
tempt2<-data.frame()
# for additive  in main
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_main"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  a<- pred.list[[trait]][ID2,]
  tempt<-rbind(tempt,a) 
}
# for GXE_additive  in Hadanberg
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_add_GXE_SH"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  b<- pred.list[[trait]][ID2,]
  tempt1<-rbind(tempt1,b)   
}
# for GXE_additive  in Vindeln
for (i in 1:dim(new.ph.df)[1]){
  r<-new.ph.df[i,]
  ID2 = r["ID2_add_GXE_SV"]
  ID2<-drop.levels(ID2[1,1])
  trait = paste0("Y", r["fold"])
  c<- pred.list[[trait]][ID2,]
  tempt2<-rbind(tempt2,c)   
}

tempt<-cbind(tempt, tempt1, tempt2)
# SV: site vindeln, SH: site Hadanberg
colnames(tempt)<-c("solution_add_main","std error_add_main", "z ratio_add_main", "solution_add_GXE_SH","std error_add_GXE_SH", "z ratio_add_GXE_SH", "solution_add_GXE_SV","std error_add_GXE_SV", "z ratio_add_GXE_SV") 
#n1<-length(diag_G)


preds<-cbind(tempt, new.ph.df[, c("Site","Genotype_id")], as.matrix(A_inbreeding[56:length(A_inbreeding)])+1)
colnames(preds)<-c("solution_add_main","std error_add_main", "z ratio_add_main", "solution_add_GXE_SH","std error_add_GXE_SH", "z ratio_add_GXE_SH", "solution_add_GXE_SV","std error_add_GXE_SV", "z ratio_add_GXE_SV",  "Site","Genotype_id","inbreeding" ) 

# this script is not right. 
# preds.obs = merge(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds), by  = 0); dim(preds.obs) # here is noly 1381 
preds.obs <- cbind(new.ph.df[,c("adj_Hjd_17","EBV_Hjd_17", "fold")], as.data.frame(preds)) 

# merge pred.obs with varcomps by fold to get the variance components for each fold
preds.obs <- merge(preds.obs, vc.mu, by = "fold")
head(preds.obs)
preds.obs$pred_add_main = preds.obs$mu + preds.obs$solution_add_main
preds.obs$pred_add_GXE_SH<-preds.obs$mu+ preds.obs$solution_add_main +preds.obs$solution_add_GXE_SH

# Compute reliabilities using Mrode's formula
names(preds.obs)[names(preds.obs) == "std error_add_main"] = "se_main"
names(preds.obs)[names(preds.obs) == "std error_add_GXE_SH"] = "se_GXE_SH"

preds.obs$reliability_main<-with(preds.obs,1-se_main^2/(Vadd_main*inbreeding))
# no meaning to do this????
preds.obs$reliability_GXE_SH<-with(preds.obs,1-se_GXE_SH^2/(Vadd_GXE_SH*inbreeding)) 


#  preds.obs = preds.obs[order(preds.obs$Site),]

# use round function to round numeric columns 
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

preds.obs <- round_df(preds.obs, digits=3) 
# head(preds.obs); dim(preds.obs)
#  View(preds.obs)

# accuracy?? , here i think  everyone use Corr(AEBV, EBV), accuracy
# one paper in GS, use reliability in ABLUP??????? 
mean(preds.obs$reliability_main)# 0.1977
sqrt(mean(preds.obs$reliability_main)) # 0.44

mean(preds.obs$reliability_GXE_SH)# 0.1977
sqrt(mean(preds.obs$reliability_GXE_SH)) # 0.44


#*******************************************
#*************** predicitive ability*******

# additive
PredAbi = with(preds.obs, cor(adj_Hjd_17, pred)) 
temp<-ddply(preds.obs, .(fold), func_PreAbi)
SE_PredAbi<-sd(temp$COR)/sqrt(10)

pred.obs_backup4<-preds.obs


#****************** accuracy*******************
Accuracy=with(pred.obs_backup4, cor(EBV_Hjd_17, pred_add_main))
temp<-ddply(preds.obs, .(fold), func)
SE_Accuracy<-sd(temp$COR)/sqrt(10)
a4<-round(mean(temp$COR),3);   b4<-round(SE_Accuracy,3)





