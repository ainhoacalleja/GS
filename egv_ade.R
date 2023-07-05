# Remove everything in the working environment.
rm(list=ls())

#Setting the working directory
setwd("~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/non_additive/EGV_traits/EGB_GBLUP_ADE")
#save.image("EGV_ade.RData")
#load("EGV_ade.RData")

####################################################################################################
########################### DENSITY TRAIT #######################

# Read in the data file and call it 'dat'.
Den1=read.table(file='pine3_ad_den.sln',header=T)
head(Den1)
tail(Den1)
str(Den1)

# Get the intercept term value from ASReml output file
# and call it 'mu'
Den_mu.idx=which(Den1$Model_Term=='mu')
Den1_mu=Den1$Effect[Den_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [dia1_add]
Den1_egv=Den1$Model_Term=='grm1(ID)'
Den1_ID=Den1$Level[Den1_egv] # select the levels
Den1_add=Den1$Effect[Den1_egv] 

####get the dominance values
Den1_egv2=Den1$Model_Term=='grm2(ID)'
Den1_ID=Den1$Level[Den1_egv2] # select the levels
Den1_dom=Den1$Effect[Den1_egv2] 

####get the epistatic values 1
Den1_egv3=Den1$Model_Term=='grm3(ID)'
Den1_ID=Den1$Level[Den1_egv3] # select the levels
Den1_ep1=Den1$Effect[Den1_egv3] 

####get the epistatic values 2
Den1_egv4=Den1$Model_Term=='grm4(ID)'
Den1_ID=Den1$Level[Den1_egv4] # select the levels
Den1_ep2=Den1$Effect[Den1_egv4] 

####get the epistatic values 3
Den1_egv5=Den1$Model_Term=='grm5(ID)'
Den1_ID=Den1$Level[Den1_egv5] # select the levels
Den1_ep3=Den1$Effect[Den1_egv5] 

###Estimate the EGV and EBV for den
EGV_den= Den1_mu + Den1_add + Den1_dom + Den1_ep1 + Den1_ep2 + Den1_ep3
EBV_den= Den1_mu + Den1_add 

################################################################


####################################################################################################
########################### DBH1 TRAIT #######################

# Read in the data file and call it 'dat'.
DBH1=read.table(file='pine3_dbh30.sln',header=T)
head(DBH1)
tail(DBH1)
str(DBH1)

# Get the intercept term value from ASReml output file
# and call it 'mu'
DBH1_mu.idx=which(DBH1$Model_Term=='mu')
DBH1_mu=DBH1$Effect[DBH1_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [dbh1_egv]
DBH1_egv=DBH1$Model_Term=='grm1(ID)'
DBH1_ID=DBH1$Level[DBH1_egv] # select the levels
DBH1_add=DBH1$Effect[DBH1_egv] 

####get the dominance values
DBH1_egv2=DBH1$Model_Term=='grm2(ID)'
DBH1_ID=DBH1$Level[DBH1_egv2] # select the levels
DBH1_dom=DBH1$Effect[DBH1_egv2] 

####get the epistatic values 1
DBH1_egv3=DBH1$Model_Term=='grm3(ID)'
DBH1_ID=DBH1$Level[DBH1_egv3] # select the levels
DBH1_ep1=DBH1$Effect[DBH1_egv3] 

####get the epistatic values 2
DBH1_egv4=DBH1$Model_Term=='grm4(ID)'
DBH1_ID=DBH1$Level[DBH1_egv4] # select the levels
DBH1_ep2=DBH1$Effect[DBH1_egv4] 

####get the epistatic values 3
DBH1_egv5=DBH1$Model_Term=='grm5(ID)'
DBH1_ID=DBH1$Level[DBH1_egv5] # select the levels
DBH1_ep3=DBH1$Effect[DBH1_egv5] 

###Estimate the EGV and EBV for den
EGV_DBH1= DBH1_mu + DBH1_add + DBH1_dom + DBH1_ep1 + DBH1_ep2 + DBH1_ep3
EBV_DBH1= DBH1_mu + DBH1_add

####################################################################################################
########################### DBH2 TRAIT #######################

# Read in the data file and call it 'dat'.
DBH2=read.table(file='pine3_dbh36.sln',header=T)
head(DBH2)
tail(DBH2)
str(DBH2)

# Get the intercept term value from ASReml output file
# and call it 'mu'
DBH2_mu.idx=which(DBH2$Model_Term=='mu')
DBH2_mu=DBH2$Effect[DBH2_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [dbh1_egv]
DBH2_egv=DBH2$Model_Term=='grm1(ID)'
DBH2_ID=DBH2$Level[DBH2_egv] # select the levels
DBH2_add=DBH2$Effect[DBH2_egv] 

####get the dominance values
DBH2_egv2=DBH2$Model_Term=='grm2(ID)'
DBH2_ID=DBH2$Level[DBH2_egv2] # select the levels
DBH2_dom=DBH2$Effect[DBH2_egv2] 

####get the epistatic values 1
DBH2_egv3=DBH2$Model_Term=='grm3(ID)'
DBH2_ID=DBH2$Level[DBH2_egv3] # select the levels
DBH2_ep1=DBH2$Effect[DBH2_egv3] 

####get the epistatic values 2
DBH2_egv4=DBH2$Model_Term=='grm4(ID)'
DBH2_ID=DBH2$Level[DBH2_egv4] # select the levels
DBH2_ep2=DBH2$Effect[DBH2_egv4] 

####get the epistatic values 3
DBH2_egv5=DBH2$Model_Term=='grm5(ID)'
DBH2_ID=DBH2$Level[DBH2_egv5] # select the levels
DBH2_ep3=DBH2$Effect[DBH2_egv5] 

###Estimate the EGV and EBV for den
EGV_DBH2= DBH2_mu + DBH2_add + DBH2_dom + DBH2_ep1 + DBH2_ep2 + DBH2_ep3
EBV_DBH2= DBH2_mu + DBH2_add

####################################################################################################
########################### HT1 TRAIT #######################

# Read in the data file and call it 'Ht1'.
Ht1=read.table(file='pine3_ht10.sln',header=T)
head(Ht1)
tail(Ht1)
str(Ht1)

# Get the intercept term value from ASReml output file
# and call it 'mu'
Ht1_mu.idx=which(Ht1$Model_Term=='mu')
Ht1_mu=Ht1$Effect[Ht1_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [Ht1_egv]
Ht1_egv=Ht1$Model_Term=='grm1(ID)'
Ht1_ID=Ht1$Level[Ht1_egv] # select the levels
Ht1_add=Ht1$Effect[Ht1_egv] 

####get the dominance values
Ht1_egv2=Ht1$Model_Term=='grm2(ID)'
Ht1_ID=Ht1$Level[Ht1_egv2] # select the levels
Ht1_dom=Ht1$Effect[Ht1_egv2] 

####get the epistatic values 1
Ht1_egv3=Ht1$Model_Term=='grm3(ID)'
Ht1_ID=Ht1$Level[Ht1_egv3] # select the levels
Ht1_ep1=Ht1$Effect[Ht1_egv3] 

####get the epistatic values 2
Ht1_egv4=Ht1$Model_Term=='grm4(ID)'
Ht1_ID=Ht1$Level[Ht1_egv4] # select the levels
Ht1_ep2=Ht1$Effect[Ht1_egv4] 

####get the epistatic values 3
Ht1_egv5=Ht1$Model_Term=='grm5(ID)'
Ht1_ID=Ht1$Level[Ht1_egv5] # select the levels
Ht1_ep3=Ht1$Effect[Ht1_egv5] 

###Estimate the EGV and EBV for den
EGV_Ht1= Ht1_mu + Ht1_add + Ht1_dom + Ht1_ep1 + Ht1_ep2 + Ht1_ep3
EBV_Ht1= Ht1_mu + Ht1_add

####################################################################################################
########################### HT2 TRAIT #######################

# Read in the data file and call it 'Ht2'.
Ht2=read.table(file='pine3_ht30.sln',header=T)
head(Ht2)
tail(Ht2)
str(Ht2)

# Get the intercept term value from ASReml output file
# and call it 'mu'
Ht2_mu.idx=which(Ht2$Model_Term=='mu')
Ht2_mu=Ht2$Effect[Ht2_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [Ht2_egv]
Ht2_egv=Ht2$Model_Term=='grm1(ID)'
Ht2_ID=Ht2$Level[Ht2_egv] # select the levels
Ht2_add=Ht2$Effect[Ht2_egv] 

####get the dominance values
Ht2_egv2=Ht2$Model_Term=='grm2(ID)'
Ht2_ID=Ht2$Level[Ht2_egv2] # select the levels
Ht2_dom=Ht2$Effect[Ht2_egv2] 

####get the epistatic values 1
Ht2_egv3=Ht2$Model_Term=='grm3(ID)'
Ht2_ID=Ht2$Level[Ht2_egv3] # select the levels
Ht2_ep1=Ht2$Effect[Ht2_egv3] 

####get the epistatic values 2
Ht2_egv4=Ht2$Model_Term=='grm4(ID)'
Ht2_ID=Ht2$Level[Ht2_egv4] # select the levels
Ht2_ep2=Ht2$Effect[Ht2_egv4] 

####get the epistatic values 3
Ht2_egv5=Ht2$Model_Term=='grm5(ID)'
Ht2_ID=Ht2$Level[Ht2_egv5] # select the levels
Ht2_ep3=Ht2$Effect[Ht2_egv5] 

###Estimate the EGV and EBV for den
EGV_Ht2= Ht2_mu + Ht2_add + Ht2_dom + Ht2_ep1 + Ht2_ep2 + Ht2_ep3
EBV_Ht2= Ht2_mu + Ht2_add

####################################################################################################
########################### MFA TRAIT #######################

# Read in the data file and call it 'MFA'.
MFA=read.table(file='pine3_MFA.sln',header=T)
head(MFA)
tail(MFA)
str(MFA)

# Get the intercept term value from ASReml output file
# and call it 'mu'
MFA_mu.idx=which(MFA$Model_Term=='mu')
MFA_mu=MFA$Effect[MFA_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [MFA_egv]
MFA_egv=MFA$Model_Term=='grm1(ID)'
MFA_ID=MFA$Level[MFA_egv] # select the levels
MFA_add=MFA$Effect[MFA_egv] 

####get the dominance values
MFA_egv2=MFA$Model_Term=='grm2(ID)'
MFA_ID=MFA$Level[MFA_egv2] # select the levels
MFA_dom=MFA$Effect[MFA_egv2] 

####get the epistatic values 1
MFA_egv3=MFA$Model_Term=='grm3(ID)'
MFA_ID=MFA$Level[MFA_egv3] # select the levels
MFA_ep1=MFA$Effect[MFA_egv3] 

####get the epistatic values 2
MFA_egv4=MFA$Model_Term=='grm4(ID)'
MFA_ID=MFA$Level[MFA_egv4] # select the levels
MFA_ep2=MFA$Effect[MFA_egv4] 

####get the epistatic values 3
MFA_egv5=MFA$Model_Term=='grm5(ID)'
MFA_ID=MFA$Level[MFA_egv5] # select the levels
MFA_ep3=MFA$Effect[MFA_egv5] 

###Estimate the EGV and EBV for den
EGV_MFA= MFA_mu + MFA_add + MFA_dom + MFA_ep1 + MFA_ep2 + MFA_ep3
EBV_MFA= MFA_mu + MFA_add

####################################################################################################
########################### MOEd TRAIT #######################

# Read in the data file and call it 'MOEd'.
MOEd=read.table(file='pine3_MOEd.sln',header=T)
head(MOEd)
tail(MOEd)
str(MOEd)

# Get the intercept term value from ASReml output file
# and call it 'mu'
MOEd_mu.idx=which(MOEd$Model_Term=='mu')
MOEd_mu=MOEd$Effect[MOEd_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [MOEd_egv]
MOEd_egv=MOEd$Model_Term=='grm1(ID)'
MOEd_ID=MOEd$Level[MOEd_egv] # select the levels
MOEd_add=MOEd$Effect[MOEd_egv] 

####get the dominance values
MOEd_egv2=MOEd$Model_Term=='grm2(ID)'
MOEd_ID=MOEd$Level[MOEd_egv2] # select the levels
MOEd_dom=MOEd$Effect[MOEd_egv2] 

####get the epistatic values 1
MOEd_egv3=MOEd$Model_Term=='grm3(ID)'
MOEd_ID=MOEd$Level[MOEd_egv3] # select the levels
MOEd_ep1=MOEd$Effect[MOEd_egv3] 

####get the epistatic values 2
MOEd_egv4=MOEd$Model_Term=='grm4(ID)'
MOEd_ID=MOEd$Level[MOEd_egv4] # select the levels
MOEd_ep2=MOEd$Effect[MOEd_egv4] 

####get the epistatic values 3
MOEd_egv5=MOEd$Model_Term=='grm5(ID)'
MOEd_ID=MOEd$Level[MOEd_egv5] # select the levels
MOEd_ep3=MOEd$Effect[MOEd_egv5] 

###Estimate the EGV and EBV for den
EGV_MOEd= MOEd_mu + MOEd_add + MOEd_dom + MOEd_ep1 + MOEd_ep2 + MOEd_ep3
EBV_MOEd= MOEd_mu + MOEd_add

####################################################################################################
########################### MOEs TRAIT #######################

# Read in the data file and call it 'MOEs'.
MOEs=read.table(file='pine3_MOEs.sln',header=T)
head(MOEs)
tail(MOEs)
str(MOEs)

# Get the intercept term value from ASReml output file
# and call it 'mu'
MOEs_mu.idx=which(MOEs$Model_Term=='mu')
MOEs_mu=MOEs$Effect[MOEs_mu.idx]


# Get the additive  values from ASReml output file
# and store them in array [MOEs_egv]
MOEs_egv=MOEs$Model_Term=='grm1(ID)'
MOEs_ID=MOEs$Level[MOEs_egv] # select the levels
MOEs_add=MOEs$Effect[MOEs_egv] 

####get the dominance values
MOEs_egv2=MOEs$Model_Term=='grm2(ID)'
MOEs_ID=MOEs$Level[MOEs_egv2] # select the levels
MOEs_dom=MOEs$Effect[MOEs_egv2] 

####get the epistatic values 1
MOEs_egv3=MOEs$Model_Term=='grm3(ID)'
MOEs_ID=MOEs$Level[MOEs_egv3] # select the levels
MOEs_ep1=MOEs$Effect[MOEs_egv3] 

####get the epistatic values 2
MOEs_egv4=MOEs$Model_Term=='grm4(ID)'
MOEs_ID=MOEs$Level[MOEs_egv4] # select the levels
MOEs_ep2=MOEs$Effect[MOEs_egv4] 

#MOEs_ID=as.data.frame(MOEs$Level[MOEs_egv4])

####get the epistatic values 3
MOEs_egv5=MOEs$Model_Term=='grm5(ID)'
MOEs_ID=MOEs$Level[MOEs_egv5] # select the levels
MOEs_ep3=MOEs$Effect[MOEs_egv5] 


###Estimate the EGV and EBV for den
EGV_MOEs= MOEs_mu + MOEs_add + MOEs_dom + MOEs_ep1 + MOEs_ep2 + MOEs_ep3
EBV_MOEs= MOEs_mu + MOEs_add

solution_EGV <- data.frame(cbind(MOEs_ID, EGV_Ht1, EGV_Ht2, EGV_DBH1, EGV_DBH2, EGV_MFA, EGV_MOEs, EGV_den, EGV_MOEd ))
solution_EBV <- data.frame(cbind(MOEs_ID, EBV_Ht1, EBV_Ht2, EBV_DBH1, EBV_DBH2, EBV_MFA, EBV_MOEs, EBV_den, EBV_MOEd))

write.csv(solution_EGV,file =  "GBLUP_ADE_EGV.csv", row.names = FALSE)
write.csv(solution_EBV, file = "GBLUP_ADE_EBV.csv", row.names = FALSE)

#save.image("EGV_ade.RData")