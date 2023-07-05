#Calculate adjusted phenotypes (BVs) 
#Remove everything from the memory
rm(list=ls())

#path of the working directory
path="~/Documents/PhD_UPSC/Project_3_Grundtjärn/analysis/EBV_extract"
setwd(path)
save.image("EBVgrowth.RData")

#Read sln datafiles created with ASReml standalone
dbh30=read.table(file='pine_S11_dbh30.sln',header=T)
head(dbh30)
dbh36=read.table(file='pine_S11_dbh36.sln', header=T)
head(dbh36)
ht10=read.table(file='pine_S11_ht10.sln',header=T)
head(ht10)
ht30=read.table(file='pine_S11_ht30.sln', header=T)
head(ht30)
tail(ht30)

#get intercept value. Call it "mu_traitname"

#Diameter 30 (spatilly adjusted previously)
mu.dbh30.idx=which(dbh30$Model_Term=='mu')
mu.dbh30=dbh30$Effect[mu.dbh30.idx]
#Get BVs from SLN file and store them in array 'fam_traitname_idx'
fam.dbh30.idx=dbh30$Model_Term=='nrm(ID)'
ID.dbh30=dbh30$Level[fam.dbh30.idx]
BLUP_d30=dbh30$Effect[fam.dbh30.idx]
EBV_d30=dbh30$Effect[fam.dbh30.idx]+mu.dbh30
head(EBV_d30)
head(BLUP_d30)

#Diameter 36 (spatially adjusted previously)
mu.dbh36.idx=which(dbh36$Model_Term=='mu')
mu.dbh36=dbh36$Effect[mu.dbh36.idx]
#Get BVs from SLN file and store them in array 
fam.dbh36.idx=dbh36$Model_Term=='nrm(ID)'
ID.dbh36=dbh36$Level[fam.dbh36.idx]
EBV_d36=dbh36$Effect[fam.dbh36.idx]+mu.dbh36
BLUP_d36=dbh36$Effect[fam.dbh36.idx]
head(EBV_d36)
head(BLUP_d36)
save.image("EBVgrowth.RData")

#Height 10 (spatially adjusted previously)
mu.ht10.idx=which(ht10$Model_Term=='mu')
mu.ht10=ht10$Effect[mu.ht10.idx]
#Get BV and store them
fam.ht10.idx=ht10$Model_Term=='nrm(ID)'
ID.ht10=ht10$Level[fam.ht10.idx]
EBV_h10=ht10$Effect[fam.ht10.idx]+mu.ht10
BLUP_h10=ht10$Effect[fam.ht10.idx]
head(EBV_h10)
head(BLUP_h10)

#Height 30 (spatially adjusted previously)
mu.ht30.idx=which(ht30$Model_Term=='mu')
mu.ht30=ht30$Effect[mu.ht30.idx] 
#Get BVs and store
fam.ht30.idx=ht30$Model_Term=='nrm(ID)'
ID.ht30=ht30$Level[fam.ht30.idx]
EBV_h30=ht30$Effect[fam.ht30.idx]+mu.ht30
BLUP_h30=ht30$Effect[fam.ht30.idx]
head(EBV_h30)
head(BLUP_h30)


#Chech if ID are equal for all the data = TRUE, they are equal
identical(ID.dbh30,ID.dbh36)
identical(ID.dbh30, ID.ht10)
identical(ID.dbh30, ID.ht30)
save.image("EBVgrowth.RData")


#merge all in a dataframe
solution<-data.frame(ID=ID.dbh30, EBV_d30=EBV_d30, BLUP_d30=BLUP_d30, EBV_d36=EBV_d36, BLUP_d36=BLUP_d36,
                     EBV_h10=EBV_h10, BLUP_h10=BLUP_h10, EBV_h30=EBV_h30, BLUP_h30=BLUP_h30)
head(solution,10)
tail(solution,10)

write.table(solution, "EBV_growth.csv",sep = ",", row.names =FALSE )
write.csv(solution, file="EBV_growth1.csv",row.names = TRUE)
save.image("EBVgrowth.RData")
