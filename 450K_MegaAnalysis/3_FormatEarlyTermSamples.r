#Create a consistent format for all RGsets you will use in the early term discovery set
#Then merge them togther to make one RGset that you will put on AWS

library(minfi)
###################################
#Dataset 1: GSE108567
###################################
GEOid <- "GSE108567"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset1 <- RGset[,which(pd$GA < 37)]
pd1 <- pd[which(pd$GA < 37),]

newpd <- DataFrame(Sex = pd1$Sex, Study = rep(GEOid, nrow(pd1)), GA = as.numeric(pd1$GA)) 
newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
pData(RGset1) <- newpd

rm(list=setdiff(ls(), "RGset1")) #dim: 622399 14

###################################
#Dataset 2: GSE100197
###################################
GEOid <- "GSE100197"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset2 <- RGset[,which(pd$Group == "PreT")]
pd2 <- pd[which(pd$Group == "PreT"),]

newpd <- DataFrame(Sex = pd2$Sex, Study = rep(GEOid, nrow(pd2)), GA = as.numeric(pd2$GA)) 
newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
pData(RGset2) <- newpd #dim: 622399 24 

rm(list=setdiff(ls(), c("RGset1","RGset2")))


###################################
#Dataset 3: GSE98224
###################################
GEOid <- "GSE98224"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset3 <- RGset[,which(pd$Title == "Cont-preT-AGA")]
pd3 <- pd[which(pd$Title == "Cont-preT-AGA"),]

newpd <- DataFrame(Sex = pd3$gender, Study = rep(GEOid, nrow(pd3)), GA = as.numeric(NA)) 
#newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
newpd$Sex <- as.character(newpd$Sex)
pData(RGset3) <- newpd #dim: 622399 5

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3")))

###################################
#Dataset 4: GSE106089
###################################
GEOid <- "GSE106089"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset4 <- RGset[,which(pd$Bacteria == "None")]
pd4 <- pd[which(pd$Bacteria == "None"),]

newpd <- DataFrame(Sex = pd4$Sex, Study = rep(GEOid, nrow(pd4)), GA = as.numeric(NA)) 
newpd$Sex <- ifelse(newpd$Sex == "male", "M", "F") 
newpd$Sex <- as.character(newpd$Sex)
pData(RGset4) <- newpd #dim: 622399 46

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4")))

###################################
#Dataset 5: GSE103413
###################################
GEOid <- "GSE103413"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset5 <- RGset[,which(pd$Tissue == "chorionic villus")]
pd5 <- pd[which(pd$Tissue == "chorionic villus"),]

newpd <- DataFrame(Sex = pd5$Sex, Study = rep(GEOid, nrow(pd5)), GA = as.numeric(NA)) 
newpd$Sex <- ifelse(newpd$Sex == "male", "M", "F") 
newpd$Sex <- as.character(newpd$Sex)
pData(RGset5) <- newpd #dim: 622399 2

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4","RGset5")))

###################################
#Dataset 6: GSE98938
###################################
GEOid <- "GSE98938"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset6 <- RGset[,which(pd$Tissue %in% c("chorion","villi","trophoblast") & pd$GA == "second trimester")]
pd6 <- pd[which(pd$Tissue %in% c("chorion","villi","trophoblast") & pd$GA == "second trimester"),]

newpd <- DataFrame(Sex = NA, Study = rep(GEOid, nrow(pd6)), GA = as.numeric(NA)) 
#newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
#newpd$Sex <- as.character(newpd$Sex)
pData(RGset6) <- newpd #dim: 622399 6

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4","RGset5","RGset6")))

###################################
#Dataset 7: GSE93208
###################################
GEOid <- "GSE93208"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset7 <- RGset

newpd <- DataFrame(Sex = as.character(NA), Study = rep(GEOid, nrow(pd)), GA = as.numeric(NA)) 
pData(RGset7) <- newpd

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4","RGset5","RGset6","RGset7")))

###################################
#Dataset 12: GSE69502
###################################
GEOid <- "GSE69502"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset12 <- RGset[,which(pd$Tissue == "chorionic villi" & pd$Phenotype == "control")]
pd12 <- pd[which(pd$Tissue == "chorionic villi" & pd$Phenotype == "control"),]

newpd <- DataFrame(Sex = pd12$Sex, Study = rep(GEOid, nrow(pd12)), GA = as.numeric(pd12$GA))
newpd$Sex <- ifelse(newpd$Sex == "male fetus", "M", "F") 
newpd$Sex <- as.character(newpd$Sex)
pData(RGset12) <- newpd

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4","RGset5","RGset6","RGset7","RGset12")))

###################################
#Dataset 14: GSE66210
###################################
GEOid <- "GSE66210"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset14 <- RGset[,which(pd$Tissue == "chorionic villus" & pd$Diagnosis == "normal pregnancy")]
pd14 <- pd[which(pd$Tissue == "chorionic villus" & pd$Diagnosis == "normal pregnancy"),]

newpd <- DataFrame(Sex = pd14$Sex, Study = rep(GEOid, nrow(pd14)), GA = as.numeric(NA))
pData(RGset14) <- newpd

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset4","RGset5","RGset6","RGset7",
	"RGset12","RGset14")))

################
RGset.early <- combine(RGset1,RGset2,RGset3,RGset4,RGset5,RGset6,RGset7,RGset12,RGset14)
#dim: 622399 145

dups <- duplicatesamples(RGset.early)

dups.ids <- (dups[[2]]) #length 12
matches <- dups[[1]]
remme <- c()
for (i in 1:nrow(matches)){
	mypair <- matches[i,]
	detp <- dups.ids[match(mypair,rownames(dups.ids)),1]
	pickme <- 1
	if(length(unique(detp))>1){
		pickme <- which(detp == max(detp))
	}
	remme <- c(remme,as.numeric(names(detp)[pickme]))
}

RGset.early.nodups <- RGset.early[,-remme]
pd.early.nodups <- pData(RGset.early.nodups)
save(RGset.early.nodups, file = "RGset.early.nodups.rda")
save(pd.early.nodups, file = "pd.early.nodups.rda")

#Classify into trimester when available

pd <- pd.early.nodups
pd$Trimester <- NA

pd$Trimester[which(pd$GA < 36)] <- "First"
pd$Trimester[which(pd$GA >= 14)] <- "Second"
pd$Trimester[which(pd$GA >= 27)] <- "Third"

#Following studies have no classification: GSE103413, GSE106089, GSE66210, 
#GSE93208, GSE98224, GSE98938

#Must go through and see if associated paper has info documented. 

######################
#GSE103413
GEOid <- "GSE103413"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pdstudy <- pData(RGset)

#Not in paper but predicted to be in second trimester. 
# > pdstudy$GApred[which(pdstudy$Tissue == "chorionic villus")]
# GSM2770903_9285451169_R01C02 GSM2770904_9285451169_R02C02 
#                     15.47176                     14.08281 
pd$Trimester[which(pd$Study == GEOid)] <- "Second"


#####################
#GSE106089
#"we invited women who gave birth before 28 weeks gestational age at one"
GEOid <- "GSE106089"
pd$Trimester[which(pd$Study == GEOid)] <- "Second"


#####################
#GSE66210
#GEO title says first trimester
GEOid <- "GSE66210"
pd$Trimester[which(pd$Study == GEOid)] <- "First"


#####################
#GSE93208
#GEO title says first trimester
GEOid <- "GSE93208"
pd$Trimester[which(pd$Study == GEOid)] <- "First"

#####################
#GSE98224
GEOid <- "GSE98224"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pdstudy <- pData(RGset)
pd3 <- pdstudy[which(pdstudy$Title == "Cont-preT-AGA"),]
#All reported as third term
# > pd3[,19]
# [1] 30 31 33 32 32
pd$Trimester[which(pd$Study == GEOid)] <- "Third"

#####################
#GSE98938
#samples labeled as second trimester
GEOid <- "GSE98938"
pd$Trimester[which(pd$Study == GEOid)] <- "Second"

# > table(pd$Trimester)

#  First Second  Third 
#     31     70     38 

pData(RGset.early.nodups) <- pd
pd.early.nodups <- pData(RGset.early.nodups)
save(RGset.early.nodups, file = "RGset.early.nodups.rda")
save(pd.early.nodups, file = "pd.early.nodups.rda")

