#Create a consistent format for all RGsets you will use in the term discovery set
#Then merge them togther to make one RGset that you will put on AWS

library(minfi)

###################################
#Dataset 1: GSE108567
###################################
GEOid <- "GSE108567"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset1 <- RGset[,which(pd$GA >= 37)]
pd1 <- pd[which(pd$GA >= 37),]

newpd <- DataFrame(Sex = pd1$Sex, Study = rep(GEOid, nrow(pd1))) 
newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
pData(RGset1) <- newpd

rm(list=setdiff(ls(), "RGset1")) #dim: 622399 45

###################################
#Dataset 2: GSE100197
###################################
GEOid <- "GSE100197"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset2 <- RGset[,which(pd$Group == "Term")]
pd2 <- pd[which(pd$Group == "Term"),]

newpd <- DataFrame(Sex = pd2$Sex, Study = rep(GEOid, nrow(pd2))) 
newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
pData(RGset2) <- newpd #dim: 622399 19 

rm(list=setdiff(ls(), c("RGset1","RGset2")))

###################################
#Dataset 3: GSE98224
###################################
GEOid <- "GSE98224"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset3 <- RGset[,which(pd$Title == "Cont-Ter-AGA")]
pd3 <- pd[which(pd$Title == "Cont-Ter-AGA"),]

newpd <- DataFrame(Sex = pd3$gender, Study = rep(GEOid, nrow(pd3))) 
#newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
newpd$Sex <- as.character(newpd$Sex)
pData(RGset3) <- newpd #dim: 622399 9 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3")))

###################################
#Dataset 6: GSE98938
###################################
GEOid <- "GSE98938"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset6 <- RGset[,which(pd$Tissue %in% c("chorion","villi","trophoblast") & pd$GA == "term")]
pd6 <- pd[which(pd$Tissue %in% c("chorion","villi","trophoblast") & pd$GA == "term"),]

newpd <- DataFrame(Sex = NA, Study = rep(GEOid, nrow(pd6))) 
#newpd$Sex <- ifelse(newpd$Sex == "MALE", "M", "F") 
#newpd$Sex <- as.character(newpd$Sex)
pData(RGset6) <- newpd #dim: 622399 6

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6"))) 

###################################
#Dataset 8: GSE71719
###################################
GEOid <- "GSE71719"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset8 <- RGset[,which(pd$DNAType == "bisulfite")]
pd8 <- pd[which(pd$DNAType == "bisulfite"),]

newpd <- DataFrame(Sex = pd8$Sex, Study = rep(GEOid, nrow(pd8))) 
newpd$Sex <- ifelse(newpd$Sex == "Male", "M", "F") 
pData(RGset8) <- newpd #dim: 622399 23 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6","RGset8"))) 

###################################
#Dataset 9: GSE71678
###################################
GEOid <- "GSE71678"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset9 <- RGset[,which(pd$GA >= 37)]
pd9 <- pd[which(pd$GA >= 37),]

newpd <- DataFrame(Sex = pd9$Sex, Study = rep(GEOid, nrow(pd9))) 
newpd$Sex <- ifelse(newpd$Sex == "Male", "M", "F") 
pData(RGset9) <- newpd #dim: 622399 323 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6","RGset8","RGset9"))) 

###################################
#Dataset 10: GSE75248
###################################
GEOid <- "GSE75248"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset10 <- RGset[,which(pd$bwg == "AGA")]
pd10 <- pd[which(pd$bwg == "AGA"),]

newpd <- DataFrame(Sex = pd10$gender, Study = rep(GEOid, nrow(pd10))) 
newpd$Sex <- ifelse(newpd$Sex == "Male", "M", "F") 
pData(RGset10) <- newpd #dim: 622399 174 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6","RGset8","RGset9","RGset10"))) 

###################################
#Dataset 11: GSE75196
###################################
GEOid <- "GSE75196"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset11 <- RGset[,which(pd$Pheno == "normal/healthy" & pd$GA >= 37)]
pd11 <- pd[which(pd$Pheno == "normal/healthy" & pd$GA >= 37),]

newpd <- DataFrame(Sex = pd11$Sex, Study = rep(GEOid, nrow(pd11))) 
newpd$Sex <- ifelse(newpd$Sex == "Male", "M", "F") 
pData(RGset11) <- newpd #dim: 622399 16 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6","RGset8","RGset9","RGset10","RGset11"))) 

###################################
#Dataset 13: GSE74738
###################################
GEOid <- "GSE74738"
load(paste0("~/Documents/UCSF/FetalSexMetaAnalysis/RawData/",GEOid,"/RGset.rda"))
pd <- pData(RGset)

RGset13 <- RGset[,which(pd$GA >= 37)]
pd13 <- pd[which(pd$GA >= 37),]

newpd <- DataFrame(Sex = pd13$Sex, Study = rep(GEOid, nrow(pd13))) 
newpd$Sex <- ifelse(newpd$Sex == "F", "F", "M") 
pData(RGset13) <- newpd #dim: 622399 22 

rm(list=setdiff(ls(), c("RGset1","RGset2","RGset3","RGset6","RGset8",
	"RGset9","RGset10","RGset11","RGset13"))) 

RGset.term <- combine(RGset1,RGset2,RGset3,RGset6,RGset8,RGset9,RGset10,
			RGset11,RGset13)

save(RGset.term, file = "RGset.term.rda")
pd.term <- pData(RGset.term)
save(pd.term, file = "pd.term.rda")

#Duplicate samples. 
duplicatesamples <- function(RGsetX, cutoff = 0.95){
  #FIND THE SNPS 
  snps <- getSnpBeta(RGsetX)
  
  #CLUSTER EACH SNP INTO ONE OF THREE GROUPS
  clusters <- matrix(rep(0, ncol(RGsetX) * nrow(snps)), nrow = nrow(snps))
  rownames(clusters) <- rownames(snps)
  colnames(clusters) <- colnames(snps)
  for(n in 1:nrow(snps)){
    kmeans1 <- kmeans(snps[n, ], 3)
    clusters[n, ] <- kmeans1$cluster
  }
  
  #PAIRWISE COMPARISON BETWEEN SNPS 
  similarity <- matrix(rep(0, ncol(RGsetX) * ncol(RGsetX)), ncol = ncol(RGsetX))
  rownames(similarity) <- colnames(clusters)
  colnames(similarity) <- colnames(clusters)
  for(s1 in 1:(ncol(RGsetX)-1)){
    for(s2 in (s1+1):ncol(RGsetX)){
      sim <- which(clusters[, s1] == clusters[, s2])
      percent <- length(sim) / nrow(snps)
      similarity[s1, s2] <- percent
      #similarity[s2, s1] <- percent
    }
  }
  simdim <- which(similarity >= cutoff, arr.ind = TRUE)
  duplicates <- simdim
  rownames(duplicates) <- NULL
  
  #NUMBER OF PROBES THAT HAVE P-VAL > 0.01
  dupvec <- unique(as.vector(duplicates))
  output <- duplicates
  if(length(dupvec)>0){
	  p <- detectionP(RGsetX[, dupvec], type = "m+u")
	  failed <- p > 0.01
	  numfailedpval <- matrix(rep(0, length(dupvec)), ncol = 1)
	  rownames(numfailedpval) <- dupvec
	  for(i in 1:length(dupvec)){
	    numfailedpval[i, 1] <- length(which(failed[, i] == TRUE))
	  }
	  output <- list(duplicates, numfailedpval)
	}
return(output)
}

dups <- duplicatesamples(RGset.term)
#Sift through dups and keep the smallest number 
#(try to minimize numbers coming from different batches if same exact idats )

dups.ids <- rownames(dups[[2]]) #length 56
matches <- dups[[1]]
remme <- c()
for (i in 1:length(dups.ids)){
	relrows <- which(matches[,1] == dups.ids[i])
	if(length(relrows) > 0)
		allsamps <- unique(c(matches[relrows,]))
		newrows <- unique(c(which(matches[,1] %in% allsamps |  matches[,2] %in% allsamps)))
		candidates <- unique(c(matches[newrows,]))	
		remme <- c(remme,candidates[-1])
}

remme <- unique(remme) #length 33
dups[[2]][-which(rownames(dups[[2]])%in%remme),]

RGset.term.nodups <- RGset.term[,-remme]
pd.term.nodups <- pData(RGset.term.nodups)
save(RGset.term.nodups, file = "RGset.term.nodups.rda")
save(pd.term.nodups, file = "pd.term.nodups.rda")
