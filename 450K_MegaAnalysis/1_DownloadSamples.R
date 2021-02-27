#Download the public datasets for use in fetal sex placenta DNAm
#meta-analysis. 

###################################
#Dataset 1: GSE108567
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE108567",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE108567",destdir=getwd())
mypheno<-(phenoData(mine$GSE108567_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE108567")
system("tar -xvf GSE108567_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),7))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-as.numeric(unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,4]<-unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,5]<-unlist(lapply(as.character(mypheno[,14]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,6]<-unlist(lapply(as.character(mypheno[,15]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,35]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,7]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","Tissue","GA","ProcessGroup","SentrixID","SentrixPos","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 2: GSE100197
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE100197",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE100197",destdir=getwd())
mypheno<-(phenoData(mine$GSE100197_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE100197")
system("tar -xvf GSE100197_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),9))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-(unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,4]<-as.numeric(unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,5]<-unlist(lapply(as.character(mypheno[,14]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,6]<-unlist(lapply(as.character(mypheno[,15]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,7]<-unlist(lapply(as.character(mypheno[,16]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,8]<-unlist(lapply(as.character(mypheno[,17]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,38]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,9]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("ID","Group","Sex","GA","Tissue","Plate","SentrixID","SentrixPos","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 3: GSE98224
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE98224",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE98224",destdir=getwd())
mypheno<-phenoData(mine[[1]])
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE98224")
system("tar -xvf GSE98224_RAW.tar")
system("gunzip *.idat.gz")

myvars <- grep("characteristics",colnames(mypheno))
pd <- lapply(myvars, function (x) {
	unlist(lapply(as.character(mypheno[,x]),function(y){strsplit(y,": ")[[1]][2]}))
	})
pd <- do.call("cbind",pd)
tmp <- unlist(lapply(as.character(mypheno[,59]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd <- data.frame(cbind(pd, unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))))
colnames(pd) <- c(unique(unlist(lapply(myvars, function (x) {
	unlist(lapply(as.character(mypheno[,x]),function(y){strsplit(y,": ")[[1]][1]}))
	}))),"Basename")
pd$Title <- unlist(lapply(as.character(mypheno[,1]),function(x){strsplit(x,",")[[1]][1]}))
pd[,3] <- as.numeric(as.character(pd[,3]))
pd[,4] <- as.numeric(as.character(pd[,4]))
pd[,10] <- as.numeric(as.character(pd[,10]))
pd[,11] <- as.numeric(as.character(pd[,11]))
pd[,12] <- as.numeric(as.character(pd[,12]))
pd[,13] <- as.numeric(as.character(pd[,13]))
pd[,14] <- as.numeric(as.character(pd[,14]))
pd[,19] <- as.numeric(as.character(pd[,19]))
pd[,20] <- as.numeric(as.character(pd[,20]))
pd[,22] <- as.numeric(as.character(pd[,22]))
pd[,23] <- as.numeric(as.character(pd[,23]))
pd[,24] <- as.numeric(as.character(pd[,24]))
pd[,25] <- as.numeric(as.character(pd[,25]))
pd[,26] <- as.numeric(as.character(pd[,26]))
pd$Basename <- as.character(pd$Basename)

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")
system("rm ../GPL6244.soft")

###################################
#Dataset 4: GSE106089
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE106089",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE106089",destdir=getwd())
mypheno<-(phenoData(mine$GSE106089_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE106089")
system("tar -xvf GSE106089_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),9))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-(unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,4]<-(unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,5]<-unlist(lapply(as.character(mypheno[,14]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,6]<-unlist(lapply(as.character(mypheno[,15]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,7]<-unlist(lapply(as.character(mypheno[,16]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,8]<-unlist(lapply(as.character(mypheno[,17]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,36]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,9]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Bacteria","Tissue","Sex","Race","Antiobiotics","VaginalInf","UTI","Csection","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 5: GSE103413
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE103413",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE103413",destdir=getwd())
mypheno<-(phenoData(mine$GSE103413_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE103413")
system("tar -xvf GSE103413_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),3))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,31]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,3]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","Tissue","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 6: GSE98938
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE98938",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE98938",destdir=getwd())
mypheno<-(phenoData(mine$GSE98938_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE98938")
system("tar -xvf GSE98938_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),3))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,32]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,3]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Tissue","GA","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 7: GSE93208: Nordor (1st Trimester, single cell type)
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE93208",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE93208",destdir=getwd())
mypheno<-(phenoData(mine$GSE93208_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE93208")
system("tar -xvf GSE93208_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),4))
pd[,1]<-unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,14]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-unlist(lapply(as.character(mypheno[,15]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,36]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,4]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("GA","DNAquant","Date","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 8: GSE71719
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE71719",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE71719",destdir=getwd())
mypheno<-(phenoData(mine$GSE71719_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE71719")
system("tar -xvf GSE71719_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),3))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,33]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,3]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","DNAType","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 9: GSE71678: Green et al (arsenic)...IDATs available
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE71678",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE71678",destdir=getwd())
mypheno<-(phenoData(mine$GSE71678_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE71678")
system("tar -xvf GSE71678_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),5))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-as.numeric(unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(as.character(x),": ")[[1]][2]})))
pd[,3]<-as.numeric(unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(as.character(x),": ")[[1]][2]})))
pd[,4]<-as.numeric(unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(as.character(x),": ")[[1]][2]})))
tmp<-unlist(lapply(as.character(mypheno[,36]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,5]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","MaternalAge","BWg","GA","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 10: GSE75248: Paquette et al (infant neurobehavioral outcomes)...IDATS available
###################################
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
library(minfi)
library(GEOquery)
getGEOSuppFiles("GSE75248",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE75248",destdir=getwd())
mypheno<-(phenoData(mine$GSE75248_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE75248")
system("tar -xvf GSE75248_RAW.tar")
system("gunzip *.idat.gz")

gender<-unlist(lapply(mypheno[,10],function(x){strsplit(sub(" ", ";", x), ";")[[1]][2]}))
gender<-factor(gender,levels=c("Male","Female"))
batch<-unlist(lapply(mypheno[,11],function(x){strsplit(sub(" ", ";", x), ";")[[1]][2]}))
batch<-factor(batch,levels=c("2011GROUP","2014GROUP"))
bwg<-unlist(lapply(mypheno[,12],function(x){strsplit(gsub(" ", ";", x), ";")[[1]][4]}))
bwg<-factor(bwg,levels=c("AGA","LGA","SGA"))
attn<-as.numeric(unlist(lapply(mypheno[,13],function(x){strsplit(sub(" ", ";", x), ";")[[1]][2]})))
qom<-as.numeric(unlist(lapply(mypheno[,14],function(x){strsplit(gsub(" ", ";", x), ";")[[1]][4]})))
arous<-as.numeric(unlist(lapply(mypheno[,15],function(x){strsplit(sub(" ", ";", x), ";")[[1]][2]})))
leth<-as.numeric(unlist(lapply(mypheno[,16],function(x){strsplit(sub(" ", ";", x), ";")[[1]][2]})))

vars<-data.frame(gender,batch,bwg,attn,qom,arous,leth)
filename<-unlist(lapply(as.character(mypheno[,36]),function(x){strsplit(gsub("/", ";", x), ";")[[1]][9]}))
filename<-unlist(lapply(filename,function(x){strsplit(x,"_G")[[1]][1]}))
slidename<-unlist(lapply(filename,function(x){strsplit(x,"_")[[1]][2]}))
arrayname<-unlist(lapply(filename,function(x){strsplit(x,"_")[[1]][3]}))
vars$Array<-arrayname
vars$Slide<-slidename
vars$Basename<-filename

RGset<-read.metharray.exp(getwd(),targets=vars,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")


###################################
#Dataset 11: GSE75196: Yeung et al (preeclampsia)
###################################

library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE75196",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE75196",destdir=getwd())
mypheno<-(phenoData(mine$GSE75196_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE75196")
system("tar -xvf GSE75196_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),4))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,"disease: ")[[1]][2]}))
#Hardcode gestational age because can't go as numeric with number there. 
pd[,3]<-c(36,32,32,35,38,36,37,35,39,38,39,38,39,38,40,39,39,39,38,38,39,40,38,40)
#pd[,3]<-as.numeric(unlist(lapply(mypheno[,12],function(x){strsplit(x,"gestation (wk): ")[[1]][2]})))
tmp<-unlist(lapply(as.character(mypheno[,32]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,4]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","Pheno","GA","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")


###################################
#Dataset 12: GSE69502: Price et al 2016 (NTB)
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE69502",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE69502",destdir=getwd())
mypheno<-(phenoData(mine$GSE69502_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE69502")
system("tar -xvf GSE69502_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),5))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,4]<-unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,37]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,5]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Phenotype","Sex","Tissue","GA","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 13: GSE74738: Hanna et al 2016 (Wendy Robinson imprinting thing)
###################################

library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE74738",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE74738",destdir=getwd())
mypheno<-(phenoData(mine$GSE74738_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE74738")
system("tar -xvf GSE74738_RAW.tar")
system("gunzip *.idat.gz")

mypheno<-mypheno[which(mypheno[,12]=="sample tissue: placental chorionic villi"),]
mypheno<-mypheno[which(mypheno[,10]=="status/group: control"),]

pd<-data.frame(matrix(0,nrow(mypheno),3))
pd[,1]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-as.numeric(unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]})))
tmp<-unlist(lapply(as.character(mypheno[,37]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,3]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","GA","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#Dataset 14: GSE66210
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/UCSF/FetalSexMetaAnalysis/RawData")
getGEOSuppFiles("GSE66210",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE66210",destdir=getwd())
mypheno<-(phenoData(mine$GSE66210_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")
rm(mine)
setwd("./GSE66210")
system("tar -xvf GSE66210_RAW.tar")
system("gunzip *.idat.gz")

pd<-data.frame(matrix(0,nrow(mypheno),4))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]}))
tmp<-unlist(lapply(as.character(mypheno[,33]),function(x){strsplit(x,"suppl/")[[1]][2]}))
pd[,4]<-unlist(lapply(tmp,function(x){strsplit(x,"_Grn")[[1]][1]}))
colnames(pd)<-c("Sex","Tissue","Diagnosis","Basename")

RGset<-read.metharray.exp(getwd(),targets=pd,verbose=TRUE)
save(RGset,file="RGset.rda")
system("rm *.idat")
system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")

###################################
#PLACENTA CELL TYPE DATASET: GSE159526
###################################
library(minfi)
library(GEOquery)
setwd("~/Documents/Misc/450k_FetalSex_Placenta/RawData")
getGEOSuppFiles("GSE159526",makeDirectory=TRUE,baseDir=getwd())
mine<-getGEO(GEO="GSE159526",destdir=getwd())

B <-  exprs(mine[[1]])

mypheno<-(phenoData(mine$GSE159526_series_matrix.txt.gz))
mypheno<-as(mypheno,"data.frame")

setwd("./GSE159526")
pd<-data.frame(matrix(0,nrow(mypheno),6))
pd[,1]<-unlist(lapply(as.character(mypheno[,10]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,2]<-unlist(lapply(as.character(mypheno[,11]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,3]<-as.numeric(unlist(lapply(as.character(mypheno[,12]),function(x){strsplit(x,": ")[[1]][2]})))
pd[,4]<-unlist(lapply(as.character(mypheno[,13]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,5]<-unlist(lapply(as.character(mypheno[,14]),function(x){strsplit(x,": ")[[1]][2]}))
pd[,6]<-unlist(lapply(as.character(mypheno[,15]),function(x){strsplit(x,": ")[[1]][2]}))
colnames(pd)<-c("CaseID","Sex","GA","Trimester","Tissue","CellType")

pd <- cbind(ID = mypheno$geo_accession, pd)
pd$ID <- as.character(pd$ID)
identical(pd$ID, colnames(B))

save(B, file = "B_GSE159526.rda")
save(pd, file = "pd_GSE159526.rda")

system("rm *.gz")
system("rm *.tar")
system("rm ../*.gz")
