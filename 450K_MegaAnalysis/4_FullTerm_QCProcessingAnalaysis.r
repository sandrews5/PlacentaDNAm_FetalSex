library(minfi)
load("RGset.term.nodups.rda")
colnames(RGset.term.nodups) <- seq(1,ncol(RGset.term.nodups))

################################
################################
####Processing and QC of data
################################
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")
noob <- preprocessNoob(RGset.term.nodups)

######
#Identify low intensity samples. 
rawobj <- preprocessRaw(RGset.term.nodups)
rawobj <- rawobj[match(rownames(noob),rownames(rawobj)),]
Meth <- log2(colMedians(getMeth(rawobj)))
Unmeth <- log2(colMedians(getUnmeth(rawobj)))
lowintes <- which(Meth<11 | Unmeth < 11) #66 samples with low overall intensity
rm(rawobj,Meth,Unmeth)

######
#Identify sample and probe detection p-value failures.
#Samples: Detection p-value > 0.01 at more than 1% of sites
#Probes: Detection p-value > 0.01 at more than 10% of samples 
detP <- detectionP(RGset.term.nodups)
detP <- detP[match(rownames(noob),rownames(detP)),]
failed <- detP > 0.01
samplesfailed <- which(colMeans(failed)>0.01) 
probesfailed <- which(rowMeans(failed)>0.1) 

######
#Remove probe failures and sample failures (either low intensity or detection p-value failures)
if(length(probesfailed) > 0){
	noob <- noob[-as.numeric(probesfailed),]
}
if(length(unique(c(lowintes,samplesfailed))) > 0){
	noob <- noob[,-unique(c(lowintes,samplesfailed))]
}
dim(noob) #Note the new dimensions of noob object.

######
#Attach positional information to noob object
noob <- mapToGenome(noob)

######
#Remove sex discrepant samples and attach predicted sex to phenotype data. 
predictedSex <- getSex(noob)
pd <- pData(noob)
#Change 'sexvar' below to whatever the sex variable is called in the pd object (can do colnames(pd) to check)
sexvar <- "Sex" ###CHANGE ME IF NECESSARY
#May need to change "pd[,sexvar] == X" depending on how variable is coded. (i.e "male" vs "M")
sexdiscrepant<-which((predictedSex$predictedSex=="M" & pd[,sexvar] == "F")|(predictedSex$predictedSex=="F" & pd[,sexvar] == "M")) 
predSex <- predictedSex$predictedSex
if(length(sexdiscrepant) > 0){
	noob <- noob[,-sexdiscrepant]
	predSex <- predSex[-sexdiscrepant]
}
dim(noob) #Note the new dimensions of noob object.

pd <- pData(noob)
pd$predSex <- predSex
pd$predSex <- factor(pd$predSex,levels=c("F","M"))
pData(noob) <- pd

######
#Remove ambiguously mapping probes. 
load("cross.probes.info.rda") #Chen et al. 2013 ambiguous mapping probe object. 
noob <- noob[-which(rownames(noob)%in%cross.probes.info$TargetID),]
dim(noob) #Note the new dimensions of noob object.

######
#Generate PCs on M values with and without sex chromosome data. 
chrnames <- as.character(seqnames(noob))
M <- getM(noob)
noob.noXY <- noob[-which(chrnames%in%c("chrX","chrY")),]
dim(noob.noXY)
M.noXY <- getM(noob.noXY)

PCs <- prcomp(t(M),center=T, scale.=F)
PCs.noXY <- prcomp(t(M.noXY),center=T, scale.=F)

#####
#Generate pairs() plots to get visual depiction of drivers of variability in data. 

pd <- pData(noob)
mycols<-c(2,5)
mynames<-c("Study","PredictedSex")
colorlist1 <- c('black','chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3',
             'mediumorchid2', 'turquoise3', 'wheat4', 'slategray2')
colorlist2<-c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
	"#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
	"#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

pdf("TotalPCA.pdf")
PropEx <- round((PCs$sdev^2/sum(PCs$sdev^2))*100,digits=2)
head(PropEx)
barplot(PropEx[1:8], las=2, xlab='', ylab='% Variance Explained')
PropEx <- round((PCs.noXY$sdev^2/sum(PCs.noXY$sdev^2))*100,digits=2)
head(PropEx)
barplot(PropEx[1:8], las=2, xlab='', ylab='% Variance Explained')
for (i in 1:length(mycols)){
#for (i in 1:5){
	print(i)
	myvar<-pd[,mycols[i]]
	if(class(myvar)=="factor"){levels(myvar)<-c(levels(myvar),"Missing")}
	myvar[which(is.na(myvar))]<-"Missing"
	myvar<-factor(myvar,levels=unique(myvar))
	colorlist<-colorlist1
	if(length(levels(myvar))>12){colorlist<-colorlist2}
	mycolors<-unlist(lapply(myvar,function(x){return(colorlist[which(levels(myvar)==x)])}))
	
	#WITH SEX CHR
	plot.new()
	par(mar=c(0, 0, 0, 0))
	legend("center",levels(myvar), 
	title=mynames[i],fill=colorlist[1:length(levels(myvar))])
	pairs(PCs$x[,1:3],labels=paste0("PC",1:3),col=mycolors,pch=20,cex=2)

	#NO SEX CHR
	pairs(PCs.noXY$x[,1:3],labels=paste0("PC",1:3),col=mycolors,pch=20,cex=2)
}
dev.off()

#Clean up workspace a bit. 
rm(M,M.noXY,PCs,PCs.noXY,noob)

######
#Normalization
norm <- preprocessQuantile(RGset)
norm <- norm[match(rownames(noob.noXY),rownames(norm)),match(colnames(noob.noXY),colnames(norm))]

#PC plots post normalization. 
M <- getM(norm)
PCs <- prcomp(t(M),center=T, scale.=F)
pd <- pData(norm)
pd$predictedSex <- factor(pd$predictedSex, levels = c("F","M"))

mycols<-c(2,5) #CHANGE ME
mynames<-c("Study","PredictedSex")
colorlist1 <- c('black','chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3',
             'mediumorchid2', 'turquoise3', 'wheat4', 'slategray2')
colorlist2<-c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
	"#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
	"#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

pdf("TotalPCA_postnorm.pdf")
PropEx <- round((PCs$sdev^2/sum(PCs$sdev^2))*100,digits=2)
head(PropEx)
barplot(PropEx[1:8], las=2, xlab='', ylab='% Variance Explained')
for (i in 1:length(mycols)){
#for (i in 1:5){
	print(i)
	myvar<-pd[,mycols[i]]
	if(class(myvar)=="factor"){levels(myvar)<-c(levels(myvar),"Missing")}
	myvar[which(is.na(myvar))]<-"Missing"
	myvar<-factor(myvar,levels=unique(myvar))
	colorlist<-colorlist1
	if(length(levels(myvar))>12){colorlist<-colorlist2}
	mycolors<-unlist(lapply(myvar,function(x){return(colorlist[which(levels(myvar)==x)])}))
	
	#WITH SEX CHR
	plot.new()
	par(mar=c(0, 0, 0, 0))
	legend("center",levels(myvar), 
	title=mynames[i],fill=colorlist[1:length(levels(myvar))])
	pairs(PCs$x[,1:3],labels=paste0("PC",1:3),col=mycolors,pch=20,cex=2)
}
dev.off()

######
#Batch effect correction and visualization
mod = model.matrix(~predictedSex, data = pd) 
library(sva)
svaout <- sva(M, mod, method="irw") 

regressme <- data.frame(svaout$sv)
colnames(regressme) <- paste0("SV",1:svaout$n.sv)
form1 <- as.formula(paste0("~",paste0(colnames(regressme),collapse="+")))
mymod <- model.matrix(form1, data = regressme)
library(limma)
fit <- lmFit(M,mymod)
resids <- residuals(fit,M)

PCs <- prcomp(t(resids),center=T, scale.=F)

mycols<-c(2,5) #CHANGE ME
mynames<-c("Study","PredictedSex")
colorlist1 <- c('black','chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3',
             'mediumorchid2', 'turquoise3', 'wheat4', 'slategray2')
colorlist2<-c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", 
	"#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
	"#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

pdf("TotalPCA_postnorm_postSVA.pdf")
PropEx <- round((PCs$sdev^2/sum(PCs$sdev^2))*100,digits=2)
head(PropEx)
barplot(PropEx[1:8], las=2, xlab='', ylab='% Variance Explained')
for (i in 1:length(mycols)){
#for (i in 1:5){
	print(i)
	myvar<-pd[,mycols[i]]
	if(class(myvar)=="factor"){levels(myvar)<-c(levels(myvar),"Missing")}
	myvar[which(is.na(myvar))]<-"Missing"
	myvar<-factor(myvar,levels=unique(myvar))
	colorlist<-colorlist1
	if(length(levels(myvar))>12){colorlist<-colorlist2}
	mycolors<-unlist(lapply(myvar,function(x){return(colorlist[which(levels(myvar)==x)])}))

	plot.new()
	par(mar=c(0, 0, 0, 0))
	legend("center",levels(myvar), 
	title=mynames[i],fill=colorlist[1:length(levels(myvar))])
	pairs(PCs$x[,1:3],labels=paste0("PC",1:3),col=mycolors,pch=20,cex=2)
}
dev.off()

################################
################################
####Single - site association analysis
################################

pheno <- "predictedSex"

#Use limma to run single site analysis. Annotate results with CHR and position. 
regressme <- data.frame(pd[,pheno],svaout$sv)
colnames(regressme) <- c(pheno,paste0("SV",1:svaout$n.sv))
form1 <- as.formula(paste0("~",paste0(colnames(regressme),collapse="+")))
mymod <- model.matrix(form1, data = regressme)
fit <- lmFit(M,mymod)
fit <- eBayes(fit)
singlesite <- topTable(fit,coef=2,number=nrow(M)) #singel site results

chrnames <- as.character(seqnames(norm))
pos <- as.numeric(start(norm))
singlesite$CHR <- chrnames[match(rownames(singlesite),rownames(norm))]
singlesite$Pos <- pos[match(rownames(singlesite),rownames(norm))]
save(singlesite, file = "singlesite.rda")

#SAVE OBJECTS
save(M, file = "M.rda")
save(B, file = "B.rda")
save(pd, file = "pd.rda")
save(regressme, file = "regressme.rda")
save(chrnames, file = "chrnames.rda")
save(pos, file = "pos.rda")

################################
################################
####SRegion-based analysis
################################

bumps <- bumphunter (M, design=mymod, coef=2,pos=pos,chr=chrnames, 
	pickCutoff=TRUE,nullMethod="bootstrap",B=1000,
	verbose=TRUE, smoothFunction=loessByCluster)

save(bumps, file = "bumps.rda")


