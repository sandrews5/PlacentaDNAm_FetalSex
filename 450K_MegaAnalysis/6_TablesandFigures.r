library(karyoploteR)
library(GenomicRanges)
library(ggplot2)
library(missMethyl)
library(topGO)
detach("package:missMethyl",unload = TRUE)
detach("package:IlluminaHumanMethylationEPICanno.ilm10b4.hg19",unload = TRUE)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Generate all the figures and tables for the analysis
load("singlesite.rda")
load("pd.rda")
load("B.rda")
B <- B[match(rownames(singlesite),rownames(B)),]
singlesite$MD <- rowMeans(B[,which(pd$predictedSex == "M")])-rowMeans(B[,which(pd$predictedSex == "F")])

studies <- c("GSE108567","GSE71678","GSE75248")
studyspec <- lapply(1:length(studies), function(x){
	pd.tmp <- pd[which(pd$Study == studies[x]),]
	B.tmp <- B[,which(pd$Study == studies[x])]

	out <- rowMeans(B.tmp[,which(pd.tmp$predictedSex == "M")]) - rowMeans(B.tmp[,which(pd.tmp$predictedSex == "F")])
	ifelse(sign(out) == sign(singlesite$MD),"Y","N")
})
studyspec <- do.call("cbind",studyspec)

singlesite$StudySpecific <- apply(studyspec,1,paste,collapse = "")

#############
#############
#Supplementary Data 1

suppdata1 <- data.frame(ID = rownames(singlesite))
suppdata1<- cbind(suppdata1, singlesite[c("CHR","Pos")])
suppdata1 <- cbind(suppdata1, singlesite[,!colnames(singlesite)%in%c("CHR","Pos")])
write.csv(suppdata1, "SupplementaryData1.csv", row.names = F, quote = F)

hits <- singlesite[singlesite$P.Value < 1E-8,]

#############
#############
#Supplementary Figure 5

data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
islandstat <- annotation.table$Relation_to_Island[match(rownames(singlesite),rownames(annotation.table))]
regfeature <- annotation.table$Regulatory_Feature_Group[match(rownames(singlesite),rownames(annotation.table))]

for(i in 1:length(unique(islandstat))){
	png(paste0("~/Documents/Misc/450k_FetalSex_Placenta/SciReports_Resubmission/Density_",unique(islandstat)[i],".png"),
		width = 10, height = 10, units = "in", res = 480)
	densityPlot(B[which(islandstat == unique(islandstat)[i]),],sampGroups = pd$predictedSex)
	dev.off()
}

myregs <- c("Gene_Associated","Gene_Associated_Cell_type_specific","Promoter_Associated","Promoter_Associated_Cell_type_specific")
for(i in 1:length(unique(myregs))){
	png(paste0("~/Documents/Misc/450k_FetalSex_Placenta/SciReports_Resubmission/Density_",unique(myregs)[i],".png"),
		width = 15, height = 10, units = "in", res = 480)
	densityPlot(B[which(regfeature == unique(myregs)[i]),],sampGroups = pd$predictedSex)
	dev.off()
}

island.plot <- data.frame(Group = rep("All\nProbes",length(unique(islandstat))),
	Count = table(islandstat))
names(island.plot) <- c("Group","IslandProx","Count")

island.plot.hits <- data.frame(Group = rep("Sig\nProbes",length(unique(islandstat))),
	Count = table(islandstat[which(singlesite$P.Value < 1E-8)]))
names(island.plot.hits) <- c("Group","IslandProx","Count")

island.plot <- rbind(island.plot,island.plot.hits)
island.plot$IslandProx <- factor(island.plot$IslandProx, 
	levels = c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","OpenSea"))

png("SuppFigure5C.png",
	width = 6, height = 4, units = "in", res = 480)
ggplot(island.plot, aes(fill=IslandProx, y=Count, x=Group)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() +
   guides(fill = guide_legend(title = "Proximity to\nIsland")) +
    ylab("Proportion") + 
    coord_flip()
dev.off()

chisq.test(rbind(table(islandstat),table(islandstat[which(singlesite$P.Value < 1E-8)])))

regfeature.plot <- data.frame(Group = rep("All\nProbes",length(unique(regfeature))),
	Count = table(regfeature))
colnames(regfeature.plot) <- c("Group","RegFeature","Count")

regfeature.plot.hits <- data.frame(Group = rep("Sig\nProbes",length(unique(regfeature))),
	Count = table(regfeature[which(singlesite$P.Value < 1E-8)]))
colnames(regfeature.plot.hits) <- c("Group","RegFeature","Count")

regfeature.plot <- rbind(regfeature.plot,regfeature.plot.hits)
regfeature.plot$RegFeature <- as.character(regfeature.plot$RegFeature)
regfeature.plot$RegFeature[which(regfeature.plot$RegFeature == "")] <- "None"


png("SuppFigure5D.png",
	width = 6, height = 4, units = "in", res = 480)
ggplot(regfeature.plot, aes(fill=RegFeature, y=Count, x=Group)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() +
   guides(fill = guide_legend(title = "Regulatory\nFeature")) +
    ylab("Proportion") + 
    coord_flip()
dev.off()

chisq.test(rbind(table(regfeature),table(regfeature[which(singlesite$P.Value < 1E-8)])))


#############
#############
#Figure 2

toplot <- data.frame(seqnames = as.numeric(unlist(lapply(singlesite$CHR,function(x){strsplit(x,"chr")[[1]][2]}))),
	start = singlesite$Pos, end = singlesite$Pos, pval = singlesite$P.Value)
toplot$seqnames <- paste0("chr", toplot$seqnames)
ewas <- makeGRangesFromDataFrame(toplot, keep.extra.columns = TRUE)
ewas <- sort(ewas)

toplot.malehyper <- toplot[which(singlesite$MD > 0),]
toplot.femalehyper <- toplot[which(singlesite$MD < 0),]
ewas.male <- makeGRangesFromDataFrame(toplot.malehyper, keep.extra.columns = TRUE)
ewas.male <- sort(ewas.male)
ewas.female <- makeGRangesFromDataFrame(toplot.femalehyper, keep.extra.columns = TRUE)
ewas.female <- sort(ewas.female)
 
png("Figure2A", width = 10, height = 4, units = "in", res = 480)
  kp <- plotKaryotype(plot.type=4, chromosomes = "autosomal")
  ticks <- seq(0,200,by = 50)
  kpAxis(kp, ymin=0, ymax=200, r0 = 0.5, tick.pos = ticks)
  kp <- kpPlotManhattan(kp, data=ewas.male, points.cex = 0.5, r0=0.5, r1=1,ymax = 200, 
  	suggestive.col = "white", genomewide.col = "darkorange", genomewideline = 8, genomewide.lwd = 1.5)
  kpAxis(kp, ymin=0, ymax=200, r0=0.5, r1=0, tick.pos = ticks)
  kp <- kpPlotManhattan(kp, data=ewas.female, points.cex = 0.5, r0=0.5, r1=0, ymax = 200, 
  	suggestive.col = "white", genomewide.col = "darkorange", genomewideline = 8, genomewide.lwd = 1.5)
dev.off()

hits <- hits[order(abs(hits$MD), decreasing = T),]
cglist <- c("cg01382982","cg22905511")
fignames <- c("Figure2B","Figure2C")
mytrunc <- B[match(cglist, rownames(B)),]
singlesite[match(cglist,rownames(singlesite)),]

for (i in 1:length(cglist)){
	toplot <- data.frame(Beta = mytrunc[i,], Sex = pd$predictedSex)

	png(paste0(fignames[i],".png"),
		width = 3.33, height = 3.33, units = 'in', res = 480)

	p1 <- ggplot(toplot, aes(x = Sex, y = Beta, color = Sex)) + 
		scale_color_manual(guide = FALSE, values = c("goldenrod","dodgerblue")) + 
		geom_boxplot(outlier.color = NA) + ylab("Percent DNA Methylation") + xlab("Fetal Sex") +
		scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),
			labels = seq(0,100,by = 25)) +
		geom_jitter(position=position_jitter(0.2)) + theme_bw()  

	print(p1)

	dev.off()
}

#Figure 2D: Histogram of mean differences in hits
hits$sign <- sign(hits$MD)
hits$sign <- ifelse(hits$sign == 1, "M","F")
hits$sign <- factor(hits$sign, levels = c("F","M"))
toplot <- with(hits, data.frame(MD = abs(MD), Sign = sign))
png("Figure2D.png",
	width = 3.33, height = 3.33, units = 'in', res = 480)

p1 <- ggplot(toplot, aes(x = MD, fill = Sign)) + geom_density(alpha=0.4) +
	scale_fill_manual(name = "Hypermethylated in:",values = c("goldenrod","dodgerblue")) + 
	theme_bw() + xlab("Absolute value of Mean Difference (%)") +
	scale_x_continuous(breaks = seq(0,0.2,by = 0.05), limits = c(0,0.23),
		labels = seq(0,20,by = 5)) +
	theme(legend.justification = c(1,1), legend.position = c(0.95, 0.95))
print(p1)
dev.off()
with(toplot, ks.test(MD[Sign == "M"], MD[Sign == "F"]))


#############
#############
#Supplementary Figure 6

hits <- hits[order(abs(hits$MD), decreasing = T),]
cglist <- c("cg01382982","cg22905511")
fignames <- c("FigureS6B","FigureS6C")
mytrunc <- B[match(cglist, rownames(B)),]
singlesite[match(cglist,rownames(singlesite)),]

mytrunc <- mytrunc[,pd$Study%in%studies]
pd.trunc <- pd[pd$Study%in%studies,]

for (i in 1:length(cglist)){
	toplot <- data.frame(Beta = mytrunc[i,], Sex = pd.trunc$predictedSex, Study = pd.trunc$Study)

	print(unlist(lapply(studies, function(x){
		blah <- toplot[which(toplot$Study == x),]
		mean(blah$Beta[blah$Sex == "M"]) - mean(blah$Beta[blah$Sex == "F"])
	})))

	png(paste0(fignames[i],".png"),
		width = 6, height = 3, units = 'in', res = 480)

	p1 <- ggplot(toplot, aes(x = Sex, y = Beta, color = Sex)) + 
		scale_color_manual(guide = FALSE, values = c("goldenrod","dodgerblue")) + 
		geom_boxplot(outlier.color = NA) + ylab("Percent DNA Methylation") + xlab("Fetal Sex") +
		scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),
			labels = seq(0,100,by = 25)) +
		geom_jitter(position=position_jitter(0.2)) + theme_bw() +
		facet_wrap(~Study)

	print(p1)

	dev.off()
}

volcano <- data.frame(MD = singlesite$MD*100, p = -log10(singlesite$P.Value),
	Significant = ifelse(singlesite$P.Value < 1E-8,1,0),Name = rownames(singlesite))
volcano$Significant <- factor(volcano$Significant,levels = c(0,1))
volcano$Name[!volcano$Name %in% cglist] <- ""

png("FigureS6A.png",
	width = 6, height = 6, units = 'in', res = 480)
p1 <- ggplot(volcano, aes(x = MD, y = p, color = Significant)) + 
geom_point(shape = 20, size = 1) + theme_bw() + 
xlab("Mean Difference (% in Males - % in Females)") +
ylab("-log10 p-value") + xlim(-25,25) + 
scale_color_manual(name = "Genome-wide\nSignificant", labels = c("No","Yes"),
	values = c("gray","purple"))
print(p1)
dev.off()

#############
#############
#Table 1

cutoff <- 1E-8
sigcpgs <- rownames(singlesite)[which(singlesite$P.Value < cutoff)]
allcpgs <- rownames(singlesite)
gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, plot.bias=FALSE, prior.prob=TRUE)
table(gst$FDR<0.05)
res <- topGSA(gst,number=Inf)
toppaths <- res[res$FDR < 0.05,]
toppaths$GO_ID <- rownames(toppaths)
toppaths <- toppaths[,c("GO_ID","ONTOLOGY","TERM","N","DE","P.DE","FDR")]
toppaths$FDR <- round(toppaths$FDR)
write.csv(toppaths, file = "Table1.csv", row.names = F, quote = F)

#############
#############
#Supplementary Data 2

res$GO_ID <- rownames(res)
res<- res[,c("GO_ID","ONTOLOGY","TERM","N","DE","P.DE","FDR")]
res$FDR <- round(res$FDR)
write.csv(res, file = "SupplementaryData2.csv", row.names = F, quote = F)
system("SupplementaryData2.csv")

#############
#############
#Table 2

#Concordance with other studies
#Martin et al. 2017 - Placenta
#2561 probes that were replicated in the analysis that adjusted for gestational age
martin <- read.csv("../LiteratureResults/Martin2017Tophits.csv", header = F, stringsAsFactors = F)
martin$chr <- Locations$chr[match(martin$V1, rownames(Locations))]
martin$Inc <- martin$V1 %in% rownames(hits) #1 CpG site, lol
martin$IncTotal <- martin$V1 %in% rownames(singlesite)
with(martin[!martin$chr == "chrX",], table(Inc))
with(martin[!martin$chr == "chrX",], table(IncTotal))
# FALSE  TRUE 
#    14     1 

#Yousefi et al. 2015 - cord blood
yousefi <- read.csv("../LiteratureResults/Yousefi.csv", header = T, stringsAsFactors = F)
table(yousefi$ProbeID %in% rownames(hits))
table(yousefi$ProbeID %in% rownames(singlesite))
# FALSE  TRUE 
#  2713   318
# TRUE 
# 3031 


#Suderman et al
suderman <- read_excel("../LiteratureResults/Suderman.xlsx", sheet = 2)
suderman$Placenta <- suderman$targetid%in%rownames(hits)
suderman$PlacentaAll <- suderman$targetid%in%rownames(singlesite)
with(suderman,table(different.0))
with(suderman,table(different.7))
with(suderman,table(different.17))
with(suderman,table(different.0, Placenta))
with(suderman,table(different.7, Placenta))
with(suderman,table(different.17, Placenta))
with(suderman,table(different.0, PlacentaAll))
with(suderman,table(different.7, PlacentaAll))
with(suderman,table(different.17, PlacentaAll))

#Price et al. 
price <- read_excel("../LiteratureResults/Price.xlsx", sheet = 1, skip = 1)
price$Placenta <- price$ID_REF%in%rownames(hits)
price$PlacentaAll <- price$ID_REF%in%rownames(singlesite)
with(price,table(XY_Hits,Placenta))
with(price,table(XY_Hits,PlacentaAll))

#Xu et al.
xu <- read_excel("../LiteratureResults/Xu.xlsx", sheet = 1, skip = 1)
xu$autosomal <- !xu$Chromosome %in% c("X","Y")
xu$Placenta <- xu$CpGs%in%rownames(hits)
xu$PlacentaAll <- xu$CpGs%in%rownames(singlesite)
with(xu,table(autosomal,Placenta))

#Spiers et al.
spiers <- read.csv("../LiteratureResults/Spiers.csv", header = T, stringsAsFactors = F)
spiers$Placenta <- spiers$Probe%in%rownames(hits)
spiers$PlacentaAll <- spiers$Probe%in%rownames(singlesite)

#Inoshita et al.
inoshita <- read_excel("../LiteratureResults/Inoshita.xlsx", sheet = 1, skip = 1)
inoshita <- data.frame(inoshita)
inoshita$Placenta <- inoshita[,2]%in%rownames(hits)
inoshita$PlacentaAll <- inoshita[,2]%in%rownames(singlesite)

#Singmann
singmann <- read_excel("../LiteratureResults/Singmann.xlsx", sheet = 1, skip = 2)
singmann.sig <- singmann[!is.na(singmann$Pearson_p.value.ALSPAC),]
singmann.sig$Placenta <- singmann.sig$CpG%in%rownames(hits)
singmann.sig$PlacentaAll <- singmann.sig$CpG%in%rownames(singlesite)

#Xia
xia <- read_excel("../LiteratureResults/Xia.xlsx", sheet = 1, skip = 2)
xia$autosomal <- !xia$CHR %in% c("X","Y")
xia <- data.frame(xia)
xia$Placenta <- xia[,1]%in%rownames(hits)
xia$PlacentaAll <- xia[,1]%in%rownames(singlesite)
with(xia,table(autosomal,Placenta))
with(xia,table(autosomal,PlacentaAll))

#Hall
hall.p1 <- read_excel("../LiteratureResults/Hall_p1.xlsx", sheet = 1, skip = 1)
hall.p2 <- read_excel("../LiteratureResults/Hall_p2.xlsx", sheet = 1, skip = 1)
hall <- data.frame(rbind(hall.p1,hall.p2))
hall$Placenta <- hall[,1]%in%rownames(hits)
hall$PlacentaAll <- hall[,1]%in%rownames(singlesite)

#############
#############
#Supplementary Figure 7

celltype <- read_excel("12864_2020_7186_MOESM4_ESM.xlsx", sheet = 1) #From Yuan et al.

#Does the distribution of MD differ in these groups?
hits$celltype <- as.numeric(rownames(hits)%in%celltype$cpg)
hits$celltype <- factor(hits$celltype, levels = c(0,1))
toplot <- with(hits, data.frame(MD = abs(MD), CellType = celltype))

pdf("SupplementaryFigure7.pdf",
	width = 3.33, height = 3.33)
p1 <- ggplot(toplot, aes(x = MD, fill = CellType)) + geom_density(alpha=0.4) +
	scale_fill_manual(name = "Cell Type Distinguishing",
		values = c("#f7931e","grey50"), labels = c("No","Yes")) + 
	theme_bw() + xlab("Absolute value of Mean Difference (%)") +
	scale_x_continuous(breaks = seq(0,0.2,by = 0.05), limits = c(0,0.23),
		labels = seq(0,20,by = 5)) +
	theme(legend.justification = c(1,1), legend.position = c(0.95, 0.95))
print(p1)
dev.off()

#############
#############
#Supplementary Data 3 and 4

allcpgs <- rownames(singlesite)
gst <- gometh(sig.cpg = rownames(hits)[hits$celltype==0], all.cpg = allcpgs, plot.bias=FALSE, prior.prob=TRUE)
table(gst$FDR<0.05)
res <- topGSA(gst,number=Inf)
toppaths <- res#[res$FDR < 0.05,]
toppaths$GO_ID <- rownames(toppaths)
toppaths <- toppaths[,c("GO_ID","ONTOLOGY","TERM","N","DE","P.DE","FDR")]
write.csv(toppaths, file = "SupplementaryData3.csv", row.names = F, quote = F)
system("gzip SupplementaryData3.csv")

allcpgs <- rownames(singlesite)
gst <- gometh(sig.cpg = rownames(hits)[hits$celltype==1], all.cpg = allcpgs, plot.bias=FALSE, prior.prob=TRUE)
table(gst$FDR<0.05)
res <- topGSA(gst,number=Inf)
toppaths <- res#[res$FDR < 0.05,]
toppaths$GO_ID <- rownames(toppaths)
toppaths <- toppaths[,c("GO_ID","ONTOLOGY","TERM","N","DE","P.DE","FDR")]
write.csv(toppaths, file = "SupplementaryData4.csv", row.names = F, quote = F)
system("gzip SupplementaryData4.csv")


#############
#############
#Figure 3

fullterm <- singlesite
pd.full <- pd
B.full <- B

term3 <- get(load("singlesite_Trim3.rda"))
B.term3 <- get(load("B_Trim3.rda"))
B.term3 <- B.term3[match(rownames(term3),rownames(B.term3)),]
pd.term3 <- get(load("pd_Trim3.rda"))
term3$MD <- rowMeans(B.term3[,which(pd.term3$predictedSex == "M")])-rowMeans(B.term3[,which(pd.term3$predictedSex == "F")])
#write.csv(term3[,c("CHR","Pos","P.Value","MD")], 
#	file = "../Figures&Tables/SingleSiteResults_ThirdTrimester.csv")
#Trim 2
term2 <- get(load("singlesite_Trim2.rda"))
B.term2 <- get(load("B_Trim2.rda"))
B.term2 <- B.term2[match(rownames(term2),rownames(B.term2)),]
pd.term2 <- get(load("pd_Trim2.rda"))
term2$MD <- rowMeans(B.term2[,which(pd.term2$predictedSex == "M")])-rowMeans(B.term2[,which(pd.term2$predictedSex == "F")])
#write.csv(term2[,c("CHR","Pos","P.Value","MD")], 
#	file = "../Figures&Tables/SingleSiteResults_SecondTrimester.csv")
#Trim 1
term1 <- get(load("singlesite_Trim1.rda"))
B.term1 <- get(load("B_Trim1.rda"))
B.term1 <- B.term1[match(rownames(term1),rownames(B.term1)),]
pd.term1 <- get(load("pd_Trim1.rda"))
term1$MD <- rowMeans(B.term1[,which(pd.term1$predictedSex == "M")])-rowMeans(B.term1[,which(pd.term1$predictedSex == "F")])

ga.term <- predictAge(B, type = 'RPC')
ga.term3 <- predictAge(B.term3, type = 'RPC')
ga.term2 <- predictAge(B.term2, type = 'RPC')
ga.term1 <- predictAge(B.term1, type = 'RPC')

toplot <- data.frame(GA = c(ga.term,ga.term3,ga.term2,ga.term1),
	AssignedGroup = c(rep("Full Term",length(ga.term)),rep("Trimester 3",length(ga.term3)),
		rep("Trimester 2",length(ga.term2)),rep("Trimester 1",length(ga.term1))))
toplot$AssignedGroup <- factor(toplot$AssignedGroup, levels = c("Full Term","Trimester 3","Trimester 2","Trimester 1"))

png("Figure3A.png",width = 13.33, height = 3, units = 'in', res = 480)

p1 <- ggplot(toplot, aes(x = GA, fill = AssignedGroup)) + geom_density(alpha=0.4) +
	scale_fill_manual(name = "Assigned Gestational Period",
		values = c("coral1","cyan1","orange1","plum1")) + 
	theme_bw() + xlab("Predicted Gestational Age") +
	scale_x_continuous(breaks = seq(0,42,by = 5), limits = c(0,42),
		labels = seq(0,42,by = 5)) 
print(p1)
dev.off()

ssobjs <- list(term3, term2, term1)
termlabs <- c("SupplementalFigure8A","SupplementalFigure8B","SupplementalFigure8C")
supplabs <- c("SupplementaryData5",'SupplementaryData6',"SupplementaryData7")
maxvals <- c(16,20,12)

for(myterm in 1:length(ssobjs)){
	singlesite <- ssobjs[[myterm]]
	toplot <- data.frame(seqnames = singlesite$CHR,
		start = singlesite$Pos, end = singlesite$Pos, pval = singlesite$P.Value)
	
	toplot.malehyper <- toplot[which(singlesite$MD > 0),]
	toplot.femalehyper <- toplot[which(singlesite$MD < 0),]
	ewas.male <- makeGRangesFromDataFrame(toplot.malehyper, keep.extra.columns = TRUE)
	ewas.male <- sort(ewas.male)
	ewas.female <- makeGRangesFromDataFrame(toplot.femalehyper, keep.extra.columns = TRUE)
	ewas.female <- sort(ewas.female)

	png(paste0(termlabs[myterm],".png"), width = 10, height = 4, units = "in", res = 480)
		kp <- plotKaryotype(plot.type=4, chromosomes = "autosomal")
		  ticks <- seq(0,maxvals[myterm],by = 2)
		  kpAxis(kp, ymin=0, ymax=maxvals[myterm], r0 = 0.5, tick.pos = ticks)
		  kp <- kpPlotManhattan(kp, data=ewas.male, points.cex = 0.5, r0=0.5, r1=1,ymax = maxvals[myterm], 
		  	suggestive.col = "white", genomewide.col = "darkorange", genomewideline = 8, genomewide.lwd = 1.5)
		  kpAxis(kp, ymin=0, ymax=maxvals[myterm], r0=0.5, r1=0, tick.pos = ticks)
		  kp <- kpPlotManhattan(kp, data=ewas.female, points.cex = 0.5, r0=0.5, r1=0, ymax = maxvals[myterm], 
		  	suggestive.col = "white", genomewide.col = "darkorange", genomewideline = 8, genomewide.lwd = 1.5)
	dev.off()	

	supptable1<-data.frame(ID = rownames(singlesite))
	supptable1 <- cbind(supptable1, singlesite[c("CHR","Pos")])
	supptable1 <- cbind(supptable1, singlesite[,!colnames(singlesite)%in%c("CHR","Pos")])
	write.csv(supptable1, paste0(supplabs[myterm],".csv"), row.names = F, quote = F)
	system(paste0("gzip ",supplabs[myterm],".csv"))

}

commoncgs <- intersect(rownames(fullterm),rownames(term1))
commoncgs <- intersect(commoncgs,rownames(term2))
commoncgs <- intersect(commoncgs,rownames(term3))

fullterm.common <- fullterm[which(rownames(fullterm)%in%commoncgs),]
dim(fullterm.common) 
term3 <- term3[which(rownames(term3)%in%commoncgs),]
term3 <- term3[match(rownames(fullterm.common),rownames(term3)),]
term2 <- term2[which(rownames(term2)%in%commoncgs),]
term2 <- term2[match(rownames(fullterm.common),rownames(term2)),]
term1 <- term1[which(rownames(term1)%in%commoncgs),]
term1 <- term1[match(rownames(fullterm.common),rownames(term1)),]

fullterm.sig <- fullterm.common[which(fullterm.common$P.Value < 1E-8),]
term3.sig <- term3[which(fullterm.common$P.Value < 1E-8),]
term2.sig <- term2[which(fullterm.common$P.Value < 1E-8),]
term1.sig <- term1[which(fullterm.common$P.Value < 1E-8),]

cutoff <- 0.05
term3.sig$Pcheck <- with(term3.sig, P.Value < cutoff)
term3.sig$Direction <- sign(fullterm.sig$MD) == sign(term3.sig$MD)
term2.sig$Pcheck <- with(term2.sig, P.Value < cutoff)
term2.sig$Direction <- sign(fullterm.sig$MD) == sign(term2.sig$MD)
term1.sig$Pcheck <- with(term1.sig, P.Value < cutoff)
term1.sig$Direction <- sign(fullterm.sig$MD) == sign(term1.sig$MD)

term3.sig$Conserve <- with(term3.sig, Pcheck + Direction)
term2.sig$Conserve <- with(term2.sig, Pcheck + Direction)
term1.sig$Conserve <- with(term1.sig, Pcheck + Direction)

totalconserve <- term3.sig$Conserve + term2.sig$Conserve + term1.sig$Conserve

png("Figure3B", width = 10, height = 4, units = "in", res = 480)
	layout(matrix(c(1,2,3),1,3,byrow = TRUE))
	quadtable <- data.frame(Term = 100*fullterm.sig$MD, Third = 100*term3.sig$MD,
		Second = 100*term2.sig$MD,First = 100*term1.sig$MD, P = term1.sig$P.Value, totalconserve = totalconserve)
	quadtable <- quadtable[order(quadtable$P, quadtable$totalconserve, decreasing = TRUE),]
	colvec <- rep("black",nrow(quadtable))
	colvec[which(quadtable$P < 0.05)] <- "deepskyblue1"
	colvec[which(quadtable$totalconserve == 6)] <- "chartreuse"
	with(quadtable, plot(Term,First, pch = 20, main=paste0("Corr = ", round(cor(Term,First),3)), cex=1.2,
		xlab="Fetal sex % DNAm difference, Full Term",ylab="Fetal sex % DNAm difference, First Trimester", xlim=range(quadtable[,1:4]), 
		ylim=range(quadtable[,1:4]), cex.lab=1.2,cex.axis=1.2,col=colvec))
	abline(v=0,lty="dashed")
	abline(h=0,lty="dashed")
	legend("bottomright",legend=c("p \u2265 0.05, First Trimester","p < 0.05, First Trimester","Conserved Across Gestation"),fill=c("black","deepskyblue1","chartreuse"))
	quadtable <- data.frame(Term = 100*fullterm.sig$MD, Third = 100*term3.sig$MD,
		Second = 100*term2.sig$MD,First = 100*term1.sig$MD, P = term2.sig$P.Value, totalconserve = totalconserve)
	quadtable <- quadtable[order(quadtable$P, quadtable$totalconserve, decreasing = TRUE),]
	colvec <- rep("black",nrow(quadtable))
	colvec[which(quadtable$P < 0.05)] <- "deepskyblue1"
	colvec[which(quadtable$totalconserve == 6)] <- "chartreuse"
	with(quadtable, plot(Term,Second, pch = 20, main=paste0("Corr = ", round(cor(Term,Second),3)), cex=1.2,
		xlab="Fetal sex % DNAm difference, Full Term",ylab="Fetal sex % DNAm difference, Second Trimester", xlim=range(quadtable[,1:4]), 
		ylim=range(quadtable[,1:4]), cex.lab=1.2,cex.axis=1.2,col=colvec))
	abline(v=0,lty="dashed")
	abline(h=0,lty="dashed")
	legend("bottomright",legend=c("p \u2265 0.05, Second Trimester","p < 0.05, Second Trimester","Conserved Across Gestation"),fill=c("black","deepskyblue1","chartreuse"))
	quadtable <- data.frame(Term = 100*fullterm.sig$MD, Third = 100*term3.sig$MD,
		Second = 100*term2.sig$MD,First = 100*term1.sig$MD, P = term3.sig$P.Value, totalconserve = totalconserve)
	quadtable <- quadtable[order(quadtable$P, quadtable$totalconserve, decreasing = TRUE),]
	colvec <- rep("black",nrow(quadtable))
	colvec[which(quadtable$P < 0.05)] <- "deepskyblue1"
	colvec[which(quadtable$totalconserve == 6)] <- "chartreuse"
	with(quadtable, plot(Term,Third, pch = 20, main=paste0("Corr = ", round(cor(Term,Third),3)), cex=1.2,
		xlab="Fetal sex % DNAm difference, Full Term",ylab="Fetal sex % DNAm difference, Third Trimester", xlim=range(quadtable[,1:4]), 
		ylim=range(quadtable[,1:4]), cex.lab=1.2,cex.axis=1.2,col=colvec))
	abline(v=0,lty="dashed")
	abline(h=0,lty="dashed")
	legend("bottomright",legend=c("p \u2265 0.05, Third Trimester","p < 0.05, Third Trimester","Conserved Across Gestation"),fill=c("black","deepskyblue1","chartreuse"))
dev.off()

#############
##Figure 3C
#############

totalconserve <- term3.sig$Conserve + term2.sig$Conserve + term1.sig$Conserve
sum(totalconserve == 6)
totalconserve <- term3.sig$Conserve + term2.sig$Conserve
sum(totalconserve == 4)
totalconserve <- term3.sig$Conserve 
sum(totalconserve == 2)

> sum(totalconserve == 6)
[1] 194
> totalconserve <- term3.sig$Conserve + term2.sig$Conserve
> sum(totalconserve == 4)
[1] 698
> totalconserve <- term3.sig$Conserve 
> sum(totalconserve == 2)
[1] 1596

#Individual site across different time periods

mycg <- "cg17612569"

myss.fullterm <- data.frame(Beta = B.full[match(mycg,rownames(B.full)),], 
	Sex = pd.full$predictedSex)
myss.fullterm$Term <- "Full Term"
myss.term1 <- data.frame(Beta = B.term1[match(mycg,rownames(B.term1)),], 
	Sex = pd.term1$predictedSex)
myss.term1$Term <- "Trimester 1"
myss.term2 <- data.frame(Beta = B.term2[match(mycg,rownames(B.term2)),], 
	Sex = pd.term2$predictedSex)
myss.term2$Term <- "Trimester 2"
myss.term3 <- data.frame(Beta = B.term3[match(mycg,rownames(B.term3)),], 
	Sex = pd.term3$predictedSex)
myss.term3$Term <- "Trimester 3"

toplot <- rbind(myss.term1, myss.term2, myss.term3, myss.fullterm)
toplot$Sex <- factor(toplot$Sex, levels = c("F","M"))
toplot$Term <- factor(toplot$Term, levels = c("Trimester 1", "Trimester 2", "Trimester 3", "Full Term"))

png("Figure3D.png",width = 12.5, height = 3.33, units = 'in', res = 480)
p1 <- ggplot(toplot, aes(x = Term, y = Beta, color = Sex)) + 
	scale_color_manual(name = "Fetal Sex", values = c("goldenrod","dodgerblue")) + 
	geom_boxplot(outlier.color = NA) + ylab("Percent DNA Methylation") + xlab("Gestational Period") +
	scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),
		labels = seq(0,100,by = 25)) + theme_bw() +
	geom_point(shape = 20, position = position_jitterdodge(), size = 2)    

print(p1)
dev.off()

#############
#############
#Table 3

pickme <- which(totalconserve == 6)

table3 <- data.frame(ID = rownames(fullterm.sig)[pickme])
table3$FullTerm.MD <- with(fullterm.sig[pickme,],round(MD*100,1))
table3$FullTerm.pval <- with(fullterm.sig[pickme,],signif(P.Value,3))
table3$ThirdTrim.MD <- with(term3.sig[pickme,],round(MD*100,1))
table3$ThirdTrim.pval <- with(term3.sig[pickme,],signif(P.Value,3))
table3$SecondTrim.MD <- with(term2.sig[pickme,],round(MD*100,1))
table3$SecondTrim.pval <- with(term2.sig[pickme,],signif(P.Value,3))
table3$FirstTrim.MD <- with(term1.sig[pickme,],round(MD*100,1))
table3$FirstTrim.pval <- with(term1.sig[pickme,],signif(P.Value,3))

table3 <- table3[order(abs(table3$FullTerm.MD), decreasing = T),]
write.csv(table3, file = "Table3.csv", row.names = F, quote = F)

#############
#############
#Supplementary Figure 9

load("bumps.rda")
load("pd.rda")
load("B.rda")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
chrnames <- Locations$chr[match(rownames(B),rownames(Locations))]
pos <- Locations$pos[match(rownames(B),rownames(Locations))]

dmrs <- bumps$table[order(bumps$table$fwerArea),]
dmrs$MD <- unlist(lapply(1:nrow(dmrs),function(x){
	if(dmrs$L[x] > 1){
		mean(colMeans(B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "M")]))-
		mean(colMeans(B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "F")]))
	} else {
		mean((B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "M")]))-
		mean((B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "F")]))
	}
  }))
dmrs$width <- with(dmrs, end - start)

dmrs <- dmrs[order(dmrs$fwerArea, decreasing = TRUE),]
dmrs$col <- factor(ifelse(dmrs$fwerArea < 0.1, 1, 0), levels = c(0,1))
dmrs.full <- dmrs
dmrs.full <- dmrs.full[order(dmrs.full$fwerArea, decreasing = FALSE),]


png("SupplementaryFigure9", ,width = 5, height = 3.75, units = 'in', res = 480)
p1 <- ggplot(dmrs.full,aes(x = MD, y = width, col = col))+ geom_point() + 
	scale_color_manual(name = "",values = c("black","deepskyblue1"), labels = c("fwerArea \u2265 0.1","fwerArea < 0.1")) + 
	theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), 
		legend.position=c(0.025,0.95), legend.justification='left',legend.direction='vertical',
    	panel.background = element_rect(fill = "white")) + 
	xlim(-0.15,0.20)+
	ylab("DMR width (bp)") + xlab("Mean DNAm in Males - Mean DNAm in Females")
print(p1)
dev.off()

#############
#Supplementary Data 8

supptab8 <- dmrs.full[,c("chr","start","end","p.value","fwer","p.valueArea","fwerArea","MD")]
write.csv(supptab8, file = "SupplementaryTable8.csv", row.names = F, quote = F)
system("gzip SupplementaryTable8.csv")


#############
#Supplementary Data 9-11

term.gr <- makeGRangesFromDataFrame(dmrs)
basedir <- "EarlyTermAnalysisObjects/"
fronts <- paste0("Trimester",1:3,"/") 
backs <- paste0("Trim",1:3,".rda")

bumps <- lapply(1:3,function(x){
	load(paste0(basedir,fronts[x],"bumps_",backs[x]))
	return(bumps$table)
})

myterm <- bumps[[3]]
myterm <- myterm[order(myterm$fwerArea, decreasing = FALSE),]
supptab <- myterm[,c("chr","start","end","p.value","fwer","p.valueArea","fwerArea")]
write.csv(supptab, file = "SupplementaryData9.csv", row.names = F, quote = F)
system("SupplementaryData9.csv")

myterm <- bumps[[2]]
myterm <- myterm[order(myterm$fwerArea, decreasing = FALSE),]
supptab <- myterm[,c("chr","start","end","p.value","fwer","p.valueArea","fwerArea")]
write.csv(supptab, file = "SupplementaryData10.csv", row.names = F, quote = F)
system("gzip SupplementaryData10.csv")

myterm <- bumps[[1]]
myterm <- myterm[order(myterm$fwerArea, decreasing = FALSE),]
supptab <- myterm[,c("chr","start","end","p.value","fwer","p.valueArea","fwerArea")]
write.csv(supptab, file = "SupplementaryData11.csv", row.names = F, quote = F)
system("gzip SupplementaryData11.csv")

#############
#Table 4

Betas <- lapply(1:3,function(x){
	load(paste0(basedir,fronts[x],"B_",backs[x]))
	return(B)
})
phenos <- lapply(1:3,function(x){
	load(paste0(basedir,fronts[x],"pd_",backs[x]))
	return(pd)
})

earlybumps.gr <- lapply(bumps,function(x){
		makeGRangesFromDataFrame(x)
})
earlymeth.gr <- lapply(Betas,function(x){
	tmp <- data.frame(chr = Locations$chr[match(rownames(x),rownames(Locations))],
				start = Locations$pos[match(rownames(x),rownames(Locations))],
				end = Locations$pos[match(rownames(x),rownames(Locations))])
	makeGRangesFromDataFrame(tmp)
	})

#Overlap of actual bumphunter region
bumpov <- lapply(earlybumps.gr, function(x){
	findOverlaps(term.gr, x)
	})

#Overlap of probes in region
matchMDs <- lapply(1:3, function(x){
	ov <- findOverlaps(term.gr, earlymeth.gr[[x]])
	 pd <- phenos[[x]]
	 B <- Betas[[x]]
	MD <- unlist(lapply(1:length(term.gr),function(y){
		indxs <- sort(subjectHits(ov)[which(queryHits(ov)==y)])
		 mean(colMeans(B[indxs,which(pd$predictedSex == "M")]))-
  		mean(colMeans(B[indxs,which(pd$predictedSex == "F")]))
	}))
	return(MD)
})

#Gene info
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb)
library(annotate)
library(org.Hs.eg.db)
genenames<-getSYMBOL(genes(txdb)$gene_id,data='org.Hs.eg')
genes.df<-data.frame(genes)
genes.df$Name <- genenames

tab4 <- dmrs[,c("chr","start","end","p.value","fwer","p.valueArea","fwerArea","MD")]
tab4$p.value <- signif(tab4$p.value,3)
tab4$p.valueArea <- signif(tab4$p.valueArea,3)
tab4$MD <- signif(tab4$MD,3)
tab4$Gene <- genes.df$Name[nearest(term.gr,genes)]
tab4$MD_Trim1 <- signif(matchMDs[[1]],3)
tab4$MD_Trim2 <- signif(matchMDs[[2]],3)
tab4$MD_Trim3 <- signif(matchMDs[[3]],3)
write.csv(tab4, file = "Table4.csv", row.names = F, quote = F)

#############
#Figure 4 and Supplementary Figures 10-17

dmrs <- dmrs[order(dmrs$fwerArea, decreasing = FALSE),]
dmrs <- dmrs[which(dmrs$fwerArea < 0.1),]

library(ggplot2)
library(gridExtra)
library(karyoploteR)
labs <- c("Full Term","1st Trimester","2nd Trimester","3rd Trimester")
for (mydmr in 2:nrow(dmrs)){
	periods <- list()
	for (i in 1:4){
		periods[[i]] <- methylationplot2(dmrs[mydmr,],Betas[[i]],phenos[[i]]$predictedSex, ylab = paste0("% DNAm, ",labs[i]),
			colchoices = c("goldenrod","dodgerblue"), legendLabel = "Fetal Sex",
			leftBufferSize = 100, rightBufferSize = 100)		
	}

	png(paste0("DMR",mydmr,".png"), width = 6, height = 8, units = "in", res = 480)
		grid.arrange(periods[[1]],periods[[4]],periods[[3]],periods[[2]],nrow = 4, ncol = 1)
	dev.off()

	res.gr <- dmrs[mydmr,]
	res.gr$start <- res.gr$start - 100
	res.gr$end <- res.gr$end + 100
	res.gr <- makeGRangesFromDataFrame(res.gr)

	png(paste0("GeneTrack_",mydmr,".png"), ,width = 6, height = 2, units = 'in', res = 480)
	  kp <- plotKaryotype(plot.type=4, zoom = range(res.gr))

	   genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
	  	genes.data <- addGeneNames(genes.data)
	     genes.data <- mergeTranscripts(genes.data)
	     kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8,
	                 gene.name.position = "left")
 	dev.off()

}

#############
#Figure 5 and Supplementary Figures 18-21

load("GSE159526/B_GSE159526.rda")
load("GSE159526/pd_GSE159526.rda")
Locations <- Locations[match(rownames(B),rownames(Locations)),]
Locations <- with(Locations, data.frame(chr = chr, start = pos, end = pos, ID = rownames(Locations)))
B <- B[!is.na(Locations$chr),]
Locations <- Locations[!is.na(Locations$chr),]
Locations <- makeGRangesFromDataFrame(Locations,keep.extra.columns = TRUE)
B <- B[order(Locations),]
Locations <- Locations[order(Locations),]

B.celltype <- B
Locations.celltype <- Locations
pd.celltype <- pd

load("bumps.rda")
load("pd.rda")
load("B.rda")
detach("package:missMethyl",unload = TRUE)
detach("package:IlluminaHumanMethylationEPICanno.ilm10b4.hg19",unload = TRUE)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
chrnames <- Locations$chr[match(rownames(B),rownames(Locations))]
pos <- Locations$pos[match(rownames(B),rownames(Locations))]

dmrs <- bumps$table[order(bumps$table$fwerArea),]
dmrs$MD <- unlist(lapply(1:nrow(dmrs),function(x){
	if(dmrs$L[x] > 1){
		mean(colMeans(B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "M")]))-
		mean(colMeans(B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "F")]))
	} else {
		mean((B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "M")]))-
		mean((B[dmrs$indexStart[x]:dmrs$indexEnd[x],which(pd$predictedSex == "F")]))
	}
  }))
dmrs$width <- with(dmrs, end - start)

dmrs <- dmrs[order(dmrs$fwerArea, decreasing = TRUE),]
dmrs <- dmrs[order(dmrs$fwerArea, decreasing = FALSE),]
dmrs <- dmrs[which(dmrs$fwerArea < 0.1),]

B.term <- B
pd.term <- pd
B.term1 <- get(load("B_Trim1.rda"))
B.term1 <- B.term1[match(rownames(B.term),rownames(B.term1)),]
pd.term1 <- get(load("pd_Trim1.rda"))


dmrs.gr <- makeGRangesFromDataFrame(dmrs, keep.extra.columns=TRUE)
#What DMRS contain probes that distinguish cell type?
celltype$end <- celltype$start
celltype.gr <- makeGRangesFromDataFrame(celltype, keep.extra.columns = TRUE)

ov <- findOverlaps(dmrs.gr,celltype.gr)
table(queryHits(ov))

#Prep cell type data
B.celltype <- B.celltype[,-grep("enz",pd.celltype$CellType)]
pd.celltype <- pd.celltype[-grep("enz",pd.celltype$CellType),]
pd.celltype$CellType <- gsub(" cs","",pd.celltype$CellType)
pd.celltype$CellType <- factor(pd.celltype$CellType, 
	levels = c("Endothelial","Hofbauer","Stromal","Trophoblasts","Villi"))
pd.celltype$Sex <- factor(pd.celltype$Sex, levels = c("F","M"))

#DMRS with buffer 
buffersize <- 100
dmrs.gr.buffer <- dmrs.gr
start(dmrs.gr.buffer) <- start(dmrs.gr.buffer) - buffersize
end(dmrs.gr.buffer) <- end(dmrs.gr.buffer) + buffersize
ov.buffer <- findOverlaps(dmrs.gr.buffer, Locations.celltype)

#Example for Fig 5 (mydrm < -2), same code used for Supplementary Figures 18-21
#DMR 2 - Full Term
p1 <- methylationplot2(data.frame(dmrs[2,]),B.term,pd.term$predictedSex, ylab = "Placenta % DNAm at Full Term",
			colchoices = c("goldenrod","dodgerblue"), legendLabel = "Fetal Sex",
			leftBufferSize = 100, rightBufferSize = 100)	

mydmr <- 2
xmin <- start(dmrs.gr.buffer)[mydmr]
xmax <- end(dmrs.gr.buffer)[mydmr]
plotme.beta <- B.celltype[,pd.celltype$Trimester == "Third"]
plotme.pd <- pd.celltype[pd.celltype$Trimester == "Third",]
beta.chr <- plotme.beta[subjectHits(ov.buffer)[queryHits(ov.buffer) == 2],]
tempbeta <- stack(t(beta.chr))
variable <- plotme.pd$CellType
myshape <- plotme.pd$Sex
mypos <- start(Locations.celltype)[subjectHits(ov.buffer)[queryHits(ov.buffer) == 2]]
colchoices <- c("purple","blue","green","gold","firebrick")

df <- data.frame("BetaValues" = tempbeta[, 3], 
	"Var" = rep(variable, times = length(unique(tempbeta$col))), 
	"Shape" = rep(myshape, times = length(unique(tempbeta$col))),
	"Pos" = rep(mypos, each = ncol(beta.chr)))

out <- ggplot(df, aes(Pos, BetaValues, col = Shape)) + geom_point() + 
	theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), legend.position = "top", 
		panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5),
		axis.text.x = element_text(size = 12),strip.text.x = element_text(size = 12)) + 
	scale_x_continuous(limits=c(xmin,xmax))+
	xlab(paste0("Genomic Position on " , dmrs$chr[mydmr])) + ylab("Placenta Cell Type % DNAm at Full Term") + 
	scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),labels = seq(0,100,by = 25)) +
	scale_color_manual(guide = FALSE, name = "Fetal Sex", values = c("goldenrod","dodgerblue")) + stat_smooth(method="loess",se=FALSE) + 
	geom_vline(xintercept = dmrs$start[mydmr], linetype = "dashed") + geom_vline(xintercept = dmrs$end[mydmr], linetype = "dashed") +
	facet_grid(cols = vars(Var)) + 
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

png("CellType_DMR2_FullTerm.png", width = 12, height = 8, units = "in", res = 480)
grid.arrange(p1,out,nrow =2)
dev.off()

##############
#DMR 2 - 1st Trim
p1 <- methylationplot2(data.frame(dmrs[2,]),B.term1,pd.term1$predictedSex, ylab = "Placenta % DNAm at 1st Trimester",
			colchoices = c("goldenrod","dodgerblue"), legendLabel = "Fetal Sex",
			leftBufferSize = 100, rightBufferSize = 100)	

mydmr <- 2
xmin <- start(dmrs.gr.buffer)[mydmr]
xmax <- end(dmrs.gr.buffer)[mydmr]
plotme.beta <- B.celltype[,pd.celltype$Trimester == "First"]
plotme.pd <- pd.celltype[pd.celltype$Trimester == "First",]
beta.chr <- plotme.beta[subjectHits(ov.buffer)[queryHits(ov.buffer) == 2],]
tempbeta <- stack(t(beta.chr))
variable <- plotme.pd$CellType
myshape <- plotme.pd$Sex
mypos <- start(Locations.celltype)[subjectHits(ov.buffer)[queryHits(ov.buffer) == 2]]
colchoices <- c("purple","blue","green","gold","firebrick")

df <- data.frame("BetaValues" = tempbeta[, 3], 
	"Var" = rep(variable, times = length(unique(tempbeta$col))), 
	"Shape" = rep(myshape, times = length(unique(tempbeta$col))),
	"Pos" = rep(mypos, each = ncol(beta.chr)))

out <- ggplot(df, aes(Pos, BetaValues, col = Shape)) + geom_point() + 
	theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), legend.position = "top", 
		panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5),
		axis.text.x = element_text(size = 12),strip.text.x = element_text(size = 12)) + 
	scale_x_continuous(limits=c(xmin,xmax))+
	xlab(paste0("Genomic Position on " , dmrs$chr[mydmr])) + ylab("Placenta Cell Type % DNAm at 1st Trimester") + 
	scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),labels = seq(0,100,by = 25)) +
	scale_color_manual(guide = FALSE, name = "Fetal Sex", values = c("goldenrod","dodgerblue")) + stat_smooth(method="loess",se=FALSE) + 
	geom_vline(xintercept = dmrs$start[mydmr], linetype = "dashed") + geom_vline(xintercept = dmrs$end[mydmr], linetype = "dashed") +
	facet_grid(cols = vars(Var)) + 
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

png("CellType_DMR2_1stTerm.png", width = 12, height = 8, units = "in", res = 480)
grid.arrange(p1,out,nrow =2)
dev.off()

#############
#Supplemental Table 2 and Figure 6

library(Biobase)  
library(GEOquery) 
library(limma) 
library(dplyr) 
library(stringr) 
library(pheatmap) 
library(readr) 
library(limma) 
library(ggplot2) 
library(ggrepel) 

# load series and platform data from GEO. 'length(gset)' checks how many platforms used
gset_orig <- getGEO("GSE75010", GSEMatrix =TRUE) #Rows: 32321 Columns: 158 
gset <- gset_orig
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))
summary(exprs(gset)) 

## extract phenotype data
sampleInfo <- pData(gset)

## select column(s) from phenotype data that contain factor(s) we need for excluding samples 
sampleInfo <- select(sampleInfo,preeclampsia=characteristics_ch1,ga_week=characteristics_ch1.18,newborn_sex=characteristics_ch1.20) ## characteristics_ch1: diagnosis: PE /no-PE, 
sampleInfo$ga_week <- gsub("ga\\ \\(week\\)\\:\\ ","",as.character(sampleInfo$ga_week)) 
sampleInfo$ga_week <- as.numeric(as.character(sampleInfo$ga_week)) 
sampleInfo$newborn_sex <- gsub("infant\\ gender\\:\\ ","",as.character(sampleInfo$newborn_sex)) 

## create filter to exclude pre-eclampsia and preterm samples 
sml <- sampleInfo$newborn_sex ## group names 
sel <- which((sampleInfo$preeclampsia=="diagnosis: non-PE")&(sampleInfo$ga_week>36)) 
sml <- sml[sel] 
gset<-gset[ ,sel] 
sampleInfo<-subset(sampleInfo,(preeclampsia=="diagnosis: non-PE")& (ga_week>36))

## summary, Check the normalisation and scales used
summary(exprs(gset)) 

## A boxplot can also be generated to see if the data have been normalised. 
exprs(gset) <- log2(exprs(gset)) 
boxplot(exprs(gset),outline=FALSE) 

## view heatmap, look to see that groups as expected 
#library(pheatmap) 
## argument use="c" stops an error if there are any missing data points 
corMatrix <- cor(exprs(gset),use="c") 
anno_colors<-list(newborn_sex=c(M="turquoise", F="salmon"))
pheatmap(corMatrix,annotation=select(sampleInfo, newborn_sex), annotation_colors = anno_colors) 

## Print the rownames of the sample information and check it matches the correlation matrix 
rownames(sampleInfo) 
colnames(corMatrix) 
rownames(sampleInfo)==rownames(sampleInfo) 

full_output <- cbind(fData(gset),exprs(gset)) 
write_csv(full_output, file="gse_full_output.csv") 

features <- fData(gset) 

### Look at the features data frame and decide the names of the columns you want to keep 
features <- select(features,ID,SPOT_ID,seqname,gene_assignment,mrna_assignment) 
full_output <- cbind(features,exprs(gset)) 
write_csv(full_output, file="gse_full_output_v2_20220222.csv") ##file="gse_full_output_v2_20220218.csv")

## The design matrix is a matrix of 0 and 1s; one row for each sample and one column for each sample group. A 1 in a particular row and column indicates that a given sample (the row) belongs to a given group (column).
design <- model.matrix(~0+sampleInfo$newborn_sex) 
## rename column names 
colnames(design) <- c("F","M") 
design 

summary(exprs(gset)) 

## The lmFit function is used to fit the model to the data. The result of which is to estimate the expression level in each of the groups that we specified.
fit <- lmFit(exprs(gset), design) 
head(fit$coefficients) 

## define contrasts for DGE 
contrasts <- makeContrasts(M - F, levels=design) 
fit2 <- contrasts.fit(fit, contrasts) 

##apply the empirical Bayes’ step to get our differential expression statistics and p-values. 
fit2 <- eBayes(fit2) 

#look at the results by using the topTable command 
topTable(fit2) 

##to know how many genes are differentially-expressed overall, use the decideTests function. 
#decideTests(fit2) 
table(decideTests(fit2)) 

anno <- fData(gset) 
anno <- select(anno,ID,SPOT_ID,seqname,gene_assignment,mrna_assignment) 
anno$gene_transcript <- strsplit2(strsplit2(anno$gene_assignment," /// ")[,1]," // ")[,1]
anno$gene <- strsplit2(strsplit2(anno$gene_assignment," /// ")[,1]," // ")[,2]
anno$gene_alias <- strsplit2(strsplit2(anno$gene_assignment," /// ")[,1]," // ")[,3]
anno$gene_loc <- strsplit2(strsplit2(anno$gene_assignment," /// ")[,1]," // ")[,4]
anno$gene_misc <- strsplit2(strsplit2(anno$gene_assignment," /// ")[,1]," // ")[,5]
anno$ensembl_ids <- str_match_all(anno[,"mrna_assignment"], "gene:\\s*(.*?)\\s* gene_biotype")
anno<-subset(anno,select=-c(gene_assignment,mrna_assignment)) ## minus sign (-) to deselect

fit2$genes <- anno 
topTable(fit2,number=50) 

## for making volcano plot 
full_results <- topTable(fit2, number=Inf) 
full_results$log10adjPval <- -log10(full_results$adj.P.Val)
full_results$log10Pvalue <- -log10(full_results$P.Value)
full_results <- tibble::rownames_to_column(full_results,"ID_study") 

p_cutoff <- 0.05 
fc_cutoff <- 0.01 #1 

topN <- 10 
options(ggrepel.max.overlaps = Inf)

## volcano plot with labels, excluding genes on sex chromosomes
filter(full_results, (seqname != "chrX")&(seqname !="chrY")) %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%  
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, gene,"")) %>% 
  ggplot(aes(x = logFC, y = log10adjPval, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black") 

topN <- nrow(filter(full_results, (seqname != "chrX")&(seqname !="chrY")&(adj.P.Val < p_cutoff)&(abs(logFC) > fc_cutoff)))

ids_of_interest <- mutate(filter(full_results, (seqname != "chrX")&(seqname !="chrY")&(adj.P.Val < p_cutoff)&(abs(logFC) > fc_cutoff)), Rank = 1:n()) %>% 
  filter(Rank <= topN) %>% 
  pull(ID)

gene_names <- mutate(filter(full_results, (seqname != "chrX")&(seqname !="chrY")&(adj.P.Val < p_cutoff)&(abs(logFC) > fc_cutoff)), Rank = 1:n()) %>% 
  filter(Rank <= topN) %>% 
  pull(gene) 

gene_matrix <- exprs(gset)[,][as.character(ids_of_interest),]

## create heatmap for genome wide significantly expressed DEGs
# Specify colors
anno_colors<-list(newborn_sex=c(M="turquoise", F="salmon"))

pheatmap(gene_matrix,
         labels_row = gene_names,
         annotation=select(sampleInfo, newborn_sex),
         annotation_colors = anno_colors,
         scale="row")

#genes_meth_names<-c("ZNF175","ZNF300","C5orf63","CDKN1C","CROT","HOXA4","SNCA","ZBED9","ERCC6L2")
genes_meth_ensembl_ids<-c("ENSG00000105497","ENSG00000145908","ENSG00000164241","ENSG00000129757","ENSG00000005469","ENSG00000197576","ENSG00000145335","ENSG00000232040","ENSG00000182150")
##ZNF175 ENSG00000105497
##ZNF300 ENSG00000145908
##C5orf63 ENSG00000164241
##CDKN1C ENSG00000129757
##CROT  ENSG00000005469
##HOXA4 ENSG00000197576
##SNCA ENSG00000145335
##ZBED9 ENSG00000232040
##ERCC6L2 ENSG00000182150 

gene_ids_meth<-c()
genes_meth_names<-c()

for (genes_meth_ensembl_id in genes_meth_ensembl_ids){
  #print(paste("searching for", genes_meth_ensembl_id))
  for (i in 1:nrow(anno)){ 
    if (genes_meth_ensembl_id %in% anno$ensembl_ids[[i]]){
      print(paste(genes_meth_ensembl_id,i,anno$ID[i],anno$gene[i]))
      gene_ids_meth <- c(gene_ids_meth, anno$ID[i])
      genes_meth_names <- c(genes_meth_names, anno$gene[i])
    }
  }
}

gene_meth_matrix <- exprs(gset)[,][as.character(gene_ids_meth),]

anno_colors<-list(newborn_sex=c(M="turquoise", F="salmon"))
pheatmap(gene_meth_matrix,
         labels_row = genes_meth_names,
         annotation=select(sampleInfo, newborn_sex),
         annotation_colors = anno_colors,
         scale="row")

## getting adjusted p-values for Shan's genes
gene_meth_pvalues <- subset(full_results,ID %in% as.character(gene_ids_meth)) %>% select('gene','adj.P.Val')

## getting UN-adjusted p-values for Shan's genes
gene_meth_pvalues <- subset(full_results,ID %in% as.character(gene_ids_meth)) %>% select('gene','P.Value','logFC','ID')

p = gene_meth_pvalues[[2]]
gene_meth_pvalues$adj.P.Val<-p.adjust(p, method = "BH", n = length(p))

## sort by gene name
## Supplementary Table 2
gene_meth_pvalues[order(gene_meth_pvalues),][(nrow(gene_meth_pvalues[order(gene_meth_pvalues),])/length(gene_meth_pvalues[order(gene_meth_pvalues),])*(length(gene_meth_pvalues[order(gene_meth_pvalues),])-1) +1):nrow(gene_meth_pvalues[order(gene_meth_pvalues),]),]  #[13:24,] or #[37:48]

s_data <- sampleInfo ## originaly from pData(gset)
s_data <- cbind(geo_accession = rownames(s_data), s_data)
dim(s_data)

e_data <- full_output
e_data$gene <- strsplit2(strsplit2(e_data$gene_assignment," /// ")[,1]," // ")[,2]
e_data<-subset(e_data,select=-c(SPOT_ID,seqname,gene_assignment,mrna_assignment)) ## minus sign (-) to deselect columns that we want to exclude
dim(e_data)

## transpose so dimensions of e_data and s_data are the same
# e_data_subset <- e_data_subset %>% tidyr::gather(geo_accession, Expression,-ID) %>% 
#   tidyr::spread(ID, Expression)
# e_data_subset
e_data <- e_data %>% 
  dplyr::select(-gene) %>% 
  tidyr::gather(geo_accession, Expression,-ID) %>% 
  tidyr::spread(ID, Expression)
dim(e_data) ## should now have same # of rows as s_data

## join sample info with expression data
#s_data <- left_join(s_data, e_data_subset)
s_data <- left_join(s_data, e_data) #by geo_accession
dim(s_data)

## for ZNF300
#names(s_data)[names(s_data) == ""] <- "ZNF300"
ggplot(s_data, aes(x = newborn_sex, y =`8115196`,fill=newborn_sex),) + ylab("ZNF300")  + geom_boxplot()
ggsave("Figure6.png")

#####################

methylationplot2 <- function(bumpstab1, beta, variable, legendLabel = "Fetal Sex", 
	leftBufferSize = 100, rightBufferSize = 100, ylab = "DNAm", colchoices = c("goldenrod","dodgerblue")){

	chr <- Locations$chr[match(rownames(beta),rownames(Locations))]
	pos <- Locations$pos[match(rownames(beta),rownames(Locations))]
  if(!is.factor(variable)){
    variable <- factor(variable)}
  if(length(leftBufferSize) == 1){
    leftBufferSize <- rep(leftBufferSize, nrow(bumpstab1))
  }
  if(length(rightBufferSize) == 1){
    rightBufferSize <- rep(rightBufferSize, nrow(bumpstab1))
  }
  if(is.vector(leftBufferSize) == TRUE & length(leftBufferSize) != nrow(bumpstab1)){
    stop("Length of leftBufferSize does not equal number of rows in bumpstab1")
  }
  if(is.vector(rightBufferSize) == TRUE & length(rightBufferSize) != nrow(bumpstab1)){
    stop("Length of rightBufferSize does not equal number of rows in bumpstab1")
  }
  dmrgenes <- rep(NA,nrow(bumpstab1))
  
  n <- 1

    beta.chr <- beta[which(chr == bumpstab1$chr[n]),]
    pos.chr <- pos[which(chr == bumpstab1$chr[n])]
    #print(head(pos.chr))
    #print(class(pos.chr))
    startpos <- which((pos.chr >= bumpstab1$start[n] - leftBufferSize[n]) == TRUE)
    start <- startpos[1]
    endpos <- which((pos.chr <= bumpstab1$end[n] + rightBufferSize[n]) == TRUE)
    end <- endpos[length(endpos)]


    xmin <- bumpstab1$start[n] - leftBufferSize[n]
    xmax <- bumpstab1$end[n] + rightBufferSize[n]

    #message("HIII")
    #print(head(beta.chr[start:end,]))
    tempbeta <- stack(t(beta.chr[start:end, ]))
    #print(head(tempbeta))

    df <- data.frame("Basename" = rep(colnames(beta.chr), times = length(pos.chr[start:end])), 
    	"BetaValues" = tempbeta[, 3], 
    	"Var" = rep(variable, times = length(pos.chr[start:end])), 
    	"Chr" = rep(beta.chr[start:end], each = ncol(beta.chr)), 
    	"Pos" = rep((pos.chr[start:end]), each = ncol(beta.chr)))
  	message("Plotting DMR ", n)
    shiftme <- levels(df$Var)[2]
    df$Pos[which(df$Var == shiftme)] <- df$Pos[which(df$Var == shiftme)] + 0

    out <- ggplot(df, aes(Pos, BetaValues, col = Var)) + geom_point() + 
    	theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), legend.position = "top", 
    		panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5),
    		axis.text.x = element_text(size = 8)) + 
    	scale_x_continuous(limits=c(xmin,xmax))+
    	xlab(paste0("Genomic Position on " , bumpstab1$chr[n])) + ylab(ylab) + 
    	scale_y_continuous(breaks = seq(0,1,by = 0.25), limits = c(0,1),labels = seq(0,100,by = 25)) +
    	scale_color_manual(name = legendLabel, values = colchoices) + stat_smooth(method="loess",se=FALSE) + 
    	geom_vline(xintercept = bumpstab1$start[n], linetype = "dashed") + geom_vline(xintercept = bumpstab1$end[n], linetype = "dashed")

    return(out)
}

