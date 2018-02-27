#' **Script: ** `3_mirna_eQTL_characterization.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  `12/7/17`
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/eQTL/paper/output/`
#' 
#' 3. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 2. `MSUPRP_miRNA.Rdata`
#' 
#' 3. `funct_eqtl.Rdata`
#'
#' 4. `6_mirna_precursor_annot_ssc11.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_mirna_names.txt`
#' 
#' 2. `4_miRNA_eQTL_local_distant_regulators.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to characterize the miRNA eQTL results. This includes:
#'
#' Using the summary file from the miRNA eQTL analysis determine the number of miRNA eQTL peaks
#' that have local regulation and distant regulation.
#' 
#' Looking only at the eQTL peaks on the same chromosome or overlapping with the mapped gene position check how many markers are before and after the gene position.
#' 
#' Create a data frame containing the results of the eQTL scan along with a column
#' specifying the eQTL regulator type, local (cis), distant (tran)
#' 
#' This analysis will be completed using functions created by DV in directory /mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z`
#' 
#' ## Install libraries

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' ## Load data
#'
rm(list=ls())
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

#' Load DV's functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

#' Load the gwa results:
load("../2_gwa_results.Rdata")

#' Load the miRNA eQTLsummary tables:
load("../3_eqtl_summary_tables_maps.Rdata")

#' Load the MSUPRP_miRNA gpdata object for the mapping information:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")

#' Load the miRNA annotation files to obtain info for miRNA eQTL with multiple precursors:
load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")

#' Load the MSUPRP gpdata object:
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
ls()

#' ## Analysis

#' eQTL peaks
head(sum.eqtl)
str(sum.eqtl)

#' ---
#' 
#' Separate the first column for use in the Target Prediction analysis, substituting the "ssc" for "hsa":
hsamir<-data.frame(V1=unique(gsub("ssc-", "", sum.eqtl$miRNA)))
hsamir

#' Save that table for use in target prediction analysis:
write.table(hsamir, file="../../6_mirna_eQTL_target_prediction/1_mirna_names.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#' ---
#' 
#' eQTL peaks with unknown position
as.character(sum.eqtl$miRNA[is.na(sum.eqtl$chr.miR)])

#' Hotspots
tail(sort(table(as.character(fullsum.eqtl$SNP))))
mk <- table(as.character(fullsum.eqtl$SNP))

table(mk)
#' Names of hotspot SNP (SNP associated with at least 4 miRNAs)
nm <- names(mk[mk >= 4])
nm

#' Check map positions of hotspot SNP:
MSUPRP_miRNA$map[nm,]
#' All my hotspot SNP are on chromosome 15; LD?
abspos<-absposmap
abspos[nm]


#' Want to extract the miRNA information for each potential hotspot SNP
htsp.mir<-list()
for(i in nm){
htsp.mir[[i]]<-fullsum.eqtl[grep(i, as.character(fullsum.eqtl$SNP)), c("miRNA", "chr.miR")]
}

htsp.mir

#' Build a matrix of qvalues for all the miRNAs.
#' 
#' This will be used in plotting the eQTL map and determining local vs distal-acting miRNA eQTL
sigmirqval<-do.call(cbind, lapply(unique(as.character(sum.eqtl$miRNA)), function(x) rst.gwa[rst.gwa$miRNA==x, "gwa.qval"]))

colnames(sigmirqval)<-unique(as.character(sum.eqtl$miRNA))
rownames(sigmirqval)<-rownames(MSUPRP_miRNA$map)
dim(sigmirqval)
head(sigmirqval)
str(sigmirqval)

#' Check for correctness in building this matrix:
names(rst.gwa)
head(rst.gwa)
#' Rownames didn't change between datasets:
if (sum(unique(rst.gwa$SNPid)!=rownames(sigmirqval)) !=0) stop("Rownames changed between datasets")
#' Print the sum of the q-values not equal between the gwa results and the subset matrix:
for (i in colnames(sigmirqval)){
	print(sum(rst.gwa[rst.gwa$miRNA==i,"gwa.qval"] != sigmirqval[,i]))
}
#' Looks good.

#' Extract the q-values for 17 significant miRNA associations with hotspot SNPs:
htsp.qv <- sigmirqval[nm,]
dim(htsp.qv)
head(htsp.qv)


#' ---
#'
#' Table of positions for significant miRNA-marker associations within a hotspot using gene-wise FDR
#' 
#' tbpos = table of positions; takes the qvalues and absolute positions of hotspot SNPs & miRNAs (htsp.qv), 
#' and puts them into a table for use in the eQTL map function later on.
threshold <- 0.05
#' ttp.htsp = table to plot. hotspot
ttp.htsp <- tbpos(qval=htsp.qv, abspos=abspos, threshold=threshold)
dim(ttp.htsp)
ttp.htsp

#' Table of positions for all gene-marker associations and for significant eQTL SNPs
#' 
#' Notice, it removed the miRNAs that don't have assembly information, miR-6782-3p and miR-9785-5p.
#' 
#' It can't add it to the eQTL map if it doesn't have mapping information. 
ttp.all <- tbpos(qval=sigmirqval, abspos=abspos, threshold=threshold)
dim(ttp.all)
head(ttp.all)
table(as.character(ttp.all$Gene))

#' Also lose 21 SNPs, which are the ones significantly associated to miR-6782-3p and miR-9785-5p:
table(sigmirqval[,"ssc-miR-6782-3p"]<threshold)
names(which(sigmirqval[,"ssc-miR-6782-3p"]<threshold))
table(as.character(ttp.all$SNP) %in% names(which(sigmirqval[,"ssc-miR-6782-3p"]<threshold)))

table(sigmirqval[,"ssc-miR-9785-5p"]<threshold)
names(which(sigmirqval[,"ssc-miR-9785-5p"]<threshold))
table(as.character(ttp.all$SNP) %in% names(which(sigmirqval[,"ssc-miR-9785-5p"]<threshold)))

#' Build a data.frame of the miRNA eQTL peak positions and their lengths
ttp.peaks <- data.frame(miRNA=sum.eqtl[,"miRNA"], SNP=sum.eqtl[,"SNP"],
	pos.snp=abspos[as.character(sum.eqtl$SNP)],
	pos.miR=abspos[as.character(sum.eqtl$miRNA)],
	diff=abs(abspos[as.character(sum.eqtl$miRNA)] - abspos[as.character(sum.eqtl$SNP)]),
	qvalue=sum.eqtl[,"qvalue"])
head(ttp.peaks)
ttp.peaks

#' Local eQTL (na.omit to remove the miRNA with no map position data)
local <- na.omit(mirpeaks[mirpeaks$chr.miR == mirpeaks$chr.snp,])
dim(local)
local

#' Distant eQTL (na.omit to remove the miRNA with no map position data)
distant <- na.omit(mirpeaks[mirpeaks$chr.miR != mirpeaks$chr.snp,])
dim(distant)
distant

#' Compute the distance between the mapped position of the gene expression and the position of the peak
#+ results='hide'
distancemir<-function(peaks){
    dist <- data.frame(start.min=peaks$start.miR - peaks$min.pos, max.end=peaks$max.pos - peaks$end.miR, diff=abs(peaks$pos.snp - ((peaks$end.miR + peaks$start.miR)/2)))
    return(dist)
}

#' Identify the local regulators:
distL <- distancemir(local)
dim(distL)
distL

#' Make sure the correct amount of miRNAs are in the "distant" category (distD not used downstream):
distD <- distancemir(distant)
dim(distD)
distD

#' miRNA expression mapped within the eQTL peak
cisI <- local[distL$start.min > 0 & distL$max.end > 0,]
nrow(cisI)
cisI

#' eQTL peak within the mapped position of the miRNA expression (only one marker within peak)
cisII <- local[distL$start.min < 0 & distL$max.end < 0,]
nrow(cisII)

#' miRNA expressions mapped in close proximity to the eQTL peak (less than 10MB)

# miRNA maps to the right hand side of its peak
cisIIIa <- local[distL$start.min > 0 & distL$max.end < 0,]
cisIIIa <- data.frame(cisIIIa, dist=abs(ifelse(cisIIIa$max.pos - cisIIIa$start.miR < 0,
	cisIIIa$max.pos - cisIIIa$start.miR, cisIIIa$pos.snp - cisIIIa$start.miR)),
	contained=ifelse(cisIIIa$max.pos - cisIIIa$start.miR < 0, "No", "Yes"))
cisIIIa<-data.frame(cisIIIa,position=rep("right",nrow(cisIIIa)))
nrow(cisIIIa)
cisIIIa

# miRNA maps to the left hand side of its peak
cisIIIb <- local[distL$start.min < 0 & distL$max.end > 0,]
cisIIIb <- data.frame(cisIIIb, dist=abs(ifelse(cisIIIb$end.miR - cisIIIb$min.pos < 0,
	cisIIIb$end.miR - cisIIIb$min.pos, cisIIIb$end.miR - cisIIIb$pos.snp)),
	contained=ifelse(cisIIIb$end.miR - cisIIIb$min.pos < 0, "No", "Yes"))
cisIIIb<-data.frame(cisIIIb,position=rep("left",nrow(cisIIIb)))
nrow(cisIIIb)

#' miRNA overlapping peak region
cisIII <- rbind(cisIIIa,cisIIIb)
cisIII <- cisIII[cisIII$contained == "Yes",]
nrow(cisIII)

#' miRNAs on the same chromosome as peak
cisIV <- rbind(cisIIIa,cisIIIb)
cisIV <- cisIV[!cisIV$contained == "Yes",]
nrow(cisIV)
cisIV

#' eQTL mapped less than 5Mb from miRNA
idx <- abs(cisIV$dist) < 5e6
nrow(cisIV[idx,])
cisIV

#' eQTL mapped on same chromosome but at a distance greater than 5MB
idx <- abs(cisIV$dist) > 5e6
nrow(cisIV[idx,])
tranI <- cisIV[idx,]
cisIV <- cisIV[!idx,]
tranI

#' eQTL mapped on a different chromosome than the mapped position of the miRNA expression
tranII <- distant[distant$range.peak > 0 & distant$chr.miR %in% 1:19,]
nrow(tranII)
tranII

#' eQTL mapped on a different chromosome than the mapped position of the miRNA expression but with only one marker in peak
tranIII <- distant[distant$range.peak == 0 & distant$chr.miR %in% 1:19,]
nrow(tranIII)
tranIII

#' ---
#'  
#' Function to plot eQTL peaks
plot.GMA <-
function (ttp, abspos, ...) {
    plot(ttp[,"pos.snp"],ttp[,"pos.miR"],ylim=c(min(abspos),
         max(abspos)),xlim=c(min(abspos),max(abspos)), pch=1,
         cex.main=2,cex.lab=1.3,
         cex.axis=1.2, ...)
    abline(0,1, col="blue", lwd=2)
}

#' eQTL Plot
#' 
#' Change map positions to Mb for mapping
absposmb<-abspos/1e6
head(absposmb)
ttp.peaks$pos.snp<-ttp.peaks$pos.snp/1e6
ttp.peaks$pos.miR<-ttp.peaks$pos.miR/1e6
head(ttp.peaks)
dim(ttp.peaks)

#+ eQTLmap-regulators, fig.align='center'
par(oma=c(5,2,2,2))
plot.GMA(ttp=ttp.peaks, abspos=absposmb, xlab="SNP position Mb",
		ylab="Transcript position Mb", main="miRNA eQTL Map")

# miRNA expression mapped within the eQTL peak
points(ttp.peaks[match(cisI$SNP,as.character(ttp.peaks$SNP)),"pos.snp"],
	ttp.peaks[match(cisI$SNP,as.character(ttp.peaks$SNP)),"pos.miR"],
	col="black", pch=20)

# miRNAs on the same chromosome as peak
points(ttp.peaks[match(cisIV$SNP,as.character(ttp.peaks$SNP)),"pos.snp"],
	ttp.peaks[match(cisIV$SNP,as.character(ttp.peaks$SNP)),"pos.miR"],
	col="yellow", pch=20)

# miR-eQTL mapped on a different chromosome than the mapped position of the miRNA expression (> 1 SNP in peak)
points(ttp.peaks[match(tranII$miRNA,as.character(ttp.peaks$miRNA)),"pos.snp"], ttp.peaks[match(tranII$miRNA, as.character(ttp.peaks$miRNA)),"pos.miR"],
	col="blue", pch=20)

# miR-eQTL mapped on a different chromosome than the mapped position of the miRNA expression but with only one marker in peak
points(ttp.peaks[match(tranIII$SNP,as.character(ttp.peaks$SNP)),"pos.snp"], ttp.peaks[match(tranIII$SNP,as.character(ttp.peaks$SNP)),"pos.miR"],
	col="forestgreen", pch=20)

# Putative hotspot SNP
points(ttp.peaks[as.character(ttp.peaks$SNP) %in% nm,"pos.snp"], ttp.peaks[as.character(ttp.peaks$SNP) %in% nm,"pos.miR"],
	col="red", pch=20)


add_legend(-0.75,-0.75, legend=c("Overlapping", "Same Chr <5Mb",
	"Different Chr", "Single SNP Different Chr", "Putative Hotspot"),
	pch=c(19,19,19,19,19), bty="o",ncol=2, cex=1.2)
	add_legend(-0.75,-0.75, legend=c("Overlapping", "Same Chr <5Mb",
	"Different Chr", "Single SNP Different Chr", "Putative Hotspot"),
	pch=c(20,20,20,20,20), col=c("black", "yellow",
	"blue", "forestgreen", "red"), bty="n",ncol=2, cex=1.2)


#' Save maps for markers_around_peaks.R analysis
#' 
#' Complete map (unfiltered SNPs; X=19, Y=20)
comp.snp.map <- MSUPRP$map
comp.snp.map$chr<-as.character(comp.snp.map$chr)
comp.snp.map$chr<-gsub("X", "19", comp.snp.map$chr)
comp.snp.map$chr<-gsub("Y", "20", comp.snp.map$chr)

dim(comp.snp.map)
table(comp.snp.map$chr)

#' Reduced map (SNPs used in this eQTL analysis)
snp.map <- MSUPRP_miRNA$map
snp.map$chr<-as.character(snp.map$chr)
dim(snp.map)
head(snp.map)
table(snp.map$chr)

#' Change the scale of the snp position to Mb
snp.map$pos <- snp.map$pos /1e6
comp.snp.map$pos <- comp.snp.map$pos /1e6

#' ----
#' 
#' Looking only at the eQTL peaks on the same chromosome or overlapping with the mapped gene position check how many markers are before and after the gene position.
#' 
#' Number of snp per chromosome
# Complete map
nsnpall <- c(table(as.numeric(comp.snp.map[,"chr"])))
nsnpall

# Filtered map
nsnpred <- c(table(as.numeric(snp.map[,"chr"])),0,0)
names(nsnpred) <- 1:20
nsnpred
# Percent
perc <- round(nsnpred / nsnpall * 100)
perc

#' Plot number of snp per chromosome
# + snpmap, fig.align='center', fig.width=16, fig.height=12
nsnp <- rbind(nsnpred, nsnpall=nsnpall - nsnpred)
bp <- barplot(nsnp, col=c("blue","black"), ylim=c(0,max(nsnpall)+10),
	main="SNP per Chromosome", xlab="chromosome", ylab="#", cex.main=2,
	cex.axis=1.5, cex.lab=1.5, legend=c("Retained", "Removed"))
text(bp, nsnp[1,]-130, labels=nsnp[1,], col="white", cex=0.75)
text(bp, colSums(nsnp)-130, labels=nsnp[2,], col="white", cex=0.75)
text(bp[1:19], 100, labels=perc[1:19], col="white", cex=0.75)
text(bp[20], 100, labels=perc[20], col="black", cex=0.75)

#' Check the position of miR-140-5p eQTL that is potentially local-acting:
mirpeaks[mirpeaks$miRNA=="ssc-miR-140-5p",]

mirpeaks[mirpeaks$miRNA=="ssc-miR-140-5p","max.pos"]<mirpeaks[mirpeaks$miRNA=="ssc-miR-140-5p","start.miR"]
#' The SNPs in the peak end before the miRNA transcript starts. 
#' Check how many SNPs lie between the end of the eQTL peak and the start of the miRNA:
mirpeaks[6,"max.pos"]
mirpeaks[6,"start.miR"]

mirpeaks[6,"start.miR"]-mirpeaks[6,"max.pos"]

map.full["ssc-miR-140-5p",]

map.full["ALGA0117081",]

ssc6<-map.full[map.full$chr=="6",]
dim(ssc6)

ssc6<-ssc6[-c(2283, 2284, 2285, 2286),]
ssc6<-ssc6[order(ssc6$pos),]

ssc6[ssc6$pos>mirpeaks[6,"max.pos"] & ssc6$pos<mirpeaks[6,"start.miR"],]


ssc6[310:330,]

#' So only two SNPs lie between the miRNA and the eQTL peak, covering 100kb of the genome. 
#' 
mirpeaks[6,"start.miR"]-mirpeaks[6,"max.pos"]

#' 
#' Investigate correlation of genotypes between the SNPs between the miR-eQTL and the miRNA:
snp.btwn<- c("ASGA0102265", "MARC0033972", "ALGA0117081")

#' Examine allele frequencies of the hotspot SNPs
#' 
#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-MSUPRP$geno[!((as.numeric(rownames(MSUPRP$geno))<=1000) | (as.numeric(rownames(MSUPRP$geno))>=6000)),]
dim(geno_f2)

#' Subset the hotspot SNPs:
geno.btwn<-geno_f2[,snp.btwn]
dim(geno.btwn)

#' Subset the 174 animals:
geno.btwn<-geno.btwn[rownames(MSUPRP_miRNA$geno),]
dim(geno.btwn)
head(geno.btwn)

if(sum(rownames(geno.btwn)!=rownames(MSUPRP_miRNA$geno)) != 0) stop ("Animal IDs not correct for genotype object")

geno.cor<-cor(geno.btwn)
geno.cor

dim(rst.gwa)
head(rst.gwa)
mir1405p<-rst.gwa[rst.gwa$miRNA=="ssc-miR-140-5p",] 
dim(mir1405p)
mir1405p[mir1405p$SNPid=="ASGA0102265",]
mir1405p[mir1405p$SNPid=="MARC0033972",]
mir1405p[mir1405p$SNPid=="ALGA0117081",]

#' --------
#' 
#' Classify eQTL Peaks
#' 
#' Local eQTL (same chr)
cis <- rbind(cisI, cisII, cisIII[cisIII$contained == "Yes",1:13])
cis$miRNA <- as.character(cis$miRNA)
cis$SNP <- as.character(cis$SNP)
nrow(cis)
cisIII <- cisIII[!cisIII$contained == "Yes",]
cisIII
cis <- data.frame(cis, regulator=rep("cis", nrow(cis)))
cis

#' Distant eQTL
# Same chromosome, < 5Mb
transc<-cisIV[,1:14]
# Different chromosome
tran <- rbind(transc, tranII, tranIII)
tran$miRNA <- as.character(tran$miRNA)
tran$SNP <- as.character(tran$SNP)
nrow(tran)

tran <- data.frame(tran, regulator=c("transc.5Mb",rep("tran", nrow(tran)-1)))
tran

regul <- rbind(cis, tran)
rownames(regul) <- NULL

#' Add eQTL peaks lacking assembly information:
regul <- rbind(regul, data.frame(mirpeaks[!mirpeaks$miRNA %in% regul$miRNA,],
	regulator=rep(NA, sum(!mirpeaks$miRNA %in% regul$miRNA))))
regul
#' eQTL plot with updated classifications:
#+ miRNA_eQTL_local_distant, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
par(oma=c(5,2,2,2))
plot.GMA(ttp=ttp.peaks, abspos=absposmb, xlab="SNP position Mb",
		ylab="Transcript position Mb")

points(ttp.peaks[rownames(cis),"pos.snp"],
	ttp.peaks[rownames(cis),"pos.miR"], pch=20, col="white")

points(ttp.peaks[rownames(tran),"pos.snp"], ttp.peaks[rownames(tran),"pos.miR"],
	col="forestgreen", pch=20)

add_legend(-.25,-0.75,legend=c("Local",
	"Distant"),
	pch=c(1,1), bty="o",ncol=2, cex=1.2)
	add_legend(-0.25, -0.75,legend=c("Local",
	"Distant"),
	pch=c(1,19), 
	col=c("black","forestgreen"), bty="n",ncol=2, cex=1.2)

#' 
#' ## Save data
#' 
#' Save the data frame containing all the regulators, defining them as cis or trans (local or distant)
save(regul, file="../4_miRNA_eQTL_local_distal_regulators.Rdata")