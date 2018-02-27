#' **Script:** `5_target_mrna-pqtl-coloc.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  12/13/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/`, `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k`, `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z`
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`
#' 
#' **Input File(s):** 
#' 
#' 1. `voom.Rdata`, `pQTL_60k.Rdata`, `inrange_function.Rdata`
#' 
#' 1. `10_mrna_mirna_corr_rst.Rdata`
#' 
#' 1. `MSUPRP_gpData_Ss11.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Output File(s):** ``
#' 
#' 1. `12_target_mrna_coloc_pqtl.Rdata`, `13_target_mrna_coloc_pqtl.txt`, `14_target_mrna_neg_coloc_pqtl.txt`
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
#' The objective of this script is to co-localize the significant target mRNAs with pQTL identified in the MSUPRP dataset
#' 
#' ## Install libraries
#' 
rm(list=ls())
library(methods)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' Session Information
sessionInfo()
#' ## Load data
#' 
#' Load required R objects 
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/10_mrna_mirna_corr_rst.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")
ls()

#' ## Analysis
#' 
#' ### Annotation of corrleated target genes
genes <- dge$genes
colnames(genes)[1] <- "chr"
annot <- do.call(rbind, lapply(names(sig.mrnaR), function(x) 
	data.frame(genes[sig.mrnaR[[x]],], miRNA=rep(x, length(sig.mrnaR[[x]])), 
		rst.corR[[x]][sig.mrnaR[[x]],c("cor", "pvalue", "qvalue")])))
annot$ID<-as.character(rownames(annot))
head(annot)

#' ### Phenotypic QTL MSUPRP
length(QTL)
names(QTL)
QTL[[1]]

#' pQTL peak information
dim(QTLpeaks)
QTLpeaks$pheno.chr<-as.factor(make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE))
rownames(QTLpeaks)<-make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE)
head(QTLpeaks)
table(QTLpeaks$pheno)
table(QTLpeaks$chr)

#' q-values pQTL GWAS
qval <- GWAS$qvalue
dim(qval)

#' p-values pQTL GWAS
pval <- GWAS$pvalue
dim(pval)

#' Standardized SNP effects pQTL GWA
sdEff <- GWAS$Estimate
dim(sdEff)

map<-MSUPRP$map
str(map)
map$chr<-gsub("X", "19", map$chr)
map$chr<-gsub("Y", "20", map$chr)
map$chr<-as.numeric(map$chr)
str(map)

#' Retain marker information for all pQTL
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(map[names(sig[[x]]),], 
	std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]
length(sig)
names(sig)

#' If a pQTL contains more than one peak split each peak
idx <- lapply(sig, function(x) as.numeric(names(table(x$chr))))
idx <- idx[unlist(lapply(idx, function(x) length(x) > 1))]

mp <- lapply(1:length(idx), function(x) lapply(1:length(idx[[x]]), function(y) sig[[names(idx)[x]]][sig[[names(idx)[x]]]$chr == idx[[x]][y],]))
names(mp) <- names(idx)

for (i in 1:length(mp)){
names(mp[[i]]) <- paste(names(mp)[[i]], 1:length(mp[[i]]), sep=".")
}

qtl <- sig[!names(sig) %in% names(mp)]
for (i in 1:length(mp)){
qtl <- c(qtl, mp[[i]])
}

#' pQTL genomic regions
qtlP <- do.call(rbind, lapply(names(qtl), function(x) data.frame(pheno=strsplit(x, "[.]")[[1]][1], 
	chr=unique(qtl[[x]]$chr), start=min(qtl[[x]]$pos), end=max(qtl[[x]]$pos))))
rownames(qtlP) <- names(qtl)
qtlP <- qtlP[order(qtlP$chr),]
tmp <- list()
for(i in unique(qtlP$chr)){
	x <- qtlP[qtlP$chr == i,]
	tmp[[i]] <- x[order(x$start),]
}
qtlP <- do.call(rbind,tmp)
dim(qtlP)


#' ### Co-localization mRNA with pQTL 
#' 
#' Check all mRNA within a eQTL
win <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single=NULL, range=c(start="start", end="end")))
names(win) <- rownames(qtlP)


#' mRNA overlaping left side of pQTL
left <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="start", range=NULL))
names(left) <- rownames(qtlP)

#' mRNA overlaping right side of pQTL
right <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="end", range=NULL))
names(right) <- rownames(qtlP)

#' Merge all mRNA co-localizing with pQTL  
# Merge within and left side
coloc <- lapply(names(win), function(x) rbind(win[[x]], 
	left[[x]][!as.character(left[[x]]$geneID) %in% as.character(win[[x]]$geneID) | 
	!as.character(left[[x]]$miRNA) == as.character(win[[x]]$miRNA),]))
names(coloc) <- names(win)

names(coloc[c(1,28)])
#' For car_bf10.1 and num_ribs, there is one more right-overlapping peak than there is left-overlapping. 
#' So, when building the merged object, the two objects will have different lengths, triggering the warning 
#' and adding the last line of the longer list to the merged object.
#' 
# Merge coloc and right side
coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]], 
	right[[x]][!as.character(rownames(right[[x]])) %in% as.character(rownames(coloc[[x]])) |
	!as.character(right[[x]]$miRNA) == as.character(coloc[[x]]$miRNA),]))
names(coloc) <- names(win)

data.frame(left=sapply(1:length(left), function(x) nrow(left[[x]])), 
	right=sapply(1:length(right), function(x) nrow(right[[x]])),
	win=sapply(1:length(win), function(x) nrow(win[[x]])), 
	coloc=sapply(1:length(coloc), function(x) nrow(coloc[[x]])))

#' Final list of mRNA targets significantly correlated with miRNA and co-localizing with a pQTL
coloc <- do.call(rbind, lapply(names(coloc), function(x) data.frame(coloc[[x]], 
	pheno=rep(strsplit(x, "[.]")[[1]][1], nrow(coloc[[x]])))))
rownames(coloc) <- NULL
dim(coloc)
head(coloc)
table(coloc$miRNA)

coloc[coloc$miRNA=="ssc-miR-874",c("genes", "miRNA", "qvalue", "pheno")]

#' Investigate the occurences of multiple rows:
sort(table(as.character(coloc$genes)))


coloc$cor.mag<-ifelse(coloc$cor<0, "-", "+")

#' Convert miRNA names to human equivalent (based on seed sequence) for use in IPA:
hsa.coloc<-coloc

miRid<-as.character(hsa.coloc$miRNA)
#' Convert ssc to hsa
miRid<-gsub("ssc", "hsa", miRid)
#' ssc-miR-6782 is hsa-miR-6821
miRid<-gsub("6782", "6821", miRid)
#' ssc-miR-874 is hsa-miR-874-3p
miRid<-gsub("-874", "-874-3p", miRid)
#' ssc-miR-95 hsa-miR-95-3p
miRid<-gsub("-95", "-95-3p", miRid)

hsa.coloc$hsa.miRNA<-miRid

#' Extract the pertinent columns for the significant negatively-associated mRNAs overlapping pQTLs:
negcoloc<-coloc[coloc$cor < 0,c("chr", "ID", "genes", "miRNA", "cor", "pheno", "cor.mag")]
negcoloc
as.character(unique(negcoloc$genes))
as.character(unique(negcoloc$pheno))
table(negcoloc$miRNA)
table(as.character(negcoloc$pheno))

table(as.character(negcoloc[negcoloc$miRNA=="ssc-miR-874" & negcoloc$chr=="2","pheno"]))

#' In my previous analysis, I had significant co-localization in 3 miR-6782-3p target genes in num_ribs traits
#' 
#' These three genes (STOML1, KHNYN, and RGMA) are present in this analysis, but their correlations did not reach 
#' statistical significance (q-val between 0.06 and 0.15).

#' ## Save data
save(annot, coloc, negcoloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/13_target_mrna_coloc_pqtl.Rdata")
write.table(hsa.coloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/14_target_mrna_coloc_pqtl.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(negcoloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/15_target_mrna_neg_coloc_pqtl.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)