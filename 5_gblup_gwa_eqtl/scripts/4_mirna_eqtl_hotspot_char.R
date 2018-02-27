#' **Script:** `4_mirna_eqtl_hotspot_char.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12/11/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects`
#' 
#' 2. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`
#'
#' **Input File(s):** 
#' 
#' 1. `6_mirna_precursor_annot_ssc11.Rdata`
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' 2. `4_normalized_dge_voom.Rdata`
#' 
#' 3. `5_Z_G_miRNA.Rdata`
#' 
#' 4. `MSUPRP_gpData_Ss11.Rdata`
#' 
#' **Output File Directory:** 
#'
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl`
#' 
#' **Output File(s):** 
#' 
#' 1. `6_hotspot_miRNA_corr_rst.Rdata`
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
#' The objective of this script is to investigate the putative miRNA eQTL hotspots to determine if they are truly hotspots, or spurious associations due to high correlation of miRNA expression or genotype data.
#' 
#' To do this, miRNA expression will be correlated between hotspot miRNAs utilizing the log-cpm (v$E) . 
#' Pearson correlation will be used. 
#' 
#' I will also investigate the genomic origins of these miRNAs and SNPs, to see if the miRNAs are coming from similar genomic regions, if they have similar seed sequences, etc.
#' 
#' ## Install libraries

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

rm(list=ls())

library(regress)
library(gwaR)
library(limma)
library(edgeR)
library(parallel)
library(qvalue)
library(corrplot)

#' ## Load data
#' 
#' Load microRNA expression data
load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
#' Load the MSUPRP gpdata object (for allele freq calculation):
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")

ls()

#' ## Analysis
#' 
#' Extract the names of the miRNA in the hotspots:
#' 
#' Two SNP ("MARC0027291" "MARC0093624") were associated with five miRNAs, overlapping to some degree with eachother:
htspt.mirna5 <- as.character(c("ssc-let-7d-5p","ssc-let-7g","ssc-miR-1468","ssc-miR-345-3p","ssc-miR-95","ssc-miR-9843-3p"))
htspt.mirna5

#' Correlate the miRNA expression profiles using the log-cpm counts of miRNA expression:
data5<-t(v$E[htspt.mirna5,])
head(data5)
cormx5<-cor(data5)
head(cormx5)

#' Function to perform correlation analysis of miRNA to miRNA espression and export a summary of correlation analysis
cor.exp <- function(var1, var2, data, ...){
	x <- cor.test(data[,as.character(var1)], data[,as.character(var2)], ...)
	rst <- data.frame(var1, var2, cor=x$estimate, 
		t=x$statistic, 
		pvalue=x$p.value)
	return(rst)
}

#' ---
#' 
#' Correlation test of 5 miRNAs:
thres<-0.05
mir5.htspt<-list()
for(i in 1:length(colnames(data5))){
mir5.htspt[[i]]<-data.frame(var1=colnames(data5)[i], var2=colnames(data5)[-i])
}

vars5<-do.call(rbind, mir5.htspt)

mir.cor5<- do.call(rbind, mapply(cor.exp, var1=vars5[,1], var2=vars5[,2], MoreArgs=c(list(data=data5), alternative="two.sided", method="p"), SIMPLIFY=FALSE))

#' Perform multiple test correction (FDR) on correlation analysis:
mir.cor5$qval<-qvalue(mir.cor5$pvalue)$qvalue
mir.cor5$pi0<-qvalue(mir.cor5$pvalue)$pi0
mir.cor5

#' Significantly correlated miRNA pairs:  
mir.cor5[mir.cor5$qval<thres,]

#' ---
#' 
#' Investigate correlation of genotypes between the 4 putative eQTL hotspots:
snp.htspt<- c("MARC0027291", "MARC0093624")

#' Examine allele frequencies of the hotspot SNPs
#' 
#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-MSUPRP$geno[!((as.numeric(rownames(MSUPRP$geno))<=1000) | (as.numeric(rownames(MSUPRP$geno))>=6000)),]
dim(geno_f2)

#' Subset the hotspot SNPs:
geno.htspt<-geno_f2[,snp.htspt]
dim(geno.htspt)

#' Subset the 174 animals:
geno.htspt<-geno.htspt[rownames(data5),]
dim(geno.htspt)
head(geno.htspt)

if(sum(rownames(geno.htspt)!=rownames(data5)) != 0) stop ("Animal IDs not correct for genotype object")

geno.cor<-cor(geno.htspt)

#' 
#' Calculate allele frequency for all the F2 animals:
allele_freq<-colMeans(geno.htspt,na.rm=T)/2
allele_freq

#' Minor allele frequency:
maf<-ifelse(allele_freq>0.5, (1- allele_freq), allele_freq)
maf

#' ---
#' 
#' Identify the associated miRNAs in the annotation file, to query miRBase and the literature for information on their effects 
head(annotation)
prec.annot<-annotation[annotation$miRNA==c("ssc-let-7d-5p", "ssc-let-7g","ssc-miR-1468", "ssc-miR-345-3p", "ssc-miR-95", "ssc-miR-9843-3p"),]

#' Information extracted from target prediction input:
#' 
#' miRNA           | Seed Sequence (Human)
#' --------------- | -------------
#' let-7-5p/98-5p  | GAGGUAG
#' miR-1468-5p	   | UCCGUUU
#' miR-345-3p	   | CCCUGAA
#' miR-95-3p	   | UCAACGG
#' 
#' Also see `miRNA_targets_common-1.tiff` to compare the targets in common expressed in this dataset between miRNAs.
#' 
#' ## Visualize
#' 
#' Correlation plots of 4 and 5 hotspot-associated miRNAs:
#' 
#' 5 miRNAs:
#+ miRNA5_hotspots, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot.mixed(cormx5, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Five hotspot-associated miRNAs")

#+ geno_cor_hotspots, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot.mixed(geno.cor, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Putative Hotspot SNP Genotype Correlation")

#' ## Save data
save(mir.cor5, geno.cor, maf, file="../6_hotspot_miRNA_corr_rst.Rdata")