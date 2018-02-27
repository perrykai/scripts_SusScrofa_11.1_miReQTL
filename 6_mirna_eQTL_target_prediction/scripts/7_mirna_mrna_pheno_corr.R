#' **Script:** `7_mirna_mrna_pheno_corr.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  12/21/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11/`
#' 
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`, `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`, `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' **Input File(s):** 
#' 
#' 1. `voom.Rdata`, `MSUPRP_gpData_Ss11.Rdata`, `G_Z_matrix_Ss11.Rdata`
#' 
#' 2. `5_filtered_targets_exp_rst.Rdata`, `13_target_mrna_coloc_pqtl.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`
#' 
#' 3. `12_mrna_mirna_resid_exp.Rdata`, `9_mireqtl_pqtl_coloc_peaks.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Output File(s):** ``
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
#' The objective of this script is to correlate residual miRNA expression and residual mRNA expression with the phenotypes of interest
#' 
#' (meaning pQTL phenotypes co-localizing with miR-eQTL target genes)
#' 
#' These results will be utilized for supporting miR-874 co-localized targets results.
#' 
#' ## Install libraries
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

library(regress)
library(gwaR)
library(limma)
library(edgeR)
library(parallel)
library(qvalue)
library(corrplot)

#' ## Load data

load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11/G_Z_matrix_Ss11.Rdata")
ls()

#' Rename R object to differentiate between mRNA and microRNA
Mdge <- dge
Mv <- v
Mwcen <- wcen
MG <- G

#' Remove R objects that will not be used in this analysis
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168", "MSUPRP")))
ls()

#' Load microRNA targets
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/5_filtered_targets_exp_rst.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/13_target_mrna_coloc_pqtl.Rdata")

#' Load microRNA expression data
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()

#' Rename R object to differentiate microRNA data from mRNA data
dge.mi <- dge
v.mi <- v
wcen.mi <- wtcen
G.mi <- G

#' Retain only the objects needed for the analysis
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168", "MSUPRP",
	"dge.mi","v.mi","wcen.mi","G.mi","MSUPRP_miRNA", "targets.exp", "coloc", "negcoloc")))

load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/12_mrna_mirna_resid_exp.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/9_mireqtl_pqtl_coloc_peaks.Rdata")
ls()

#' ## Analysis
#' 
#' Extract the phenotype object and subset the 166 animals:
pheno166<-MSUPRP$pheno[,,1]
pheno166<-pheno166[colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)],]
dim(pheno166)

#' Subset the phenotypes colocalized:
nm.pheno166<-as.character(unique(negcoloc$pheno))

colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)]

rfit.misub<-rfit.mi[colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)],]
dim(rfit.misub)

#' Build list of miRNAs and associated pQTL of interest (either miR-eQTL colocalizes, or target genes colocalize):
mirqtl<-c("juiciness", "driploss", "ph_24h", "protein", "cook_yield", "WBS", "tenderness", "overtend")
pheno.mir.list<-list("last_lum", "num_ribs", mirqtl, mirqtl, mirqtl, mirqtl, mirqtl, nm.pheno166)
names(pheno.mir.list)<-c("ssc-miR-190b","ssc-miR-184","ssc-miR-345-3p","ssc-let-7d-5p","ssc-let-7g","ssc-miR-95","ssc-miR-1468","ssc-miR-874")
rm(mirqtl)

#' Build list of mRNAs and associated pQTL of interest (negatively correlated target genes that co-localize with pQTL phenotypes)
#' 
#' Target genes co-localized with pQTL:
unique(negcoloc$ID)

#' Extract these genes from the rMfit object:
sub.rMfit<-rMfit[,colnames(rMfit)%in%unique(negcoloc$ID)]
dim(sub.rMfit)

#' Subset the miR-874 target genes:
pheno.mir874<-pheno.mir.list$`ssc-miR-874`

#' Function to perform correlation analysis of mRNA and miRNA espression and export a summary of correlation analysis
cor.exp <- function(target, transc, pheno, ...){
	x <- cor.test(transc, pheno, ...)
	rst <- data.frame(cor=x$estimate, 
		z=x$statistic, 
		pvalue=x$p.value)
	rownames(rst) <- target
	return(rst)
}

#' ---
#' Run corrlation analysis
#' 
#' Correlation between residual expression of miRNA and phenotypes of interest:
rst.cor<-lapply(names(pheno.mir.list), function(x) do.call(rbind, lapply(pheno.mir.list[[x]], 
	function(y) cor.exp(target=y, transc=rfit.misub[,x], pheno=pheno166[,y], alternative="two.sided", method="pearson"))))
names(rst.cor)<-names(pheno.mir.list)
rst.cor
rst.corfil<-lapply(rst.cor, function(x) x[x$pvalue<0.05,])
rst.corfil

#' Correlation between residual expression of mRNA and phenotypes of interest:
mrst.cor<-lapply(colnames(sub.rMfit), function(x) do.call(rbind, lapply(pheno.mir874, 
	function(y) cor.exp(target=y, transc=sub.rMfit[,x], pheno=pheno166[,y], alternative="two.sided", method="pearson"))))
names(mrst.cor)<-colnames(sub.rMfit)
mrst.cor

mrst.corfil<-lapply(mrst.cor, function(x) x[x$pvalue<0.05,])
names(mrst.corfil)<-paste(names(mrst.corfil),as.character(negcoloc[match(names(mrst.corfil), negcoloc$ID),"genes"]), sep=".")
mrst.corfil

#' Summarize results; split negcoloc by XLOC ID in mrst.cor:
splt.negcoloc<-lapply(names(mrst.cor), function(x) negcoloc[negcoloc$ID==x,])
names(splt.negcoloc)<-names(mrst.cor)
splt.negcoloc

#' ## Save data
save(rst.cor, rst.corfil, mrst.cor, mrst.corfil, splt.negcoloc, file="../16_mirna_mrna_pheno_corr.R")