#' **Script:** `10_peak_SNP_geno_exp.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  `1/29/18`
#' 
#' **Input File Directory:**  
#' 
#' 1. `../../4_dge_G_objects/`
#' 
#' 2. `../../6_mirna_eQTL_target_prediction/`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`
#' 
#' 2. `12_mrna_mirna_resid_exp.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts/10_peak_SNP_geno_exp_literate/figure/`
#' 
#' **Output File(s):** `geno_exp_mir874.tiff`
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
#' ## Install libraries
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)
library(ggplot2)

#' ## Load data
#' 
#' Load DV's eqtl functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

#' Load required data:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
load("../../6_mirna_eQTL_target_prediction/12_mrna_mirna_resid_exp.Rdata")
ls()

#' ## Analysis

#' ---
#' 
#' Following up on pQTL SNP fixed for 7 phenotypes on SSC15: how highly correlated are MARC0093624 and H3GA0052416?
cor.test(Z[,"MARC0093624"],Z[,"H3GA0052416"])

#' ---
#' 
#' Investigate peak SNP for miR-874: ALGA0016550. Do we see differences in miRNA expression segragating between genotypes?
#' 
#' How many animals have each genotype?
MSUPRP_miRNA$map["ALGA0016550",]

table(MSUPRP_miRNA$geno[,"ALGA0016550"])

summary(v$E["ssc-miR-874",])

sum(names(MSUPRP_miRNA$geno[,"ALGA0016550"])!=names(v$E["ssc-miR-874",]))

mir874exp<-data.frame("geno"=as.character(MSUPRP_miRNA$geno[,"ALGA0016550"]), "exp"=v$E["ssc-miR-874",], check.rows=TRUE, stringsAsFactors=FALSE)

median(mir874exp[mir874exp$geno=="0","exp"])
median(mir874exp[mir874exp$geno=="1","exp"])
median(mir874exp[mir874exp$geno=="2","exp"])

mir874exp$geno<-gsub("0", "AA", mir874exp$geno)
mir874exp$geno<-gsub("1", "AB", mir874exp$geno)
mir874exp$geno<-gsub("2", "BB", mir874exp$geno)
head(mir874exp)

#' ---
#' 
#' Investigate peak SNP for miR-429: ALGA0118516. Do we see differences in miRNA expression segragating between genotypes?
#' 
#' How many animals have each genotype?
MSUPRP_miRNA$map["ALGA0118516",]

table(MSUPRP_miRNA$geno[,"ALGA0118516"])

summary(v$E["ssc-miR-429",])

sum(names(MSUPRP_miRNA$geno[,"ALGA0118516"])!=names(v$E["ssc-miR-429",]))

mir429exp<-data.frame("geno"=as.character(MSUPRP_miRNA$geno[,"ALGA0118516"]), "exp"=v$E["ssc-miR-429",], check.rows=TRUE, stringsAsFactors=FALSE)

median(mir429exp[mir429exp$geno=="0","exp"])
median(mir429exp[mir429exp$geno=="1","exp"])
median(mir429exp[mir429exp$geno=="2","exp"])

mir429exp$geno<-gsub("0", "AA", mir429exp$geno)
mir429exp$geno<-gsub("1", "AB", mir429exp$geno)
mir429exp$geno<-gsub("2", "BB", mir429exp$geno)
head(mir429exp)


#' ## Visualize
#+ geno_exp_mir874, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
mir874genos<-ggplot(data=mir874exp, aes(x=geno,y=exp, col=geno)) +
geom_jitter(width=0.25) +
geom_segment(aes(x=0.75, y=median(mir874exp[mir874exp$geno=="AA","exp"]), xend=1.25, yend=median(mir874exp[mir874exp$geno=="AA","exp"])), color="black") +
geom_segment(aes(x=1.75, y=median(mir874exp[mir874exp$geno=="AB","exp"]), xend=2.25, yend=median(mir874exp[mir874exp$geno=="AB","exp"])), color="black") +
geom_segment(aes(x=2.75, y=median(mir874exp[mir874exp$geno=="BB","exp"]), xend=3.25, yend=median(mir874exp[mir874exp$geno=="BB","exp"])), color="black") +
labs(x="Genotype", y="Normalized miR-874 Expression", title="ALGA0016550 [T/C]", colour="Genotype") +
scale_color_manual(values=c("purple","orange","darkgreen"))+
theme_bw()
mir874genos

#+ geno_exp_mir429, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
mir429genos<-ggplot(data=mir429exp, aes(x=geno,y=exp, col=geno)) +
geom_jitter(width=0.25) +
geom_segment(aes(x=0.75, y=median(mir429exp[mir429exp$geno=="AA","exp"]), xend=1.25, yend=median(mir429exp[mir429exp$geno=="AA","exp"])), color="black") +
geom_segment(aes(x=1.75, y=median(mir429exp[mir429exp$geno=="AB","exp"]), xend=2.25, yend=median(mir429exp[mir429exp$geno=="AB","exp"])), color="black") +
geom_segment(aes(x=2.75, y=median(mir429exp[mir429exp$geno=="BB","exp"]), xend=3.25, yend=median(mir429exp[mir429exp$geno=="BB","exp"])), color="black") +
labs(x="Genotype", y="Normalized miR-429 Expression", title="ALGA0118516 [T/C]", colour="Genotype") +
scale_color_manual(values=c("purple","orange","darkgreen"))+
theme_bw()
mir429genos

#' ## Save data
