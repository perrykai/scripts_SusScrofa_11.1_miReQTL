#' **Script:** `3_filter_target_genes_exp.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  12/13/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' 2. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom`
#' 
#' **Input File(s):** 
#' 
#' 1. `4_filtered_targetscan_rst.Rdata`
#' 
#' 2. `voom.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Output File(s):** 
#' 
#' 1. `5_filtered_targets_exp_rst.Rdata`
#' 
#' 2. `6_DAVID_bkgd_gene_names_exp.txt`
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
#' The objective of this script is to filter the target prediction results based on genes expressed in the 168 LD samples DV obtained in her analysis.
#' 
#' The filtered list will dictate which genes will be included in the correlation analysis and pQTL co-localization.
#' 
#' ## Install libraries
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

rm(list=ls())

library(limma)
library(edgeR)
library(methods)
library(biomaRt)

#' ## Load data
#' 
#' Load the list of filtered target prediction results:
load("../4_filtered_targetscan_rst.Rdata")
#' Load DV's dge object to obtain her annotation file:
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")

#' ## Analysis
#' 
genes<-dge$genes
genes$geneXID<-as.character(rownames(genes))
dim(genes)
head(genes)

#' How many genes have no name?
sum(is.na(genes$genes))

#' How many genes have ensembl IDs, but no gene symbol?
length(grep("ENSSSCG", as.character(genes$genes)))

#' Filter out the NA genes and the ensembl IDs:
genes<-genes[!is.na(genes$genes),]
dim(genes)
genes<-genes[!grepl("ENSSSCG", as.character(genes$genes)),]

#' Dimensions of retained genes object:
dim(genes)
head(genes)
#' How many unique genes:
length(unique(genes$genes))

tst<-targets.unique$"miR-128-3p"
head(tst)
sum(tst$external_gene_name %in% genes$genes)
dim(tst[tst$external_gene_name %in% genes$genes,])


#' Create a list of miRNA names to loop through:
miRnm<-as.character(names(targets.unique))
miRnm
#' Filter the genes targets based on those expressed in the dataset:
targets.exp<-list()
targets.exp.sum<-list()
for(i in miRnm){
targets.exp[[i]]<-targets.unique[[i]][targets.unique[[i]]$external_gene_name %in% genes$genes,]
# Create summary file:
targets.exp.sum[[i]]$summary<-data.frame(
	# How many targets were input into biomart? 
	gene.input=nrow(targets.unique[[i]]),
	# How many targets were expressed?
	gene.output=nrow(targets.exp[[i]]),
	# What fraction of the targets input were expressed/retained?
	gene.prop=nrow(targets.exp[[i]])/nrow(targets.unique[[i]]))
}

str(targets.exp)
targets.exp.sum

#' ---
#' 
#' Extract the gene symbols from the full expressed dataset to convert them to their ensembl gene IDs for use in DAVID:
#' 
length(unique(genes$genes))

#' ## Visualize

#' ## Save data
save(targets.exp, targets.exp.sum, file="../5_filtered_targets_exp_rst.Rdata")
write.table(unique(genes$genes), file="../6_DAVID_bkgd_gene_names_exp.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)