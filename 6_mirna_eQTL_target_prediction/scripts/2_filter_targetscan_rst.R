#' **Script:** `2_filter_targetscan_rest.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  `12/13/17`
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Input File(s):** `3_targetscan_mirna_eqtl_rst.txt`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Output File(s):** `4_filtered_targetscan_rst.Rdata`
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
#' The objective of this script is to filter the miRNA target prediction results 
#' based on target site type, retaining 8mer-1a and 7mer-m8, the two most preferentially conserved and effective classes of target sites.  
#' The results will then be filtered to one miRNA-gene target pair to create a list of putative target genes for correlation analysis and co-localization analysis with pQTL.
#' 
#' Will need to translate transcript IDs back to gene names once filtering is complete. 
#' Also, final data will be organized as nested list by miRNA.
#' 
#' THIS ANALYSIS CONDUCTED IN R/3.1.0
#' 
#' ## Install libraries
library(biomaRt)

#' ## Load data
rm(list=ls())
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

tprst<-read.table("../3_targetscan_mirna_eqtl_rst.txt", 
	sep="\t",
	header=TRUE)
head(tprst)
dim(tprst)
#' ## Analysis
#' 
#' First, filter the target prediction results based on site type: retain 8mer and 7mer-m8:
table(tprst$Site_type)
tp8m<-tprst[tprst$Site_type==c("8mer-1a", "7mer-m8"),c(1,2,9)]
tp8m$miRNA_family_ID<-as.character(tp8m$miRNA_family_ID)
tp8m$Site_type<-as.character(tp8m$Site_type)
#' 
#' The ensembl id's have a decimal at the end, indicating which version the ID comes from. In this case, the vesion number will be removed.
#' 
#' Remove the version numbers (".1", ".2" or ".3") after the ensembl ids to another column for each gene to make it compatible with biomaRt
#' 
#' See website http://useast.ensembl.org/Help/View?id=181 for details on the stability of ensembl IDs
tp8m$ensid_version <- as.character(lapply(strsplit(as.character(tp8m$a_Gene_ID), '.', fixed = TRUE), "[", 2))
tp8m$a_Gene_ID<-as.character(gsub("\\..*","",tp8m$a_Gene_ID))

head(tp8m)
table(tp8m$ensid_version)

#' How many putative targets are retained?
dim(tp8m)

#' What fraction of results remain?
nrow(tp8m)/nrow(tprst)
#' Excellent, much easier to work with!
#' 
#' Do I still have putative targets for all my miRNA families?
length(unique(tp8m$miRNA_family_ID))
miRnm<-unique(as.character(tp8m$miRNA_family_ID))

#' How many of each miRNA do I have targets for?
for (i in miRnm){
	print(sum(tp8m$miRNA_family_ID == i))
}
#' If our concern in a gene list of targets for each miRNA, all that needs to be done is to create a list 
#' containing the transcript IDs for each putative target.

targets.unique<-list()
for(i in miRnm){
	targets.unique[[i]]<-tp8m[tp8m$miRNA_family_ID==i,]
	targets.unique[[i]]<-subset(targets.unique[[i]], !duplicated(a_Gene_ID))
}
str(targets.unique)

#' ---
#' 
#' Convert the ensembl ID to human gene ID:
#' 
#' Use biomaRt package (v/2.20.0) to obtain annotation information for Ensembl dataset BLAST hits:
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="jul2016.archive.ensembl.org", dataset="hsapiens_gene_ensembl")

filters<-listFilters(mart)
dim(filters)
filters[60:70,]

attributes<-listAttributes(mart)
dim(attributes)

targets.rst <- list()
system.time(
for (i in miRnm){
targets.rst[[i]]<- getBM(attributes=c("ensembl_transcript_id", 
		"ensembl_gene_id", 
		"external_gene_name",
		"transcript_status", 
		"transcript_version", 
		"ens_hs_transcript",
		"hgnc_symbol",
		"hgnc_transcript_name",
		"gene_biotype"),
		filters="ensembl_transcript_id", 
		values=targets.unique[[i]]$a_Gene_ID, 
		mart=mart)
# Use match function to obtain the gene id, transcript id, gene name, hgnc gene symbol(check for equality), and gene status for each ensembl match and add it as a column to targets.unique:
targets.unique[[i]]$ensembl_gene_id<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "ensembl_gene_id"])
targets.unique[[i]]$ensembl_transc_id<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "ensembl_transcript_id"])
targets.unique[[i]]$external_gene_name<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "external_gene_name"])
targets.unique[[i]]$hgnc_symbol<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "hgnc_symbol"])
targets.unique[[i]]$status<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "transcript_status"])
})

str(targets.unique)

targets.rst.sum<-list()
for (i in miRnm){
targets.rst.sum[[i]]$summary<-data.frame(
	# How many transcript IDs were input into biomart? 
	transc.id.input=nrow(targets.unique[[i]]),
	# How many transcript IDs are NA (possibly retired)?
	transc.id.na=sum(is.na(targets.unique[[i]]$ensembl_transc_id)),
	# What fraction of the input transcript IDs were NA?
	transc.id.na.prop=sum(is.na(targets.unique[[i]]$ensembl_transc_id)) / nrow(targets.unique[[i]]),
	# Of the transcripts remaining (not NA), check for any differences between the input transcript ID and the output ensembl_transc_id:
	transc.id.diff=sum(targets.unique[[i]]$a_Gene_ID != targets.unique[[i]]$ensembl_transc_id, na.rm=TRUE),
	# Do any of the gene names differ between the "external gene name" and the "HGNC symbol" columns?
	gene.name.diff=sum(targets.unique[[i]]$external_gene_name != targets.unique[[i]]$hgnc_symbol, na.rm=TRUE))
}

targets.rst.sum

#' ## Visualize


#' ## Save data
#' 
#' Save the list of unique targets for each miRNA, and the summary file 
save(targets.unique, targets.rst, targets.rst.sum, file="../4_filtered_targetscan_rst.Rdata")