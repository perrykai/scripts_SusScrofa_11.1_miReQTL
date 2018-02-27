#' **Script:** `4_filter_novel_mirna_hsa_blast_results.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts/`
#' 
#' **Date:**  `11/29/17`
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_novel_mirna_precursor_blast_hsa_e5.txt`
#' 2. `1_novel_mirna_mature_blast_hsa_e5.txt`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/`
#' 
#' **Output File(s):** 
#' 
#' 1. `2_filtered_novel_mirna_precursor_abundance_e5.txt`
#' 2. `2_novel_mirna_mature_abundance_e5.txt`
#' 3. `3_filtered_novel_mirna_precursor_ids_e5.txt`
#' 4. `3_filtered_novel_mirna_mature_ids_e5.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to summarize the output from the BLAST query in order to characterize the candidate novel miRNAs present in the small RNA seq data.
#' 
#' So, need to load the precursor BLAST dataset and the full dataset of candidate novel miRNA and compare the names of sequences to see if the most abundant miRNA had BLAST results.
#' 
#' ## Install libraries
rm(list=ls())

#' ## Load data
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts")

hsa.blast.precursor<-read.table("../1_blast_novel_mirna_output/1_novel_mirna_precursor_blast_hsa_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
mature.hsa.blast<-read.table("../1_blast_novel_mirna_output/1_novel_mirna_mature_blast_hsa_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
#' --------------------
dim(hsa.blast.precursor)
head(hsa.blast.precursor)
str(hsa.blast.precursor)

hsa.blast.precursor$dbseq_id<-as.character(hsa.blast.precursor$dbseq_id)
hsa.blast.precursor$query_id<-as.character(hsa.blast.precursor$query_id)
#' --------------------
dim(mature.hsa.blast)
head(mature.hsa.blast)
str(mature.hsa.blast)

mature.hsa.blast$dbseq_id<-as.character(mature.hsa.blast$dbseq_id)
mature.hsa.blast$query_id<-as.character(mature.hsa.blast$query_id)
#' --------------------
load("../5_novel_miRNA_filtered.Rdata")

dim(novelmir10sigrandfoldmincounts)
names(novelmir10sigrandfoldmincounts)
#' Make the consensus.mature.sequence column into a character vector to count the length of the strings
novelmir10sigrandfoldmincounts$consensus.mature.sequence<-as.character(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
rownames(novelmir10sigrandfoldmincounts)<- novelmir10sigrandfoldmincounts$provisional.id
head(novelmir10sigrandfoldmincounts)
novelmir10sigrandfoldmincounts$consensus.mature.sequence.length<-nchar(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
head(novelmir10sigrandfoldmincounts)
#' ## Analysis
#' 
#' ### 1. The goal is to compare the candidate novel miRNAs with BLAST results to the most abundant candidate novel miRNAs.
#' 
#' First, return the name of the sequences to the way they were before BLASTing at e-value = 1x10^-5.

hsa.blast.precursor$seqname<-sapply(strsplit(hsa.blast.precursor$query_id, "|", fixed=TRUE),'[',2)
head(hsa.blast.precursor)

#' Identify the unique number of sequences with BLAST hits at eval = 1x10^-5
seqids<-unique(hsa.blast.precursor$seqname)
length(seqids)
seqids

#' Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset
match(seqids, rownames(novelmir10sigrandfoldmincounts))
rownames(novelmir10sigrandfoldmincounts)%in%seqids

#' Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.
novelmir.abundance<-novelmir10sigrandfoldmincounts[seqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count", "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]

#' Is this object ordered by miRDeep2 score or by total.read.count?
sum(rownames(novelmir.abundance[order(novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(novelmir.abundance))
sum(rownames(novelmir.abundance[order(novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(novelmir.abundance))

#' Combine the information; are the predicted miRNA precursors with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?
#' 
#' This object (candidate novel precursors with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.
novelmir.abundance<-cbind(novelmir.abundance[match(hsa.blast.precursor$seqname, rownames(novelmir.abundance)),], hsa.blast.precursor)
sum(novelmir.abundance$seqname != novelmir.abundance$provisional.id)
head(novelmir.abundance)

novelmir.abundance<-novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(novelmir.abundance)
head(novelmir.abundance)
novelmir.provisional.ids<-unique(novelmir.abundance$provisional.id)

#' ---------------------------------------
#' 
#' Now to repeat the same analysis on the candidate novel mature miRNAs with BLAST results at an e-value of 1x10-5
#' 
#' First, return the name of the mature sequences to the way they were before BLASTing at e-value = 1x10^-5.

mature.hsa.blast$seqname<-sapply(strsplit(mature.hsa.blast$query_id, "|", fixed=TRUE),'[',2)
head(mature.hsa.blast)

#' Identify the unique number of sequences with BLAST hits at eval = 1x10^-5
matureseqids<-unique(mature.hsa.blast$seqname)
length(matureseqids)
matureseqids

#' Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset
match(matureseqids, rownames(novelmir10sigrandfoldmincounts))
rownames(novelmir10sigrandfoldmincounts)%in%matureseqids

#' Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.
mature.novelmir.abundance<-novelmir10sigrandfoldmincounts[matureseqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]

#' Is this object ordered by miRDeep2 score or by total.read.count?
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))

#' Combine the information; are the predicted miRNAs with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?
#' 
#' This object (candidate novels with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.
mature.novelmir.abundance<-cbind(mature.novelmir.abundance[match(mature.hsa.blast$seqname, rownames(mature.novelmir.abundance)),], mature.hsa.blast)
sum(mature.novelmir.abundance$seqname != mature.novelmir.abundance$provisional.id)
head(mature.novelmir.abundance)

mature.novelmir.abundance<-mature.novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(mature.novelmir.abundance)
head(mature.novelmir.abundance)
mature.novelmir.provisional.ids<-unique(mature.novelmir.abundance$provisional.id)

#' ## Save data
write.table(novelmir.abundance, file="../1_blast_novel_mirna_output/2_filtered_novel_mirna_precursor_abundance_e5.txt")
write.table(mature.novelmir.abundance, file="../1_blast_novel_mirna_output/2_novel_mirna_mature_abundance_e5.txt")
#' Create a list of the pertinent precursor provisional.ids to extract the correct pdfs for a supplemental figure
write.table(novelmir.provisional.ids, file="../1_blast_novel_mirna_output/3_filtered_novel_mirna_precursor_ids_e5.txt", row.names=FALSE, col.names=FALSE)
write.table(mature.novelmir.provisional.ids, file="../1_blast_novel_mirna_output/3_filtered_novel_mirna_mature_ids_e5.txt", row.names=FALSE, col.names=FALSE)