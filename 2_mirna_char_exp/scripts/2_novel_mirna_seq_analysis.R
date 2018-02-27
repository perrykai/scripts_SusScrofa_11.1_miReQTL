#' **Script:** `2_novel_mirna_seq_analysis.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`
#' 
#' **Date:**  11/29/17
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`
#' 
#' **Input File(s):** `2_predicted_novel_mirna.csv`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`
#' 
#' **Output File(s):** 
#' 
#' 1. `5_novel_miRNA_filtered.Rdata`
#' 2. `6_novel_mature_mir.fa`
#' 3. `7_novel_precursor_mir.fa`
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
#' The objective of this script is to filter the extracted putative novel miRNAs for the following characteristics:
#'  
#' 1. Novel miRNA candidate must have > 90% probability of being a true positive (filter by miRDeep2 score and then estimated probability that the miRNA candidate is a true positive)
#' 2. Hairpins must have significant Randfold p-values
#' 3. Minimum read counts for putative mature and star strand sequences
#' 
#' Additionally, the novel candidate miRNAs will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)
#'
#' ## Install libraries
library("ShortRead")

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts")
#' ## Load data
#' 
novelmir<-read.table("../2_predicted_novel_mirna.csv", sep="\t", header=TRUE)
dim(novelmir)
colnames(novelmir)
head(novelmir)

#' ## Analysis
#' 

#' ### 1. Filter by the miRDeep2 score and the estimated probability of the novel miRNA being true positives:
sum(novelmir$miRDeep2.score >= 10)
table(novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)
novelmir[novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive == "91 +/- 1%", "miRDeep2.score"]

#' First filter by miRDeep2 score, then by estimated probability:
novelmir10<-novelmir[novelmir$miRDeep2.score >= 10,]
dim(novelmir10)

table(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)

sum(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive != "91 +/- 2%")
#' Filtering by miRDeep2 score removed any miRNA with estimated probability < "91 +/- 2%"
#' 

#' ### 2. Hairpins must have a significant Randfold p-value
sum(novelmir10$significant.randfold.p.value == "yes")

#' 271 potential miRNA candidates have significant Randfold p-value
novelmir10sigrandfold<-novelmir10[novelmir10$significant.randfold.p.value == "yes",]
dim(novelmir10sigrandfold)
head(novelmir10sigrandfold)

#' Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences
table(novelmir10sigrandfold$rfam.alert)
sum(novelmir10sigrandfold$rfam.alert != "-")


#' ### 3. Minimum read counts for putative mature and star strand sequences: require at least 10 counts for each
summary(novelmir10sigrandfold$mature.read.count)
summary(novelmir10sigrandfold$star.read.count)

sum(novelmir10sigrandfold$mature.read.count <= 10)

novelmir10sigrandfoldmincounts<-novelmir10sigrandfold[novelmir10sigrandfold$mature.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)

sum(novelmir10sigrandfoldmincounts$star.read.count <= 10)
novelmir10sigrandfoldmincounts<-novelmir10sigrandfoldmincounts[novelmir10sigrandfoldmincounts$star.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)

sum(novelmir10sigrandfoldmincounts$star.read.count <=10)
sum(novelmir10sigrandfoldmincounts$mature.read.count <=10)


#' Investigate the final output:
dim(novelmir10sigrandfoldmincounts)
head(novelmir10sigrandfoldmincounts)

#' How many of the novel miRNA candidates have a human miRBase miRNA with the same seed sequence
sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed != "-")

sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed == "-")

#' What is the summary of the mature and star read counts now?
summary(novelmir10sigrandfoldmincounts$mature.read.count)
summary(novelmir10sigrandfoldmincounts$star.read.count)

#' Subset the provisional.id, consensus.mature.sequence and consensus.precursor.sequence for use in BLASTN against known human and mouse miRBase sequences:
novelmircandidateBLAST<-novelmir10sigrandfoldmincounts[,c("provisional.id", "consensus.mature.sequence", "consensus.precursor.sequence")]
dim(novelmircandidateBLAST)
head(novelmircandidateBLAST)

#' ### 4. The novel candidate miRNAs (both precursor and mature) will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)
#' 
#' #### First, prepare the sequence names for the candidate novel mature miR
seqnamesmature<-paste("seq|",novelmircandidateBLAST$provisional.id, "|candidate novel mature miR")
seqnamesmature<-gsub(" ", "", seqnamesmature, fixed = TRUE)
head(seqnamesmature)

#' Create the BStringSet object with the mature sequence names
matureids<-BStringSet(seqnamesmature)
head(matureids)

#' Prepare the candidate novel mature sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet
novelmircandidateBLAST$consensus.mature.seq.tadjust <- gsub("u","t", novelmircandidateBLAST$consensus.mature.sequence)

#' Create the DNAStringSet object with the candidate novel mature sequence reads
novelmatureseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.mature.seq.tadjust)
head(novelmatureseqreads)

#' Use the ShortRead command to combine the candidate novel mature sequence IDs and the sequences
candidatenovelmaturefasta<-ShortRead(sread=novelmatureseqreads ,id=matureids)


#' #### Then, prepare the sequence names for the candidate novel precursor miR:
seqnamesprecursor<-paste("seq|", novelmircandidateBLAST$provisional.id, "|candidate novel percursor miR")
seqnamesprecursor<-gsub(" ", "", seqnamesprecursor, fixed = TRUE)

head(seqnamesprecursor)

#' Create the BStringSet object with the precursor sequence names
precursorids<-BStringSet(seqnamesprecursor)
head(precursorids)

#' Prepare the candidate precursor sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet
novelmircandidateBLAST$consensus.precursor.seq.tadjust<- gsub("u","t", novelmircandidateBLAST$consensus.precursor.sequence)

#' Create the DNAStringSet object with the novel candidate precursor sequence reads
novelprecursorseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.precursor.seq.tadjust)
head(novelprecursorseqreads)

#' Use the ShortRead command to combine the candidate novel precursor sequence IDs and the sequences
candidatenovelprecursorfasta<-ShortRead(sread=novelprecursorseqreads ,id=precursorids)

#' ## Visualize
#' 

#' ## Save data
#' 
save(novelmir10sigrandfoldmincounts, file="../5_novel_miRNA_filtered.Rdata")
writeFasta(candidatenovelmaturefasta, file="../6_novel_mature_mir.fa")
writeFasta(candidatenovelprecursorfasta, file="../7_novel_precursor_mir.fa")