#!/bin/sh  -login
#PBS -l nodes=1:ppn=1,walltime=00:05:00,mem=1Gb
#PBS -N 3_novel_miRNA_precursor_blast_hsa
#PBS -j oe
#PBS -o /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/OutputsErrors
#PBS -m a
#PBS -M perrykai@msu.edu

#' **Script:** `3_novel_miRNA_precursor_blast_hsa.sh`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`
#' 
#' **Date:**  11/29/17
#' 
#' **Input File Directory:**  
#'
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/reference_sequences`
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`
#' 
#' **Input File(s):** 
#' 
#' 1. `6_novel_mature_mir.fa`
#' 2. `7_novel_precursor_mir.fa`
#' 3. `hsa_mature_mir.fa`
#' 4. `hsa_hairpin_mir.fa`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_novel_mirna_precursor_blast_hsa_e5.txt`
#' 2. `1_novel_mirna_mature_blast_hsa_e5.txt`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Analysis](#analysis)
#' 
#' ## Objectives
#' 
#' The objective of this script is to BLAST the candidate novel precursor miRNA sequences against the known human miRNA precursor sequences.
#' The table output format will be utilized for filtering the output in R.
#' 
#' ## Install libraries
#' 
module load BLAST+/2.2.30

#' ## Analysis

cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/reference_sequences/

#'
#' BLASTing the candidate novel miRNA precursor sequences against the human miRBase precursor databases:
#'
#' Parameters used: (See websites http://www.ncbi.nlm.nih.gov/BLAST/Why.shtml, http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ, and documentation http://www.ncbi.nlm.nih.gov/books/NBK279690/ for details)
#'
#' "blastn" option used for nucleotide BLAST:
#'
#' word-size: 7 (length of initial exact match);
#' reward: 1 (reward for a nt match);
#' penalty: -3 (penalty for a nt mismatch);
#'
#' "query" defines the input dataset (in this case, the small RNA seq data)
#'
#' "task" option defines which type of blastn to run; using blastn-short option which is optimized for sequences shorter than 50 bases
#'
#' "db" defines the databases used for the BLAST query: multiple databases can be utilized by including space-separated database names in quotes
#'
#' "out" defines the name of the output file
#'
#' "evalue" defines the 'expect value' used in the analysis, and can be used as a measure for the quality of the match between query and database.
#' The "expect value" is a parameter describing the number of hits one could expect "by chance" when searching a database of a particular size. Essentially, it describes the random background noise.
#'
#' "For example, an E value of 1 assigned to a hit can be interpreted as meaning that in a database of the current size one might expect to see 1 match with a similar score simply by chance."
#' So, for this analysis, the e-value threshold is set to 0.00001, meaning we can expect that a hit with that e-value would have 0.00001 match with a similar score by chance.
#' The lower the e-value, the more "significant" the match is.
#'
#' Tabular results:

blastn -query /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/7_novel_precursor_mir.fa -task blastn-short -db /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/reference_sequences/"hsa_hairpin_mir.fa" -out /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/1_novel_mirna_precursor_blast_hsa_e5.txt -evalue 0.00001 -outfmt 6

blastn -query /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/6_novel_mature_mir.fa -task blastn-short -db /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/reference_sequences/"hsa_mature_mir.fa" -out /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/1_novel_mirna_mature_blast_hsa_e5.txt -evalue 0.00001 -outfmt 6

qstat -f ${PBS_JOBID}