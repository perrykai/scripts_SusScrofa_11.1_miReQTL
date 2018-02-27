#!/bin/sh  -login
#PBS -l nodes=1:ppn=1,walltime=00:05:00,mem=1Gb
#PBS -N 5_copy_novel_precursor_pdfs_e5
#PBS -j oe
#PBS -o /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/OutputsErrors
#PBS -m a
#PBS -M perrykai@msu.edu

#' **Script:** `5_copy_novel_precursor_pdfs_e5.sh`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`
#' 
#' **Date:**  11/30/17
#' 
#' **Input File Directory:**  
#'
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output`
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/2_quant_predict_output/pdfs_27_11_2017_t_18_08_27`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_filtered_novel_mirna_precursor_ids_e5.txt`
#' 2. One predicted secondary structure pdf file for each of the 27 unique candidate novel precursors. 
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/`
#' 
#' **Output File(s):** 
#' 
#' 1. The predicted secondary structure pdf files for each of the 27 unique candidate novel precursors will be copied to the output directory from pdfs_27_11_2017_t_18_08_27
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Analysis](#analysis)
#' 3. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' The objective of this script is to copy the predicted secondary structure pdf files for each of the 27 unique candidate novel precursors in order to build a figure containing graphics of the secondary structures of potential novel miRNAs.
#' 
#' ## Analysis
#' 
#' 
#' ## Save data

cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output

cut -d '"' -f2 3_filtered_novel_mirna_precursor_ids_e5.txt > 4_cut_filtered_ids.txt
sed -e 's/$/.pdf/' -i 4_cut_filtered_ids.txt

f1=`cat /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/4_cut_filtered_ids.txt`
cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/2_quant_predict_output/pdfs_27_11_2017_t_18_08_27

for file in `cat /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/4_cut_filtered_ids.txt`;
do cp "$file" /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/2_novel_mirna_precursor_pdf_e5/;
done

qstat -f ${PBS_JOBID}