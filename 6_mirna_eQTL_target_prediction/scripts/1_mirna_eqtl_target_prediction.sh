#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel16,walltime=00:20:00,mem=4Gb
#PBS -N miRNA_eQTL_targetscan
#PBS -j oe
#PBS -o /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/OutputsErrors/
#PBS -m a
#PBS -M perrykai@msu.edu

#' **Script:** `1_mirna_eqtl_target_prediction.sh`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  `12/12/2017`
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Input File(s):** 
#' 
#' 1. `miR_Family_Info.txt.zip`
#' 
#' 2. `UTR_Sequences.txt.zip`
#' 
#' **Output File Directory:** `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/6_miRNA_eQT_target_prediction`
#' 
#' **Output File(s):** 
#' 
#' 1. `3_targetscan_miRNA_eQTL_output.txt`
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
#' The objective of this analysis is to conduct miRNA target prediction for the 17 miRNA with significant eQTL.
#' This will be conducted using TargetScan_70 and human miRNA-mRNA interactions. 
#' The 17 miRNA seed sequences will be compared between human and pig to ensure proper interactions.
#' The results of target prediction will be filtered to retain only those genes expressed in the LD of the MSUPRP pigs.
#' That information will be obtained from DV. 
#' 
#' 1. Download the required files from the TargetScan database.
#' 
#' 2. Filter the miRNA file for the 17 miRNAs in my dataset.
#' 
#' 3. Run the first script in the TargetScan pipeline, targetscan_70.pl
#' 
#' ## Load data
cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/

#' Download the miRNA family info database from targetscan.org:
# wget http://www.targetscan.org/vert_71/vert_71_data_download/miR_Family_Info.txt.zip

# --2017-12-12 14:39:39--  http://www.targetscan.org/vert_71/vert_71_data_download/miR_Family_Info.txt.zip
# Resolving www.targetscan.org... 18.4.1.128
# Connecting to www.targetscan.org|18.4.1.128|:80... connected.
# HTTP request sent, awaiting response... 200 OK
# Length: 178768 (175K) [application/zip]
# Saving to: “miR_Family_Info.txt.zip”

# 100%[=====================================================================================>] 178,768     1013K/s   in 0.2s

# 2017-12-12 14:39:40 (1013 KB/s) - “miR_Family_Info.txt.zip” saved [178768/178768]


#' Unzip the miR_Family file:
unzip miR_Family_Info.txt.zip


#' Download the 3'UTR database from targetscan.org:
# wget http://www.targetscan.org/vert_71/vert_71_data_download/UTR_Sequences.txt.zip
# --2017-12-12 14:40:23--  http://www.targetscan.org/vert_71/vert_71_data_download/UTR_Sequences.txt.zip
# Resolving www.targetscan.org... 18.4.1.128
# Connecting to www.targetscan.org|18.4.1.128|:80... connected.
# HTTP request sent, awaiting response... 200 OK
# Length: 795617238 (759M) [application/zip]
# Saving to: “UTR_Sequences.txt.zip”

# 100%[=====================================================================================>] 795,617,238 3.13M/s   in 1m 57s

# 2017-12-12 14:42:22 (6.47 MB/s) - “UTR_Sequences.txt.zip” saved [795617238/795617238]

#' Unzip the UTR sequences:
unzip UTR_Sequences.txt.zip

#' ## Analysis
#' 
#' This line puts the UTR file into the correct format for running the analysis:
sed '1,1d' UTR_Sequences.txt | cut -f1,4,5 > UTR_sequences_all.txt

#' Then, grep the human species ID, which is 9606 (http://www.targetscan.org/vert_61/docs/species.html)
#' out of the file to extract the human 3' UTRs.
grep -w "9606" UTR_sequences_all.txt > UTR_sequences_human.txt

#' How many lines are in this file?
wc -l UTR_sequences_human.txt
# 28352 UTR_sequences_human.txt

#' Check that the sequence ID is 9606 for all of these lines:
cut -f2 UTR_sequences_human.txt |uniq
# 9606

#' Remove the (giant) original 3' UTR files:
rm UTR_Sequences.txt
rm UTR_sequences_all.txt

#' ------------

#' This line puts the miRNA family info file into the correct format for the analysis:
sed '1,1d' miR_Family_Info.txt| cut -f1,2,3 > miR_Family_Info_all.txt

#' Then, grep the human species ID, which is 9606 (http://www.targetscan.org/vert_61/docs/species.html)
#' out of the file to extract the human miRNAs.
grep -w "9606" miR_Family_Info_all.txt > miR_Family_Info_Human.txt

#' How many lines are in this file?
wc -l miR_Family_Info_Human.txt
# 2603 miR_Family_Info_Human.txt

cut -f3 miR_Family_Info_Human.txt | uniq
# 9606

#' -------------
#' Filter the miRNA file for the families contained in the eQTL results:
#' 
#' Manually created the file 1_hsa_mirna_names.txt based on seed sequence matching between pig and human miRNA sequences.
#' 
#' Notes:
#' 
#' ssc-let-7d-5p and ssc-let-7g have the same seed seq; listed as let-7-5p/98-5p in human dataset.
#' 
#' ssc-miR-128 listed as hsa-miR-128-3p
#' 
#' ssc-miR-140-5p listed as hsa-miR-140-5p, but seed sequence has 1 nt shift in human compared to ssc.
#' 
#' ssc-miR-1468 listed as hsa-miR-1468-5p
#' 
#' ssc-miR-190b listed as hsa-miR-190-5p
#' 
#' ssc-miR-345-3p listed as hsa-miR-345-3p, but seed sequence has 1 nt shift in human compared to ssc.
#' 
#' ssc-miR-429 listed as hsa-miR-200-3p/429
#' 
#' ssc-miR-6782-3p listed as hsa-miR-6821-3p (matching seed seq to ssc miR)
#' 
#' ssc-miR-7135-3p listed as hsa-miR-6888-3p (matching seed seq to ssc miR)
#' 
#' ssc-miR-874 listed as hsa-miR-874-3p
#' 
#' ssc-miR-95 listed as hsa-miR-95-3p
#' 
#' ssc-miR-9785-5p listed as hsa-miR-6072/6891-3p; seed sequence has 1 nt shift in human compared to ssc.
#' 
#' ssc-miR-9810-3p and ssc-miR-9843-3p have no matching seed sequences in the human miR database from targetscan and are therefore not included in this analysis.

cat 1_hsa_mirna_names.txt
# let-7-5p/98-5p
# miR-128-3p
# miR-1306-3p
# miR-140-5p
# miR-1468-5p
# miR-184
# miR-190-5p
# miR-345-3p
# miR-200-3p/429
# miR-6821-3p
# miR-6888-3p
# miR-874-3p
# miR-95-3p
# miR-6072/6891-3p

grep -f 1_hsa_mirna_names.txt miR_Family_Info_Human.txt | uniq > 2_hsa_mirna_targetscan_input.txt

cat 2_hsa_mirna_targetscan_input.txt
# miR-6888-3p	UCUGUCU	9606
# miR-6821-3p	GACCUCU	9606
# miR-6072/6891-3p	CCUCAUC	9606
# miR-1468-5p	UCCGUUU	9606
# miR-1306-3p	CGUUGGC	9606
# miR-874-3p	UGCCCUG	9606
# miR-345-3p	CCCUGAA	9606
# miR-200-3p/429	AAUACUG	9606
# miR-190-5p	GAUAUGU	9606
# miR-184	GGACGGA	9606
# miR-140-5p	AGUGGUU	9606
# miR-128-3p	CACAGUG	9606
# miR-95-3p	UCAACGG	9606
# let-7-5p/98-5p	GAGGUAG	9606

#' This file is complete for use in target prediction analysis.
#' 
#' Remove the original miRNA info files:
rm miR_Family_Info.txt
rm miR_Family_Info_all.txt

#' -------------
cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts

/mnt/home/perrykai/targetscan_70/targetscan_70.pl ../2_hsa_mirna_targetscan_input.txt ../UTR_sequences_human.txt ../3_targetscan_mirna_eqtl_rst.txt

qstat -f ${PBS_JOBID}