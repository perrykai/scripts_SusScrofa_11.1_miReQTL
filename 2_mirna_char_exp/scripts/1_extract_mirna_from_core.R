#' **Script:** `1_extract_mirna_from_core.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`
#' 
#' **Date:**  11/29/17
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/2_quant_predict_output`
#' 
#' **Input File(s):**  `result_27_11_2017_t_18_08_27.csv`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_core_output_stats.csv`
#' 
#' 2. `2_predicted_novel_mirna.csv`
#' 
#' 3. `3_mature_mirna_detected.csv`
#' 
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
#' This code extracts the miRDeep2 core module output from the `result_27_11_2017_t_18_08_27.csv` file.
#' This .csv file contains the miRDeep2 score distribution for the novel miRNA prediction step, 
#' the novel miRNAs predicted by miRDeep2, the mature miRBase miRNAs detected by miRDeep2, and the 
#' miRBase miRNAs not detected by miRDeep2. The first three of these items are extracted from this .csv file using this script.
#' The objective here is to characterize the known and novel miRNAs present in this dataset, isolate the sequences meeting miRDeep2 score of 10 or greater,
#' to isolate the sequences at that score cutoff having homologous seed sequence with a human miRBase miRNA,
#' and to estimate the false discovery rate of the miRDeep2 prediction step at miRDeep2 score 10. 

#' 
#' ## Install libraries
#' 
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts/")

#' ## Load Data
#' 
#' To isolate the first section of the csv containing the miRDeep2 distribution scores:
sts<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", nrows=21, header = TRUE)

head(sts)
tail(sts)

#' To isolate the second section of the csv containing the novel predicted miRNAs:
novelmirna<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", skip=26, nrow=1201,  header = TRUE, fill=TRUE)

head(novelmirna)
tail(novelmirna)

#' To isolate the third section of the csv containing the miRBase miRNAs detected by miRDeep2:
md<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", skip = 1232, nrow = 289, header = TRUE, fill = TRUE)

head(md)
tail(md)

#' ## Analysis
#' 
#' Now I can open the `novel_mirna_predicted.csv` file and filter by various thresholds
#' 
#' ### 1. miRDeep2 score of 10 or more
colnames(novelmirna)

dim(novelmirna)	

score10novelmirna<-novelmirna[novelmirna$miRDeep2.score >= 10, ]

dim(score10novelmirna)

nrow(score10novelmirna)

#' So, there are 369 predicted miRNA precursors with a miRDeep2 score > or = 10
#' 
#' Estimated false positives is 34 +/- 6 (obtained from `1_extracted_mirdeep2_core_output_stats.csv`)

34 - 6

34 + 6

28/369

40/369


#' ### 2. Significant Randfold p-value
#' 
#' Now, subset this again into those that had a significant Randfold p value, indicating ability of secondary structure formation
head(score10novelmirna$significant.randfold.p.value)

sum(score10novelmirna$significant.randfold.p.value=="yes")
randfoldsigpval<-score10novelmirna[score10novelmirna$significant.randfold.p.value == "yes", ]
dim(randfoldsigpval)

#' Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences
rfam<-randfoldsigpval$rfam.alert
sum(rfam=="-")
sum(rfam !="-")


#' ### 3. Do the putative novel sequences have a homologous human miRNA seed sequence?
sum(randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-")

#' This indicates that 114 of the sequences have a homologous human seed sequence
homologseed<-randfoldsigpval[randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-",]

#' Subset of the 100 sequences
dim(homologseed)


write.table(sts, file="../1_core_output_stats.csv", sep = "\t", col.names = TRUE)
write.table(novelmirna, file="../2_predicted_novel_mirna.csv", sep = "\t", col.names=TRUE)
write.table(md, file = "../3_mature_mirna_detected.csv", sep = "\t", col.names=TRUE)
write.table(homologseed, "../4_predicted_novel_sigranfold_homologseed.txt", sep = "\t", col.names=TRUE)