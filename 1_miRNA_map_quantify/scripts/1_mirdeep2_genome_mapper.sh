#============================================
#  File: 1_mirdeep2_genome_mapper.sh
#  Directory Code:  /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/scripts
#  Date:  11/27/17 #UPDATE: Using miRDeep2/0.0.5 due to errors received when trying to run miRDeep2/0.0.7
#  Description:  This code maps the 174 smallRNA libraries to the sus scrofa reference genome 10.2.79 using miRDeep2 mapper module.
#                1. Uses config file to feed all 174 libraries to mapper module
#                2. Uses bowtie-index of reference genome to map reads using: bowtie –f –n 0 –e 80 –l 18 –a –m 5 –best –strata
#                3. Outputs a file of processed reads, and an .arf format alignment file for use with miRDeep2 core module
#============================================
#  Input File Directory:  /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/
#                         /mnt/research/ernstc_lab/references/Sscrofa_11.1/

#  Input File(s):   config_for_mapper.txt
#                   (6) genome index files, each using the prefix Sscrofa11.1_all_chr

#  Output File Directory: /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/1_mapper_output/
#

#  Output File(s): 1_mapper_output/mapper.log
#                  1_mapper_output/bowtie.log
#                  1_mapper_output/174_library_reads_processed.fa
#                  1_mapper_output/174_library_reads_mapped.arf

#============================================
#  Module used: miRDeep2/0.0.5
#============================================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=1:intel14,walltime=2:00:00,mem=3Gb
#PBS -N 1_miRDeep2_genome_mapper
#PBS -j oe
#PBS -o /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/OutputsErrors
#PBS -m a
#PBS -M perrykai@msu.edu

cd /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/

module load miRDeep2/0.0.5

mapper.pl /mnt/research/pigeqtl/analyses/microRNA/1_preprocess_fastq_files/8_collapsed_fasta_output_expression_matrix/config_for_mapper.txt -d -c -p /mnt/research/ernstc_lab/references/Sscrofa_11.1/bowtie1_index/Sscrofa11.1_all_chr -s /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/1_mapper_output/174_library_reads_processed.fa -t /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/1_mapper_output/174_library_reads_mapped.arf -m -n -v

mv mapper.log /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/1_mapper_output/mapper.log

mv bowtie.log /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/1_mapper_output/bowtie.log

qstat -f ${PBS_JOBID}
