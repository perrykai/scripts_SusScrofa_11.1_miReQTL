
#=====================================================================
# This script runs: 4_mRNA_miRNA_corr.R
# Submited on: Mon Dec 18 10:54:02 EST 2017
#=====================================================================

#!/bin/sh  -login
#PBS -l nodes=1:ppn=11,walltime=00:15:00,mem=20G
#PBS -N 4_mRNA_miRNA_corr
#PBS -j oe
#PBS -m abe
#PBS -o /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts/4_mRNA_miRNA_corr

cd /mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts/4_mRNA_miRNA_corr

R -e 'library("knitr");knitr::spin ("4_mRNA_miRNA_corr.R")'

qstat -f ${PBS_JOBID}
