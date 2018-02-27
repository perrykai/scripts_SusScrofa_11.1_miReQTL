#' **Script:** `1_miRNA_gblup_gwas.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12/5/17
#' 
#' **Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`
#' 
#' **Input File(s):** 
#' 
#' 1. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' 2. `4_normalized_dge_voom.Rdata`
#' 
#' 3. `5_Z_G_miRNA.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl`
#' 
#' **Output File(s):** 
#' 
#' 1. `1_gblup_full_results.Rdata`, `1_gblup_results_summary.Rdata`
#' 
#' 2. `2_gwa_results.Rdata`
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
#' The objective of this script is to conduct a gblup and gwa scan for the 295 MSUPRP F2 pig miRNA expression profiles.
#' 
#' This analysis will utilize the gwaR package developed by our lab group: <https://github.com/steibelj/gwaR>
#' 
#' ## Install libraries

rm(list=ls())
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts/")

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(plyr)

#' ## Load data
#' 
#' The final miRNA gpdata object with the voom-adjusted counts as phenotype data:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
#' The normalized dge object and voom centered precision weights:
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
#' The standardized Z and G matrices:
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()

#' ## Analysis
#' 
#' Create the design for the GBLUP model:
design_1<-c(~sex + growth_group)

#' Create miRnames, the list of the miRNAs to be input into the gblup analysis:
miRnames<-colnames(MSUPRP_miRNA$pheno)
length(miRnames)

#' ## Run gblup:
system.time({
rst.gblup<-lapply(miRnames, function(x) gblup(rsp=x, data=MSUPRP_miRNA, design=design_1, G=G, vdata=NULL, wt=wtcen, pos=c(T,T)))

names(rst.gblup)<-miRnames
})


#' Check standard error of G and wt, if NA eliminate from analysis:
#' 
#' First, what does this output look like:
rst.gblup$`ssc-let-7a`

#' coefm[6:7,2] is the standard errors of G and the wt
rst.gblup$`ssc-let-7a`$coefm[6:7,2]

std<-do.call(rbind, lapply(rst.gblup,function(x) x$coefm[6:7,2]))

#' Check how many NAs in standard error:
sum(is.na(rowSums(std)))

#' Retain only those miRNAs with a non-NA Standard Error
miRnames<-miRnames[!is.na(rowSums(std))]
length(miRnames)

#' Perform Likelihood Ratio Test on gblup results:
system.time({
like<-lapply(rst.gblup, lrt)
names(like)<-miRnames
})

#' See what these results look like:
like$`ssc-let-7a`

#' Multiple test corrections: FDR:
#' 
#' Obtain p-values from LRT
gblup.pval<-unlist(lapply(like, function(x) x$pvalue))

#' Create histogram of LRT p-values to see how much they deviate from expectation
#' 
hist(gblup.pval)

#' Calculate q-values from LRT p-values:
system.time({
	gblup.qval<-qvalue(gblup.pval,lambda=0)$qvalue
})

#' Total number of miRNAs with significant h2 at FDR < 0.05
sum(gblup.qval<0.05)

#' Matrix of standard errors of GBLUP:
std<-std[miRnames,]
colnames(std)<-paste("stdEr", colnames(std), sep="-")
head(std)

#' Matrix of gblup results:
summary.rst.gblup<-do.call(rbind, 
	lapply(miRnames, function(x) cbind(t(like[[x]]$vars[1]),
	 h2=like[[x]]$vars[1,3],
		like[[x]]$llik,
	lrtpvalue=like[[x]]$pvalue)))

rownames(summary.rst.gblup)<-miRnames

head(summary.rst.gblup)

#' Join it with the standard errors from the GBLUP:
if(sum(rownames(summary.rst.gblup) != rownames(std)) != 0) stop ("SNP ids not in same order between gblup summary mx and std object")

summary.rst.gblup<-cbind(summary.rst.gblup[,1:3], std, summary.rst.gblup[,4:ncol(summary.rst.gblup)], qvalue=gblup.qval)
head(summary.rst.gblup)

#' ---
#' 
#' ## Run GWA:
#' 
#' First, transpose standardized Z matrix:
Zt<-t(Z)
dim(Zt)

#' The following function performs these functions:
#' 
#' 1. Conducts the GWA analysis, using the results of the gblup:
#' 
#' 2. Calculates z-scores and calculates p-values from ghat and var(ghat) (Gualdrón Duarte 2014):
#' 
#' 3. Conducts Multiple Test Correction for the GWA via FDR:
#' 
#' 4. Returns the sign of the SNP effect (+ or -)
#' 
gwasum <- function(gbrst, ztp) {
	# Conduct GWA (returns vectors of ghat and var(ghat) for each miRNA):
	rst.gwa <- gwas(gbrst, x=ztp)
	# Calculate z-scores and get pvalues from ghat and var(ghat) (Gualdrón Duarte 2014):
	gwa.pval<- getpvalue(rst.gwa, log.p=F, is.z=F)
	# Multiple test correction for the GWA--FDR:
	gwa.qval<- qvalue(gwa.pval,lambda=0)
	# Determine if the SNP effect was + or -:
	sign<-ifelse(rst.gwa$ghat < 0, "-", "+") 

	return(data.frame(SNPid=rownames(rst.gwa), gwa.ghat=rst.gwa$ghat, gwa.pval, gwa.qval=gwa.qval$qvalues, SNP.sign=sign))

}

#' Apply the gwasum function to the list of gblup results, returning a data.frame by using ldply:
system.time({
rst.gwa<- ldply(rst.gblup, gwasum, ztp=Zt)
})

colnames(rst.gwa)[1] <- "miRNA"

#' How many rows should be in this data.frame?
length(unique(rst.gwa$SNPid)) * length(unique(rst.gwa$miRNA))

#' Check the dimensions of the data.frame and see results:
dim(rst.gwa)
head(rst.gwa)

#' ---

#' Check significance per gene at FDR < 0.05:
threshold <- 0.05
sig5<-(rst.gwa$gwa.qval < threshold)
length(sig5[sig5!=0])
sum(sig5)

gblup.h2.se.nms<-unique(rst.gwa[which(sig5!=0),"miRNA"])

gblup.h2.se<-lapply(gblup.h2.se.nms, function(x) varcomp(rst.gblup[[x]]))
names(gblup.h2.se)<-gblup.h2.se.nms

gblup.h2.se

#' The next step will be to further characterize the eQTL results.
#' 
#' ## Save data
#' 
save(rst.gblup, file="../1_gblup_full_results.Rdata")
save(summary.rst.gblup, gblup.h2.se, file="../1_gblup_results_summary.Rdata")
save(rst.gwa, file="../2_gwa_results.Rdata")