#' **Script:** `6_mirna_eqtl_var_comps_peaksnp.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12/12/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`
#' 
#' 3. `/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/`, `/mnt/research/pigeqtl/analyses/eQTL/paper/output/`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`
#' 
#' **Input File(s):** 
#' 
#' 1. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 2. `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`
#' 
#' 3. `funct_eqtl.Rdata`, `funct_eqtl.Rdata`, `MSUPRP_gpData_Ss11.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. ``
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
#' Testing variance components of GBLUP model incorporating peak SNP as fixed effect:
#' 
#' Original GBLUP model underestimates heritability, especially in cases where SNPs have large effects due to model assumptions that all SNPs have equal, small effect
#' 
#' Thus, when you include the SNPs that have a large effect on the miRNA expression (like an eQTL SNP), 
#' the variance components would shift when accounted for as an effect in the model and could result in an 'over-estimation' compared to the heritability.
#' 
#' Remember too that heritability is an estimation fully dependent on the model! 
#' 
#' To investigate this discrepancy in the variance components, separate the peak SNP genotypes prior to Z matrix standardization, then 
#' create the G_peak and G_bkg manually, then complete the gblup adding the G_peak as a fixed effect to the model.
#' 
#' This will allow us to see the estimate of b (different SNP effect), as it falls out in the GBLUP, and compare to the total variance explained by multiplying b^2 by variance of peak SNP
#' 
#' We concluded that the differences between methods in this estimation are nonsignificant
#' 
#' ## Install libraries
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)
library(synbreed)

#' ## Load data
#' 
#' Load Yeni's function for standardizing the Z matrix:
# load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")
#' Load DV's eqtl functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
#' Load required data, including MSUPRP gpdata object:
load("../2_gwa_results.Rdata")
load("../3_eqtl_summary_tables_maps.Rdata")
load("../../3_miRNA_expression_mx/1_mean_mature_mirna_exp.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
#' Load the MSUPRP gpdata object (for allele freq calculation):
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
ls()

#' ## Analysis
#' 
#' First, isolate which miRNA peak SNP(s) we want to utilize in the function:
mirpeaks[11,]

#' I will select miR-184 SNP ASGA0034057 
#' 
#' ---
#' 
#' ### Genotype Filtering
#' 
#' First, extract the F2 population genotypes from the MSUPRP$geno object and calculate MAFs:
dim(MSUPRP$geno)

#' Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:
geno_f2<-MSUPRP$geno[!((as.numeric(rownames(MSUPRP$geno))<=1000) | (as.numeric(rownames(MSUPRP$geno))>=6000)),]
dim(geno_f2)

#' Calculate allele frequency for all the F2 animals:
allele_freq<-colMeans(geno_f2,na.rm=T)/2
length(allele_freq)
summary(allele_freq)

allele_freq[c("ASGA0034057")]

#' Filter the subset of 174 animals' genotypes for MAF < 0.10 (being sure to take both ends of the distribution), removal of fixed SNPs, and those on sex chromosomes.
#' 
#' Create vector of animal IDs from the column names of the expression matrix for subsetting the gpdata object:
pigid <- rownames(dge$samples)
length(pigid)
head(pigid)
#' Remove covariate data from animals not in this analysis:
todisc <- MSUPRP$covar$id[!MSUPRP$covar$id %in% pigid]
length(todisc)
redMSU <- discard.individuals(MSUPRP, todisc)

dim(redMSU$pedigree)
dim(redMSU$geno)

#' Extract the genotype object from the reduced gpdata object (the 174 animals and all the SNP markers):
genomat<-redMSU$geno
dim(genomat)

#' Filter SNPs for minor allele frequency calculated using entire F2s (maf >= 0.10 retained):
#' How many SNPs have a af < 0.10 and >0.90?
#' 
#' We filter this way since we don't know which allele is "minor" in each SNP based on allele frequency (if A or B is the minor allele):

length(which(allele_freq<0.10)) + length(which(allele_freq>0.90))

length(names(which(allele_freq>=0.10 & allele_freq<=0.90)))

length(allele_freq) - length(names(which(allele_freq>=0.10 & allele_freq<=0.90)))

#' Retain all SNPs with af >= 0.10 and <=0.90 (maf < 0.10 discarded):
genomat <- genomat[,allele_freq>=0.10 & allele_freq <=0.90]

#' Dimensions of remaining SNP marker matrix:
dim(genomat)

#' Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):
sdv <- apply(genomat, 2, sd)
length(sdv)
summary(sdv)

sdv[c("ASGA0034057")]
sdv[c("ALGA0026452")]
sdv[c("ALGA0124095")]

if(sum(names(sdv) != colnames(genomat)) != 0) stop ("SNP names not the same between redMSU$geno and sdv")

sum(sdv == 0)

#' To make sure the proper number of SNPs were deleted, calculate the difference between the two datasets:
ncol(genomat) - sum(sdv == 0)

#' Remove fixed SNPs from genotype matrix:
genomat <- genomat[,sdv>0]

dim(genomat)

#' Eliminate markers on sex chromosomes:
#' 
str(redMSU$map)
redMSU$map$chr<-gsub("X", "19", redMSU$map$chr)
redMSU$map$chr<-gsub("Y", "20", redMSU$map$chr)
redMSU$map$chr<-as.numeric(as.character(redMSU$map$chr))
table(redMSU$map$chr)

#' Identify SNPs on sex chr (chr == 19)
sexchr <- rownames(redMSU$map)[redMSU$map$chr == 19]
length(sexchr)
ychr <- rownames(redMSU$map)[redMSU$map$chr == 20]
length(ychr)
#' How many SNPs in my dataset are on the sex chr:
sum((colnames(genomat) %in% sexchr))
sum((colnames(genomat) %in% ychr))

#' Filter out the SNPs on the sex chromosome
genomatfil <- genomat[,!(colnames(genomat) %in% sexchr)]
genomatfil <- genomatfil[,!(colnames(genomatfil) %in% ychr)]

#' This object contains the markers that are not fixed, have a maf > 0.10, and are not mapped to the sex chromosomes:
dim(genomatfil)

if (sum(colnames(genomatfil) %in% sexchr) != 0) stop ("sex chromosome filter did not work correctly")
if (sum(colnames(genomatfil) %in% ychr) != 0) stop ("sex chromosome filter did not work correctly")

if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

ncol(MSUPRP$geno)
ncol(genomatfil)

#' Check on how many SNPs are being deleted:
ncol(MSUPRP$geno) - ncol(genomatfil)

#' Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:
#' 
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)
#' That means that, in total, we have removed 7165 SNPs from the gpdata geno object. 
#' 
#' Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.
MSUPRP_miRNA <- discard.markers(redMSU, todel)
dim(MSUPRP_miRNA$geno)

#' The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 
#' 
summary_MSUPRP_miRNA<-summary(MSUPRP_miRNA)
summary_MSUPRP_miRNA$geno
summary_MSUPRP_miRNA$pedigree
summary_MSUPRP_miRNA$covar
str(summary_MSUPRP_miRNA$pheno)
head(summary_MSUPRP_miRNA$pheno)

#' ---
#' 
#' ### Standardize the matrix of SNPs and create G matrix:
#' 
#' Using zstandard, the function created by Yeni Bernal Rubio:
#' 
#' This Z matrix is constructed using the allele frequencies of the entire F2 population and extracting the markers present in this analysis:
#'
#' Need to filter the allele_freq object to retain the correct SNP markers.
length(allele_freq)
head(colnames(MSUPRP_miRNA$geno))
allele_freq<-allele_freq[colnames(MSUPRP_miRNA$geno)]
length(allele_freq)

if(sum(names(allele_freq) != colnames(MSUPRP_miRNA$geno)) !=0) stop ("allele freq don't match SNP matrix")

#' Separate the Z object into the objects required for this test:
#'
#' What I need: 
#'
#' 1. Z peak miR-184: column of SNP ASGA0034057
#' 
#' 2. Z background miR-184: Z mx excluding SNP ASGA0034057
#' 
#' 3. Z for full miR-eQTL peak range (2Mb)
#' 
#' 4. Z background for full miR-eQTL peak range (2Mb)
#'  
#' ### miR-184:
#' 
#' Extract peak SNP:
Z184<-as.matrix(MSUPRP_miRNA$geno[,"ASGA0034057"])
table(Z184)

bkgZ184<-MSUPRP_miRNA$geno[,-which(colnames(MSUPRP_miRNA$geno)=="ASGA0034057")]
dim(bkgZ184)
sum(colnames(bkgZ184)=="ASGA0034057")

#' Isolate the 2Mb region surrounding miR-184 peak
str(map.full)
map.full$chr<-as.numeric(as.character(map.full$chr))
map.nomir<-head(map.full, n=-17)
eqtl.2Mbregion<-inrange(chr=mirpeaks$chr.snp[11], start=mirpeaks$pos.snp[11]-1E6, end=mirpeaks$pos.snp[11]+1E6, map=map.nomir, single="pos", range=NULL)

Zpkrg<-MSUPRP_miRNA$geno[,rownames(eqtl.2Mbregion)]
dim(Zpkrg)

bkgZpkrg<-MSUPRP_miRNA$geno[,!colnames(MSUPRP_miRNA$geno)%in%colnames(Zpkrg)]
dim(bkgZpkrg)

#' Standardize matrices:
#' 
#' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
#' 
#' Peak Marker:
Z184<-zstandard(Z184, alfreq=allele_freq["ASGA0034057"], procedure="heterogeneous")
table(Z184)
colnames(Z184)<-"ASGA0034057"

bkgZ184<-zstandard(bkgZ184, alfreq=allele_freq[-which(names(allele_freq)=="ASGA0034057")], procedure="heterogeneous")

#' Peak Range:
Zpkrg<-zstandard(Zpkrg, alfreq=allele_freq[colnames(Zpkrg)], procedure="heterogeneous")

bkgZpkrg<-zstandard(bkgZpkrg, alfreq=allele_freq[!names(allele_freq)%in%colnames(Zpkrg)], procedure="heterogeneous")

#' Create G matrices:
#'
#' Peak Marker:
G184<-Z184%*%t(Z184)
dim(G184)
summary(diag(G184))

IQR(diag(G184))

bkgG184<-bkgZ184%*%t(bkgZ184)
dim(bkgG184)
summary(diag(bkgG184))

IQR(diag(bkgG184))

#' Peak Range:
Gpkrg<-Zpkrg%*%t(Zpkrg)
dim(Gpkrg)
summary(diag(Gpkrg))

IQR(diag(Gpkrg))

bkgGpkrg<-bkgZpkrg%*%t(bkgZpkrg)
dim(bkgGpkrg)
summary(diag(bkgGpkrg))

IQR(diag(bkgGpkrg))

#' ### MOVING ON to the GBLUP:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")

#' Add Z184 to the covariate matrix to include it as a fixed effect
sum(MSUPRP_miRNA$covar$id != rownames(Z184))
sum(MSUPRP_miRNA$covar$id != rownames(Zpkrg))

MSUPRP_miRNA$covar$"ASGA0034057" <- Z184

#' ## Run gblup:
#' 
#' First analysis: Include G matrix from peak range as random effect in gblup model, and compare results with output from test.peak function:

design_1<-list("~sex + growth_group", "~Gpkrg")

rst.gblup.pkrg<-gblup(rsp="ssc-miR-184", data=MSUPRP_miRNA, design=design_1, G=bkgGpkrg, vdata=list(Gpkrg=Gpkrg), wt=wtcen, pos=c(T,T,T))

rst.gblup.pkrg
varcomp(rst.gblup.pkrg)

#' Second analysis: Include peak marker as fixed effect in model, then extract the b statistic (SNP effect) and use to compare variances of eQTL. 
design_2<-c("~sex + growth_group + ASGA0034057")
rst.gblup.184<-gblup(rsp="ssc-miR-184", data=MSUPRP_miRNA, design=design_2, G=bkgG184, vdata=NULL, wt=wtcen, pos=c(T,T))

rst.gblup.184
varcomp(rst.gblup.184)

varq184<-as.numeric(var(Z184))*(rst.gblup.184$coefm[6,1]^2)
varq184

varq184/sum(varq184, rst.gblup.184$sigma[1], rst.gblup.184$sigma[2])

varcomp(rst.gblup.184)

#' ------------
#' 
#' 
#' Separate the Z object into the objects required for this test:
#'
#' What I need: 
#'
#' 1. Z peak miR-7135-3p: column of SNP ALGA0124095
#' 
#' 2. Z background miR-7135-3p: Z mx excluding SNP ALGA0124095
#' 
#' 3. Z for full miR-eQTL peak range (2Mb)
#' 
#' 4. Z background for full miR-eQTL peak range (2Mb)
#'  
#' ### miR-7135-3p:
mirpeaks[17,]

allele_freq[c("ALGA0124095")]


Z7135<-as.matrix(MSUPRP_miRNA$geno[,"ALGA0124095"])
table(Z7135)

bkgZ7135<-MSUPRP_miRNA$geno[,-which(colnames(MSUPRP_miRNA$geno)=="ALGA0124095")]
dim(bkgZ7135)
sum(colnames(bkgZ7135)=="ALGA0124095")

#' Isolate the 2Mb region surrounding miR-184 peak
map.nomir<-head(map.full, n=-17)
eqtl.2Mbregion<-inrange(chr=mirpeaks$chr.snp[17], start=mirpeaks$pos.snp[17]-1E6, end=mirpeaks$pos.snp[17]+1E6, map=map.nomir, single="pos", range=NULL)
dim(eqtl.2Mbregion)

Z7135pkrg<-MSUPRP_miRNA$geno[,rownames(eqtl.2Mbregion)]
dim(Z7135pkrg)

bkgZ7135pkrg<-MSUPRP_miRNA$geno[,!colnames(MSUPRP_miRNA$geno)%in%colnames(Z7135pkrg)]
dim(bkgZ7135pkrg)

#' Standardize matrices:
#' 
#' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
#' 
#' Peak Marker:
Z7135<-zstandard(Z7135, alfreq=allele_freq["ALGA0124095"], procedure="heterogeneous")
table(Z7135)
colnames(Z7135)<-"ALGA0124095"

bkgZ7135<-zstandard(bkgZ7135, alfreq=allele_freq[-which(names(allele_freq)=="ALGA0124095")], procedure="heterogeneous")
dim(bkgZ7135)

#' Peak Range:
Z7135pkrg<-zstandard(Z7135pkrg, alfreq=allele_freq[colnames(Z7135pkrg)], procedure="heterogeneous")

bkgZ7135pkrg<-zstandard(bkgZ7135pkrg, alfreq=allele_freq[!names(allele_freq)%in%colnames(Z7135pkrg)], procedure="heterogeneous")

#' Create G matrices:
#'
#' Peak Marker:
bkgG7135<-bkgZ7135%*%t(bkgZ7135)
dim(bkgG7135)
summary(diag(bkgG7135))

IQR(diag(bkgG7135))

#' Peak Range:
G7135pkrg<-Z7135pkrg%*%t(Z7135pkrg)
dim(G7135pkrg)
summary(diag(G7135pkrg))

IQR(diag(G7135pkrg))

bkgG7135pkrg<-bkgZ7135pkrg%*%t(bkgZ7135pkrg)
dim(bkgG7135pkrg)
summary(diag(bkgG7135pkrg))

IQR(diag(bkgG7135pkrg))

#' ### MOVING ON to the GBLUP:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")

#' Add Z184 to the covariate matrix to include it as a fixed effect
sum(MSUPRP_miRNA$covar$id != rownames(Z7135))
sum(MSUPRP_miRNA$covar$id != rownames(Z7135pkrg))

MSUPRP_miRNA$covar$"ALGA0124095" <- Z7135

#' ## Run gblup:
#' 
#' First analysis: Include G matrix from peak range as random effect in gblup model, and compare results with output from test.peak function:

design_1<-list("~sex + growth_group", "~G7135pkrg")

rst.gblup.7135pkrg<-gblup(rsp="ssc-miR-7135-3p", data=MSUPRP_miRNA, design=design_1, G=bkgG7135pkrg, vdata=list(G7135pkrg=G7135pkrg), wt=wtcen, pos=c(T,T,T))

rst.gblup.7135pkrg

varcomp(rst.gblup.7135pkrg)

#' Second analysis: Include peak marker as fixed effect in model, then extract the b statistic (SNP effect) and use to compare variances of eQTL. 

design_2<-c("~sex + growth_group + ALGA0124095")
rst.gblup.7135<-gblup(rsp="ssc-miR-7135-3p", data=MSUPRP_miRNA, design=design_2, G=bkgG7135, vdata=NULL, wt=wtcen, pos=c(T,T))

rst.gblup.7135

varq<-var(Z7135)*(rst.gblup.7135$coefm[6,1]^2)

varcomp(rst.gblup.7135)

# #' Separate the Z object into the objects required for this test:
# #'
# #' What I need: 
# #'
# #' 1. Z peak miR-184: column of SNP ALGA0041993
# #' 
# #' 2. Z background miR-184: Z mx excluding SNP ALGA0041993
# #' 
# #' 3. Z for full miR-eQTL peak range (2Mb)
# #' 
# #' 4. Z background for full miR-eQTL peak range (2Mb)
# #'  
# #' ### miR-184:
# #' 
# #' Extract peak SNP:
# Z184<-as.matrix(MSUPRP_miRNA$geno[,"ALGA0041993"])
# table(Z184)

# bkgZ184<-MSUPRP_miRNA$geno[,-which(colnames(MSUPRP_miRNA$geno)=="ALGA0041993")]
# dim(bkgZ184)
# sum(colnames(bkgZ184)=="ALGA0041993")

# #' Isolate the 2Mb region surrounding miR-184 peak
# map.nomir<-head(map.full, n=-17)
# eqtl.2Mbregion<-inrange(chr=mirpeaks$chr.snp[11], start=mirpeaks$pos.snp[11]-1E6, end=mirpeaks$pos.snp[11]+1E6, map=map.nomir, single="pos", range=NULL)

# Zpkrg<-MSUPRP_miRNA$geno[,rownames(eqtl.2Mbregion)]
# dim(Zpkrg)

# bkgZpkrg<-MSUPRP_miRNA$geno[,!colnames(MSUPRP_miRNA$geno)%in%colnames(Zpkrg)]
# dim(bkgZpkrg)

# #' Standardize matrices:
# #' 
# #' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
# #' 
# #' Peak Marker:
# Z184<-zstandard(Z184, alfreq=allele_freq["ALGA0041993"], procedure="heterogeneous")
# table(Z184)
# colnames(Z184)<-"ALGA0041993"

# bkgZ184<-zstandard(bkgZ184, alfreq=allele_freq[-which(names(allele_freq)=="ALGA0041993")], procedure="heterogeneous")

# #' Peak Range:
# Zpkrg<-zstandard(Zpkrg, alfreq=allele_freq[colnames(Zpkrg)], procedure="heterogeneous")

# bkgZpkrg<-zstandard(bkgZpkrg, alfreq=allele_freq[!names(allele_freq)%in%colnames(Zpkrg)], procedure="heterogeneous")

# #' Create G matrices:
# #'
# #' Peak Marker:
# G184<-Z184%*%t(Z184)
# dim(G184)
# summary(diag(G184))

# IQR(diag(G184))

# bkgG184<-bkgZ184%*%t(bkgZ184)
# dim(bkgG184)
# summary(diag(bkgG184))

# IQR(diag(bkgG184))

# #' Peak Range:
# Gpkrg<-Zpkrg%*%t(Zpkrg)
# dim(Gpkrg)
# summary(diag(Gpkrg))

# IQR(diag(Gpkrg))

# bkgGpkrg<-bkgZpkrg%*%t(bkgZpkrg)
# dim(bkgGpkrg)
# summary(diag(bkgGpkrg))

# IQR(diag(bkgGpkrg))

# #' ### MOVING ON to the GBLUP:
# load("../../3_build_dge_object_for_eqtl/3_msuprp_mirna_gpdata_pheno_counts.Rdata")
# load("../../3_build_dge_object_for_eqtl/4_normalized_dge_object_and_voom_output.Rdata")

# #' Add Z184 to the covariate matrix to include it as a fixed effect
# sum(MSUPRP_miRNA$covar$id != rownames(Z184))
# sum(MSUPRP_miRNA$covar$id != rownames(Zpkrg))

# MSUPRP_miRNA$covar$"ALGA0041993" <- Z184

# #' ## Run gblup:
# #' 
# #' First analysis: Include G matrix from peak range as random effect in gblup model, and compare results with output from test.peak function:

# design_1<-list("~sex + growth_group", "~Gpkrg")

# rst.gblup.pkrg<-gblup(rsp="ssc-miR-184", data=MSUPRP_miRNA, design=design_1, G=bkgGpkrg, vdata=list(Gpkrg=Gpkrg), wt=wtcen, pos=c(T,T,T))

# rst.gblup.pkrg
# varcomp(rst.gblup.pkrg)

# #' Second analysis: Include peak marker as fixed effect in model, then extract the b statistic (SNP effect) and use to compare variances of eQTL. 
# design_2<-c("~sex + growth_group + ALGA0041993")
# rst.gblup.184<-gblup(rsp="ssc-miR-184", data=MSUPRP_miRNA, design=design_2, G=bkgG184, vdata=NULL, wt=wtcen, pos=c(T,T))

# rst.gblup.184

# varq184<-as.numeric(var(Z184))*(rst.gblup.184$coefm[6,1]^2)
# varq184

# varq184/sum(varq184, rst.gblup.184$sigma[1], rst.gblup.184$sigma[2])

# varcomp(rst.gblup.184)
