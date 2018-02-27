#' **Script:** `8_Z_full.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12/19/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/3_miRNA_expression_mx`
#' 2. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`
#' 
#' **Input File(s):** 
#' 
#' 1. `1_mean_mature_mirna_exp.Rdata`
#'
#' 2. `MSUPRP_gpData_Ss11.Rdata`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `13_Z_full.Rdata`
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
#' The Objectives of this script are to: 
#' 
#' Create an updated Z matrix using maf filter of 1% to include peak SNP from pQTL analysis
#' 
#' ---
#' 
#' The genotype data will be compiled as follows:
#' 
#' 1. First, extract the F2 population genotypes from the MSUPRP$geno object and calculate AFs
#' 
#' 2. Subset the genotype data to include the 174 F2 animals in this analysis
#' 
#' 3. Filter those genotypes for MAF < 0.01 (being sure to take both ends of the distribution), removal of fixed SNPs, those located on sex chromosomes
#' 
#' 4. Remove these markers from both the SNP dataset and the map data.
#' 
#' ## Install libraries
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)

#' ## Load data
rm(list=ls())

#' Load Yeni's function for standardizing the Z matrix:
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")
ls()
#' Load the miRNA expression data:
load("../../3_miRNA_expression_mx/1_mean_mature_mirna_exp.Rdata")
#' Load the MSUPRP gpdata object:
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
ls()

#' ## Analysis
#' 
#' no.zero.dfmeanrcround is a matrix of gene expression, with dimensions genes x samples (335 x 174)
dim(no.zero.dfmeanrcround)
head(colnames(no.zero.dfmeanrcround))
head(rownames(no.zero.dfmeanrcround))

#' Create vector of animal IDs from the column names of the expression matrix for subsetting the gpdata object:
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
head(pigid)

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

#' Dimensions samples x categories (174 x categories), including growth_group as a factor:
#' 
#' 
#' Remove covariate data from animals not in this analysis:
todisc <- MSUPRP$covar$id[!MSUPRP$covar$id %in% pigid]
length(todisc)
redMSU <- discard.individuals(MSUPRP, todisc)

dim(redMSU$pedigree)
dim(redMSU$geno)

#' Filter the subset of 174 animals' genotypes for MAF < 0.01 (being sure to take both ends of the distribution), removal of fixed SNPs, and those on sex chromosomes.
#' 
#' Extract the genotype object from the reduced gpdata object (the 174 animals and all the SNP markers):
genomat<-redMSU$geno
dim(genomat)

#' Filter SNPs for minor allele frequency calculated using entire F2s (maf >= 0.01 retained):
#' How many SNPs have a af < 0.01 and >0.99?
#' 
#' We filter this way since we don't know which allele is "minor" in each SNP based on allele frequency (if A or B is the minor allele):

length(which(allele_freq<0.01)) + length(which(allele_freq>0.99))

length(names(which(allele_freq>=0.01 & allele_freq<=0.99)))

length(allele_freq) - length(names(which(allele_freq>=0.01 & allele_freq<=0.99)))

#' Retain all SNPs with af >= 0.01 and <=0.99 (maf < 0.01 discarded):
genomat <- genomat[,allele_freq>=0.01 & allele_freq <=0.99]

#' Dimensions of remaining SNP marker matrix:
dim(genomat)

#' Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):
sdv <- apply(genomat, 2, sd)
length(sdv)

if(sum(names(sdv) != colnames(genomat)) != 0) stop ("SNP names not the same between redMSU$geno and sdv")

sum(sdv == 0)

#' To make sure the proper number of SNPs were deleted, calculate the difference between the two datasets:
ncol(genomat) - sum(sdv == 0)

#' Remove fixed SNPs from genotype matrix:
genomat <- genomat[,sdv>0]

dim(genomat)

#' Eliminate markers on sex chromosomes:
#' 
#' JUNE 13, 2017: Error in synbreed package code; it saves the map object as a list when using "add.Markers"
#' function, which was used by DV to add the PRKAG3 markers.
#' 
#' So, need to extract the map object, convert to data.frame, and replace in the MSUPRP gpdata object, prior to filtering out SNPs.
map <- data.frame(chr=redMSU$map[[1]],pos=redMSU$map[[2]])
rownames(map) <- colnames(redMSU$geno)
head(map)
redMSU$map <- map
str(redMSU$map)

table(as.character(redMSU$map$chr))

#' Identify SNPs on sex chr (chr == 19)
xchr <- rownames(redMSU$map)[redMSU$map$chr == "X"]
ychr <- rownames(redMSU$map)[redMSU$map$chr == "Y"]
length(xchr)
length(ychr)

#' How many SNPs in my dataset are on the sex chr:
#' 
#' On X:
sum((colnames(genomat) %in% xchr))
#' On Y:
sum((colnames(genomat) %in% ychr))

#' Filter out the SNPs on the sex chromosome
genomatfil <- genomat[,!(colnames(genomat) %in% xchr)]
genomatfil <- genomatfil[,!(colnames(genomatfil) %in% ychr)]


#' This object contains the markers that are not fixed, have a maf > 0.10, and are not mapped to the sex chromosomes:
dim(genomatfil)

if (sum(colnames(genomatfil) %in% xchr) != 0) stop ("sex chromosome filter did not work correctly")
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

MSUPRP_miRNA$map$chr<-as.numeric(as.character(MSUPRP_miRNA$map$chr))
str(MSUPRP_miRNA$map)

#' The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 
#' 
summary_MSUPRP_miRNA<-summary(MSUPRP_miRNA)
summary_MSUPRP_miRNA$geno

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

#' Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:
Zfull<-zstandard(MSUPRP_miRNA$geno, alfreq=allele_freq, procedure="heterogeneous")
dim(Zfull)

G<-Zfull%*%t(Zfull)
dim(G)
summary(diag(G))

IQR(diag(G))

#' 
#' ## Save data
#' 
#' Save the full Z matrix:
save(Zfull, file = "../13_Z_full.Rdata")