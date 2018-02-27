**Script:** `8_Z_full.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`

**Date:**  12/19/17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/3_miRNA_expression_mx`
2. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`

**Input File(s):** 

1. `1_mean_mature_mirna_exp.Rdata`

2. `MSUPRP_gpData_Ss11.Rdata`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`

**Output File(s):** 

1. `13_Z_full.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The Objectives of this script are to: 

Create an updated Z matrix using maf filter of 1% to include peak SNP from pQTL analysis

---

The genotype data will be compiled as follows:

1. First, extract the F2 population genotypes from the MSUPRP$geno object and calculate AFs

2. Subset the genotype data to include the 174 F2 animals in this analysis

3. Filter those genotypes for MAF < 0.01 (being sure to take both ends of the distribution), removal of fixed SNPs, those located on sex chromosomes

4. Remove these markers from both the SNP dataset and the map data.

## Install libraries


```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(synbreed)
library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)
```

## Load data


```r
rm(list=ls())
```

Load Yeni's function for standardizing the Z matrix:


```r
load("/mnt/research/pigsnp/GBLUP-based_GWA_Package/test_functions/testpkg/funct_eqtl.Rdata")
ls()
```

```
## [1] "absmap"     "add_legend" "AddPosGene" "manhpt"     "plot.GMA"  
## [6] "stb"        "tbpos"      "zstandard"
```

Load the miRNA expression data:


```r
load("../../3_miRNA_expression_mx/1_mean_mature_mirna_exp.Rdata")
```

Load the MSUPRP gpdata object:


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
ls()
```

```
##  [1] "absmap"                "add_legend"           
##  [3] "AddPosGene"            "manhpt"               
##  [5] "MSUPRP"                "MSUPRP168"            
##  [7] "no.zero.dfmeanrcround" "plot.GMA"             
##  [9] "stb"                   "tbpos"                
## [11] "zstandard"
```

## Analysis

no.zero.dfmeanrcround is a matrix of gene expression, with dimensions genes x samples (335 x 174)


```r
dim(no.zero.dfmeanrcround)
```

```
## [1] 335 174
```

```r
head(colnames(no.zero.dfmeanrcround))
```

```
## [1] "1034" "1036" "1041" "1049" "1058" "1060"
```

```r
head(rownames(no.zero.dfmeanrcround))
```

```
## [1] "ssc-let-7a"    "ssc-let-7c"    "ssc-let-7d-3p" "ssc-let-7d-5p"
## [5] "ssc-let-7e"    "ssc-let-7f"
```

Create vector of animal IDs from the column names of the expression matrix for subsetting the gpdata object:


```r
pigid <- colnames(no.zero.dfmeanrcround)
length(pigid)
```

```
## [1] 174
```

```r
head(pigid)
```

```
## [1] "1034" "1036" "1041" "1049" "1058" "1060"
```

---

### Genotype Filtering

First, extract the F2 population genotypes from the MSUPRP$geno object and calculate MAFs:


```r
dim(MSUPRP$geno)
```

```
## [1]  1015 43130
```

Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:


```r
geno_f2<-MSUPRP$geno[!((as.numeric(rownames(MSUPRP$geno))<=1000) | (as.numeric(rownames(MSUPRP$geno))>=6000)),]
dim(geno_f2)
```

```
## [1]   940 43130
```

Calculate allele frequency for all the F2 animals:


```r
allele_freq<-colMeans(geno_f2,na.rm=T)/2
length(allele_freq)
```

```
## [1] 43130
```

```r
summary(allele_freq)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.01011 0.28720 0.51060 0.50810 0.73240 0.98990
```

Dimensions samples x categories (174 x categories), including growth_group as a factor:


Remove covariate data from animals not in this analysis:


```r
todisc <- MSUPRP$covar$id[!MSUPRP$covar$id %in% pigid]
length(todisc)
```

```
## [1] 863
```

```r
redMSU <- discard.individuals(MSUPRP, todisc)

dim(redMSU$pedigree)
```

```
## [1] 174   5
```

```r
dim(redMSU$geno)
```

```
## [1]   174 43130
```

Filter the subset of 174 animals' genotypes for MAF < 0.01 (being sure to take both ends of the distribution), removal of fixed SNPs, and those on sex chromosomes.

Extract the genotype object from the reduced gpdata object (the 174 animals and all the SNP markers):


```r
genomat<-redMSU$geno
dim(genomat)
```

```
## [1]   174 43130
```

Filter SNPs for minor allele frequency calculated using entire F2s (maf >= 0.01 retained):
How many SNPs have a af < 0.01 and >0.99?

We filter this way since we don't know which allele is "minor" in each SNP based on allele frequency (if A or B is the minor allele):


```r
length(which(allele_freq<0.01)) + length(which(allele_freq>0.99))
```

```
## [1] 0
```

```r
length(names(which(allele_freq>=0.01 & allele_freq<=0.99)))
```

```
## [1] 43130
```

```r
length(allele_freq) - length(names(which(allele_freq>=0.01 & allele_freq<=0.99)))
```

```
## [1] 0
```

Retain all SNPs with af >= 0.01 and <=0.99 (maf < 0.01 discarded):


```r
genomat <- genomat[,allele_freq>=0.01 & allele_freq <=0.99]
```

Dimensions of remaining SNP marker matrix:


```r
dim(genomat)
```

```
## [1]   174 43130
```

Calculate standard deviation of SNP markers (if equal to 0, marker is fixed and must be removed):


```r
sdv <- apply(genomat, 2, sd)
length(sdv)
```

```
## [1] 43130
```

```r
if(sum(names(sdv) != colnames(genomat)) != 0) stop ("SNP names not the same between redMSU$geno and sdv")

sum(sdv == 0)
```

```
## [1] 32
```

To make sure the proper number of SNPs were deleted, calculate the difference between the two datasets:


```r
ncol(genomat) - sum(sdv == 0)
```

```
## [1] 43098
```

Remove fixed SNPs from genotype matrix:


```r
genomat <- genomat[,sdv>0]

dim(genomat)
```

```
## [1]   174 43098
```

Eliminate markers on sex chromosomes:

JUNE 13, 2017: Error in synbreed package code; it saves the map object as a list when using "add.Markers"
function, which was used by DV to add the PRKAG3 markers.

So, need to extract the map object, convert to data.frame, and replace in the MSUPRP gpdata object, prior to filtering out SNPs.


```r
map <- data.frame(chr=redMSU$map[[1]],pos=redMSU$map[[2]])
rownames(map) <- colnames(redMSU$geno)
head(map)
```

```
##             chr    pos
## MARC0044150   1 205163
## ASGA0000014   1 261794
## H3GA0000026   1 309120
## ASGA0000021   1 289363
## ALGA0000009   1 408640
## ALGA0000014   1 381075
```

```r
redMSU$map <- map
str(redMSU$map)
```

```
## 'data.frame':	43130 obs. of  2 variables:
##  $ chr: Factor w/ 173 levels "0","1","10","11",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ pos: int  205163 261794 309120 289363 408640 381075 373626 239523 455096 477400 ...
```

```r
table(as.character(redMSU$map$chr))
```

```
## 
##    1   10   11   12   13   14   15   16   17   18    2    3    4    5    6 
## 5092 1449 1577 1347 3312 3288 2325 1591 1408 1005 2830 2309 2975 1899 2724 
##    7    8    9    X    Y 
## 2665 2257 2674  394    9
```

Identify SNPs on sex chr (chr == 19)


```r
xchr <- rownames(redMSU$map)[redMSU$map$chr == "X"]
ychr <- rownames(redMSU$map)[redMSU$map$chr == "Y"]
length(xchr)
```

```
## [1] 394
```

```r
length(ychr)
```

```
## [1] 9
```

How many SNPs in my dataset are on the sex chr:

On X:


```r
sum((colnames(genomat) %in% xchr))
```

```
## [1] 394
```

On Y:


```r
sum((colnames(genomat) %in% ychr))
```

```
## [1] 9
```

Filter out the SNPs on the sex chromosome


```r
genomatfil <- genomat[,!(colnames(genomat) %in% xchr)]
genomatfil <- genomatfil[,!(colnames(genomatfil) %in% ychr)]
```

This object contains the markers that are not fixed, have a maf > 0.10, and are not mapped to the sex chromosomes:


```r
dim(genomatfil)
```

```
## [1]   174 42695
```

```r
if (sum(colnames(genomatfil) %in% xchr) != 0) stop ("sex chromosome filter did not work correctly")
if (sum(colnames(genomatfil) %in% ychr) != 0) stop ("sex chromosome filter did not work correctly")
if (sum(rownames(genomatfil) != colnames(no.zero.dfmeanrcround)) != 0) stop ("rownames of genotype data and count data not the same")

ncol(MSUPRP$geno)
```

```
## [1] 43130
```

```r
ncol(genomatfil)
```

```
## [1] 42695
```

Check on how many SNPs are being deleted:


```r
ncol(MSUPRP$geno) - ncol(genomatfil)
```

```
## [1] 435
```

Create the todel vector, containing the names of the SNPs in the gpdata geno object NOT found in the filtered genotype matrix:



```r
todel <- colnames(redMSU$geno)[!colnames(redMSU$geno) %in% colnames(genomatfil)]
length(todel)
```

```
## [1] 435
```

That means that, in total, we have removed 7165 SNPs from the gpdata geno object. 

Using discard.markers allows us to remove both the markers we don't want, and the map information we don't want, all in one step.


```r
MSUPRP_miRNA <- discard.markers(redMSU, todel)
dim(MSUPRP_miRNA$geno)
```

```
## [1]   174 42695
```

```r
MSUPRP_miRNA$map$chr<-as.numeric(as.character(MSUPRP_miRNA$map$chr))
str(MSUPRP_miRNA$map)
```

```
## 'data.frame':	42695 obs. of  2 variables:
##  $ chr: num  1 1 1 1 1 1 1 1 1 1 ...
##  $ pos: int  205163 261794 309120 289363 408640 381075 373626 239523 455096 477400 ...
```

The filtering of genotypes is complete. This will be used in the GBLUP and GWAS for the G matrix (genomic relationship matrix) after standardizing the SNPs. 



```r
summary_MSUPRP_miRNA<-summary(MSUPRP_miRNA)
summary_MSUPRP_miRNA$geno
```

```
## $nMarkers
## [1] 42695
## 
## $genotypes
##     0     1     2 
## 0.302 0.380 0.318 
## 
## $nNA
## [1] 0
## 
## $markerChr
## 
##    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 5092 2826 2309 2971 1898 2717 2664 2257 2667 1449 1577 1347 3310 3284 2325 
##   16   17   18 
## 1591 1407 1004 
## 
## $mappedMarkers
## [1] 42695
```

---

### Standardize the matrix of SNPs and create G matrix:

Using zstandard, the function created by Yeni Bernal Rubio:

This Z matrix is constructed using the allele frequencies of the entire F2 population and extracting the markers present in this analysis:

Need to filter the allele_freq object to retain the correct SNP markers.


```r
length(allele_freq)
```

```
## [1] 43130
```

```r
head(colnames(MSUPRP_miRNA$geno))
```

```
## [1] "MARC0044150" "ASGA0000014" "H3GA0000026" "ASGA0000021" "ALGA0000009"
## [6] "ALGA0000014"
```

```r
allele_freq<-allele_freq[colnames(MSUPRP_miRNA$geno)]
length(allele_freq)
```

```
## [1] 42695
```

```r
if(sum(names(allele_freq) != colnames(MSUPRP_miRNA$geno)) !=0) stop ("allele freq don't match SNP matrix")
```

Use the heterogeneous protocol from Yeni's function, computing Z with VanRaden et al 2008 option 2, like Jose Luis paper:


```r
Zfull<-zstandard(MSUPRP_miRNA$geno, alfreq=allele_freq, procedure="heterogeneous")
dim(Zfull)
```

```
## [1]   174 42695
```

```r
G<-Zfull%*%t(Zfull)
dim(G)
```

```
## [1] 174 174
```

```r
summary(diag(G))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.6971  0.9029  0.9472  0.9507  0.9965  1.1740
```

```r
IQR(diag(G))
```

```
## [1] 0.09360084
```


## Save data

Save the full Z matrix:


```r
save(Zfull, file = "../13_Z_full.Rdata")
```

