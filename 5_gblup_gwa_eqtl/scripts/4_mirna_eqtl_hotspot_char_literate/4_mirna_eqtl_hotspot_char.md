**Script:** `4_mirna_eqtl_hotspot_char.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`

**Date:**  12/11/17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects`

2. `/mnt/research/pigeqtl/analyses/eQTL/paper/code/corrected-Z/PRKAG3_eQTL/`

**Input File(s):** 

1. `6_mirna_precursor_annot_ssc11.Rdata`

1. `3_msuprp_mirna_gpdata.Rdata`

2. `4_normalized_dge_voom.Rdata`

3. `5_Z_G_miRNA.Rdata`

4. `MSUPRP_gpData_Ss11.Rdata`

**Output File Directory:** 

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl`

**Output File(s):** 

1. `6_hotspot_miRNA_corr_rst.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to investigate the putative miRNA eQTL hotspots to determine if they are truly hotspots, or spurious associations due to high correlation of miRNA expression or genotype data.

To do this, miRNA expression will be correlated between hotspot miRNAs utilizing the log-cpm (v$E) . 
Pearson correlation will be used. 

I will also investigate the genomic origins of these miRNAs and SNPs, to see if the miRNAs are coming from similar genomic regions, if they have similar seed sequences, etc.

## Install libraries


```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

rm(list=ls())

library(regress)
library(gwaR)
library(limma)
library(edgeR)
library(parallel)
library(qvalue)
library(corrplot)
```

```
## corrplot 0.84 loaded
```

## Load data

Load microRNA expression data


```r
load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
```

Load the MSUPRP gpdata object (for allele freq calculation):


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")

ls()
```

```
##  [1] "annotation"           "dge"                  "G"                   
##  [4] "MSUPRP"               "MSUPRP168"            "MSUPRP_miRNA"        
##  [7] "summary_MSUPRP_miRNA" "v"                    "wtcen"               
## [10] "Z"
```

## Analysis

Extract the names of the miRNA in the hotspots:

Two SNP ("MARC0027291" "MARC0093624") were associated with five miRNAs, overlapping to some degree with eachother:


```r
htspt.mirna5 <- as.character(c("ssc-let-7d-5p","ssc-let-7g","ssc-miR-1468","ssc-miR-345-3p","ssc-miR-95","ssc-miR-9843-3p"))
htspt.mirna5
```

```
## [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-1468"    "ssc-miR-345-3p" 
## [5] "ssc-miR-95"      "ssc-miR-9843-3p"
```

Correlate the miRNA expression profiles using the log-cpm counts of miRNA expression:


```r
data5<-t(v$E[htspt.mirna5,])
head(data5)
```

```
##      ssc-let-7d-5p ssc-let-7g ssc-miR-1468 ssc-miR-345-3p ssc-miR-95
## 1034      11.63362   13.39516     8.860415       6.601181   11.44954
## 1036      11.10306   13.15001     7.920425       6.476953   10.90727
## 1041      11.28869   13.24917     8.282074       6.288982   10.94674
## 1049      11.33382   13.15911     8.066415       6.763551   10.89172
## 1058      11.16527   12.88285     8.728817       6.340149   11.20300
## 1060      11.86027   13.64016     8.511007       5.896578   11.26415
##      ssc-miR-9843-3p
## 1034        7.893082
## 1036        7.694085
## 1041        7.390966
## 1049        7.804735
## 1058        8.222093
## 1060        6.215755
```

```r
cormx5<-cor(data5)
head(cormx5)
```

```
##                 ssc-let-7d-5p ssc-let-7g ssc-miR-1468 ssc-miR-345-3p
## ssc-let-7d-5p      1.00000000  0.8536184  0.742557670    0.081568516
## ssc-let-7g         0.85361842  1.0000000  0.624759973    0.111621966
## ssc-miR-1468       0.74255767  0.6247600  1.000000000    0.008478271
## ssc-miR-345-3p     0.08156852  0.1116220  0.008478271    1.000000000
## ssc-miR-95         0.28317349  0.3019469  0.424386351   -0.045211686
## ssc-miR-9843-3p   -0.43635674 -0.3532372 -0.291866596    0.453703788
##                   ssc-miR-95 ssc-miR-9843-3p
## ssc-let-7d-5p    0.283173491    -0.436356740
## ssc-let-7g       0.301946906    -0.353237197
## ssc-miR-1468     0.424386351    -0.291866596
## ssc-miR-345-3p  -0.045211686     0.453703788
## ssc-miR-95       1.000000000     0.001440106
## ssc-miR-9843-3p  0.001440106     1.000000000
```

Function to perform correlation analysis of miRNA to miRNA espression and export a summary of correlation analysis


```r
cor.exp <- function(var1, var2, data, ...){
	x <- cor.test(data[,as.character(var1)], data[,as.character(var2)], ...)
	rst <- data.frame(var1, var2, cor=x$estimate, 
		t=x$statistic, 
		pvalue=x$p.value)
	return(rst)
}
```

---

Correlation test of 5 miRNAs:


```r
thres<-0.05
mir5.htspt<-list()
for(i in 1:length(colnames(data5))){
mir5.htspt[[i]]<-data.frame(var1=colnames(data5)[i], var2=colnames(data5)[-i])
}

vars5<-do.call(rbind, mir5.htspt)

mir.cor5<- do.call(rbind, mapply(cor.exp, var1=vars5[,1], var2=vars5[,2], MoreArgs=c(list(data=data5), alternative="two.sided", method="p"), SIMPLIFY=FALSE))
```

Perform multiple test correction (FDR) on correlation analysis:


```r
mir.cor5$qval<-qvalue(mir.cor5$pvalue)$qvalue
mir.cor5$pi0<-qvalue(mir.cor5$pvalue)$pi0
mir.cor5
```

```
##                  var1            var2          cor           t
## cor     ssc-let-7d-5p      ssc-let-7g  0.853618417 21.49189018
## cor1    ssc-let-7d-5p    ssc-miR-1468  0.742557670 14.53988103
## cor2    ssc-let-7d-5p  ssc-miR-345-3p  0.081568516  1.07333770
## cor3    ssc-let-7d-5p      ssc-miR-95  0.283173491  3.87228309
## cor4    ssc-let-7d-5p ssc-miR-9843-3p -0.436356740 -6.36022669
## cor5       ssc-let-7g   ssc-let-7d-5p  0.853618417 21.49189018
## cor6       ssc-let-7g    ssc-miR-1468  0.624759973 10.49369013
## cor7       ssc-let-7g  ssc-miR-345-3p  0.111621966  1.47311424
## cor8       ssc-let-7g      ssc-miR-95  0.301946906  4.15387996
## cor9       ssc-let-7g ssc-miR-9843-3p -0.353237197 -4.95189209
## cor10    ssc-miR-1468   ssc-let-7d-5p  0.742557670 14.53988103
## cor11    ssc-miR-1468      ssc-let-7g  0.624759973 10.49369013
## cor12    ssc-miR-1468  ssc-miR-345-3p  0.008478271  0.11119548
## cor13    ssc-miR-1468      ssc-miR-95  0.424386351  6.14675903
## cor14    ssc-miR-1468 ssc-miR-9843-3p -0.291866596 -4.00204752
## cor15  ssc-miR-345-3p   ssc-let-7d-5p  0.081568516  1.07333770
## cor16  ssc-miR-345-3p      ssc-let-7g  0.111621966  1.47311424
## cor17  ssc-miR-345-3p    ssc-miR-1468  0.008478271  0.11119548
## cor18  ssc-miR-345-3p      ssc-miR-95 -0.045211686 -0.59355266
## cor19  ssc-miR-345-3p ssc-miR-9843-3p  0.453703788  6.67704916
## cor20      ssc-miR-95   ssc-let-7d-5p  0.283173491  3.87228309
## cor21      ssc-miR-95      ssc-let-7g  0.301946906  4.15387996
## cor22      ssc-miR-95    ssc-miR-1468  0.424386351  6.14675903
## cor23      ssc-miR-95  ssc-miR-345-3p -0.045211686 -0.59355266
## cor24      ssc-miR-95 ssc-miR-9843-3p  0.001440106  0.01888683
## cor25 ssc-miR-9843-3p   ssc-let-7d-5p -0.436356740 -6.36022669
## cor26 ssc-miR-9843-3p      ssc-let-7g -0.353237197 -4.95189209
## cor27 ssc-miR-9843-3p    ssc-miR-1468 -0.291866596 -4.00204752
## cor28 ssc-miR-9843-3p  ssc-miR-345-3p  0.453703788  6.67704916
## cor29 ssc-miR-9843-3p      ssc-miR-95  0.001440106  0.01888683
##             pvalue         qval pi0
## cor   0.000000e+00 0.000000e+00   1
## cor1  0.000000e+00 0.000000e+00   1
## cor2  2.846232e-01 3.557790e-01   1
## cor3  1.529773e-04 2.294660e-04   1
## cor4  1.756288e-09 5.268864e-09   1
## cor5  0.000000e+00 0.000000e+00   1
## cor6  0.000000e+00 0.000000e+00   1
## cor7  1.425479e-01 1.943835e-01   1
## cor8  5.140065e-05 9.637623e-05   1
## cor9  1.744514e-06 3.738245e-06   1
## cor10 0.000000e+00 0.000000e+00   1
## cor11 0.000000e+00 0.000000e+00   1
## cor12 9.115910e-01 9.767046e-01   1
## cor13 5.349823e-09 1.337456e-08   1
## cor14 9.319294e-05 1.553216e-04   1
## cor15 2.846232e-01 3.557790e-01   1
## cor16 1.425479e-01 1.943835e-01   1
## cor17 9.115910e-01 9.767046e-01   1
## cor18 5.535911e-01 6.387590e-01   1
## cor19 3.231275e-10 1.211728e-09   1
## cor20 1.529773e-04 2.294660e-04   1
## cor21 5.140065e-05 9.637623e-05   1
## cor22 5.349823e-09 1.337456e-08   1
## cor23 5.535911e-01 6.387590e-01   1
## cor24 9.849533e-01 9.849533e-01   1
## cor25 1.756288e-09 5.268864e-09   1
## cor26 1.744514e-06 3.738245e-06   1
## cor27 9.319294e-05 1.553216e-04   1
## cor28 3.231275e-10 1.211728e-09   1
## cor29 9.849533e-01 9.849533e-01   1
```

Significantly correlated miRNA pairs:  


```r
mir.cor5[mir.cor5$qval<thres,]
```

```
##                  var1            var2        cor         t       pvalue
## cor     ssc-let-7d-5p      ssc-let-7g  0.8536184 21.491890 0.000000e+00
## cor1    ssc-let-7d-5p    ssc-miR-1468  0.7425577 14.539881 0.000000e+00
## cor3    ssc-let-7d-5p      ssc-miR-95  0.2831735  3.872283 1.529773e-04
## cor4    ssc-let-7d-5p ssc-miR-9843-3p -0.4363567 -6.360227 1.756288e-09
## cor5       ssc-let-7g   ssc-let-7d-5p  0.8536184 21.491890 0.000000e+00
## cor6       ssc-let-7g    ssc-miR-1468  0.6247600 10.493690 0.000000e+00
## cor8       ssc-let-7g      ssc-miR-95  0.3019469  4.153880 5.140065e-05
## cor9       ssc-let-7g ssc-miR-9843-3p -0.3532372 -4.951892 1.744514e-06
## cor10    ssc-miR-1468   ssc-let-7d-5p  0.7425577 14.539881 0.000000e+00
## cor11    ssc-miR-1468      ssc-let-7g  0.6247600 10.493690 0.000000e+00
## cor13    ssc-miR-1468      ssc-miR-95  0.4243864  6.146759 5.349823e-09
## cor14    ssc-miR-1468 ssc-miR-9843-3p -0.2918666 -4.002048 9.319294e-05
## cor19  ssc-miR-345-3p ssc-miR-9843-3p  0.4537038  6.677049 3.231275e-10
## cor20      ssc-miR-95   ssc-let-7d-5p  0.2831735  3.872283 1.529773e-04
## cor21      ssc-miR-95      ssc-let-7g  0.3019469  4.153880 5.140065e-05
## cor22      ssc-miR-95    ssc-miR-1468  0.4243864  6.146759 5.349823e-09
## cor25 ssc-miR-9843-3p   ssc-let-7d-5p -0.4363567 -6.360227 1.756288e-09
## cor26 ssc-miR-9843-3p      ssc-let-7g -0.3532372 -4.951892 1.744514e-06
## cor27 ssc-miR-9843-3p    ssc-miR-1468 -0.2918666 -4.002048 9.319294e-05
## cor28 ssc-miR-9843-3p  ssc-miR-345-3p  0.4537038  6.677049 3.231275e-10
##               qval pi0
## cor   0.000000e+00   1
## cor1  0.000000e+00   1
## cor3  2.294660e-04   1
## cor4  5.268864e-09   1
## cor5  0.000000e+00   1
## cor6  0.000000e+00   1
## cor8  9.637623e-05   1
## cor9  3.738245e-06   1
## cor10 0.000000e+00   1
## cor11 0.000000e+00   1
## cor13 1.337456e-08   1
## cor14 1.553216e-04   1
## cor19 1.211728e-09   1
## cor20 2.294660e-04   1
## cor21 9.637623e-05   1
## cor22 1.337456e-08   1
## cor25 5.268864e-09   1
## cor26 3.738245e-06   1
## cor27 1.553216e-04   1
## cor28 1.211728e-09   1
```

---

Investigate correlation of genotypes between the 4 putative eQTL hotspots:


```r
snp.htspt<- c("MARC0027291", "MARC0093624")
```

Examine allele frequencies of the hotspot SNPs

Remember that all the animals with ID bigger than 6000 are F0 and smaller than 1000 are F1, thus remove those animals to retain the F2s:


```r
geno_f2<-MSUPRP$geno[!((as.numeric(rownames(MSUPRP$geno))<=1000) | (as.numeric(rownames(MSUPRP$geno))>=6000)),]
dim(geno_f2)
```

```
## [1]   940 43130
```

Subset the hotspot SNPs:


```r
geno.htspt<-geno_f2[,snp.htspt]
dim(geno.htspt)
```

```
## [1] 940   2
```

Subset the 174 animals:


```r
geno.htspt<-geno.htspt[rownames(data5),]
dim(geno.htspt)
```

```
## [1] 174   2
```

```r
head(geno.htspt)
```

```
##      MARC0027291 MARC0093624
## 1034           2           2
## 1036           2           2
## 1041           2           2
## 1049           2           2
## 1058           2           2
## 1060           1           1
```

```r
if(sum(rownames(geno.htspt)!=rownames(data5)) != 0) stop ("Animal IDs not correct for genotype object")

geno.cor<-cor(geno.htspt)
```


Calculate allele frequency for all the F2 animals:


```r
allele_freq<-colMeans(geno.htspt,na.rm=T)/2
allele_freq
```

```
## MARC0027291 MARC0093624 
##   0.8735632   0.8965517
```

Minor allele frequency:


```r
maf<-ifelse(allele_freq>0.5, (1- allele_freq), allele_freq)
maf
```

```
## MARC0027291 MARC0093624 
##   0.1264368   0.1034483
```

---

Identify the associated miRNAs in the annotation file, to query miRBase and the literature for information on their effects 


```r
head(annotation)
```

```
##             miRNA       Ensembl.Gene           Transcript chr0    start
## 1   ssc-let-7d-5p ENSSSCG00000019364 ENSSSCT00000020959.2    3 43468501
## 2      ssc-let-7g ENSSSCG00000019306 ENSSSCT00000020901.2   13 34406135
## 3     ssc-miR-128 ENSSSCG00000019503 ENSSSCT00000021098.2   15 16273557
## 4     ssc-miR-128 ENSSSCG00000019835 ENSSSCT00000021430.2   13 21038442
## 5 ssc-miR-1306-3p ENSSSCG00000019500 ENSSSCT00000021095.3   14 51473918
## 6  ssc-miR-140-5p ENSSSCG00000028710 ENSSSCT00000030513.2    6 17077517
##        end width strand                     type     Alias
## 1 43468596    96      + miRNA_primary_transcript MI0022120
## 2 34406222    88      - miRNA_primary_transcript MI0013087
## 3 16273638    82      - miRNA_primary_transcript MI0002451
## 4 21038525    84      + miRNA_primary_transcript MI0013094
## 5 51473997    80      + miRNA_primary_transcript MI0013148
## 6 17077610    94      - miRNA_primary_transcript MI0002437
```

```r
prec.annot<-annotation[annotation$miRNA==c("ssc-let-7d-5p", "ssc-let-7g","ssc-miR-1468", "ssc-miR-345-3p", "ssc-miR-95", "ssc-miR-9843-3p"),]
```

Information extracted from target prediction input:

miRNA           | Seed Sequence (Human)
--------------- | -------------
let-7-5p/98-5p  | GAGGUAG
miR-1468-5p	   | UCCGUUU
miR-345-3p	   | CCCUGAA
miR-95-3p	   | UCAACGG

Also see `miRNA_targets_common-1.tiff` to compare the targets in common expressed in this dataset between miRNAs.

## Visualize

Correlation plots of 4 and 5 hotspot-associated miRNAs:

5 miRNAs:


```r
corrplot.mixed(cormx5, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Five hotspot-associated miRNAs")
```

![plot of chunk miRNA5_hotspots](figure/miRNA5_hotspots-1.tiff)

```r
corrplot.mixed(geno.cor, mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title="Putative Hotspot SNP Genotype Correlation")
```

![plot of chunk geno_cor_hotspots](figure/geno_cor_hotspots-1.tiff)

## Save data


```r
save(mir.cor5, geno.cor, maf, file="../6_hotspot_miRNA_corr_rst.Rdata")
```

