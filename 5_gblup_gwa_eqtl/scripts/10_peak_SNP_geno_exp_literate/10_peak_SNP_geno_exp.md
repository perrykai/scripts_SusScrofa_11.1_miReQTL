**Script:** `10_peak_SNP_geno_exp.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`

**Date:**  `1/29/18`

**Input File Directory:**  

1. `../../4_dge_G_objects/`

2. `../../6_mirna_eQTL_target_prediction/`

**Input File(s):** 

1. `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`

2. `12_mrna_mirna_resid_exp.Rdata`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts/10_peak_SNP_geno_exp_literate/figure/`

**Output File(s):** `geno_exp_mir874.tiff`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives
## Install libraries


```r
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(regress)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)
```

```
## corrplot 0.84 loaded
```

```r
library(ggplot2)
```

## Load data

Load DV's eqtl functions:


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
```

Load required data:


```r
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
load("../../6_mirna_eQTL_target_prediction/12_mrna_mirna_resid_exp.Rdata")
ls()
```

```
##  [1] "absmap"               "add_legend"           "AddPosGene"          
##  [4] "dge"                  "distance"             "G"                   
##  [7] "inrange"              "manhpt"               "MSUPRP_miRNA"        
## [10] "peakrng"              "plot.GMA"             "rfit.mi"             
## [13] "rMfit"                "sigpval"              "stb"                 
## [16] "stb.nm"               "summary_MSUPRP_miRNA" "tbpos"               
## [19] "v"                    "wtcen"                "Z"                   
## [22] "zstandard"
```

## Analysis
---

Investigate peak SNP for miR-874: ALGA0016550. Do we see differences in miRNA expression segragating between genotypes?

How many animals have each genotype?


```r
MSUPRP_miRNA$map["ALGA0016550",]
```

```
##             chr       pos   snp
## ALGA0016550   2 139741258 [T/C]
```

```r
table(MSUPRP_miRNA$geno[,"ALGA0016550"])
```

```
## 
##  0  1  2 
## 16 84 74
```

```r
summary(v$E["ssc-miR-874",])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.588   4.194   4.733   4.747   5.269   7.905
```

```r
sum(names(MSUPRP_miRNA$geno[,"ALGA0016550"])!=names(v$E["ssc-miR-874",]))
```

```
## [1] 0
```

```r
mir874exp<-data.frame("geno"=as.character(MSUPRP_miRNA$geno[,"ALGA0016550"]), "exp"=v$E["ssc-miR-874",], check.rows=TRUE, stringsAsFactors=FALSE)

median(mir874exp[mir874exp$geno=="0","exp"])
```

```
## [1] 2.980244
```

```r
median(mir874exp[mir874exp$geno=="1","exp"])
```

```
## [1] 4.450585
```

```r
median(mir874exp[mir874exp$geno=="2","exp"])
```

```
## [1] 5.215782
```

```r
mir874exp$geno<-gsub("0", "AA", mir874exp$geno)
mir874exp$geno<-gsub("1", "AB", mir874exp$geno)
mir874exp$geno<-gsub("2", "BB", mir874exp$geno)
head(mir874exp)
```

```
##      geno      exp
## 1034   AB 5.422843
## 1036   BB 5.204700
## 1041   BB 5.334211
## 1049   AB 4.301332
## 1058   AB 6.464207
## 1060   AB 4.539884
```

---

Investigate peak SNP for miR-429: ALGA0118516. Do we see differences in miRNA expression segragating between genotypes?

How many animals have each genotype?


```r
MSUPRP_miRNA$map["ALGA0118516",]
```

```
##             chr      pos   snp
## ALGA0118516   6 63743210 [T/C]
```

```r
table(MSUPRP_miRNA$geno[,"ALGA0118516"])
```

```
## 
##   0   1   2 
## 108  59   7
```

```r
summary(v$E["ssc-miR-429",])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -1.6680 -0.7011  0.5851  0.5758  1.6580  3.7400
```

```r
sum(names(MSUPRP_miRNA$geno[,"ALGA0118516"])!=names(v$E["ssc-miR-429",]))
```

```
## [1] 0
```

```r
mir429exp<-data.frame("geno"=as.character(MSUPRP_miRNA$geno[,"ALGA0118516"]), "exp"=v$E["ssc-miR-429",], check.rows=TRUE, stringsAsFactors=FALSE)

median(mir429exp[mir429exp$geno=="0","exp"])
```

```
## [1] -0.2628393
```

```r
median(mir429exp[mir429exp$geno=="1","exp"])
```

```
## [1] 1.457471
```

```r
median(mir429exp[mir429exp$geno=="2","exp"])
```

```
## [1] 2.545934
```

```r
mir429exp$geno<-gsub("0", "AA", mir429exp$geno)
mir429exp$geno<-gsub("1", "AB", mir429exp$geno)
mir429exp$geno<-gsub("2", "BB", mir429exp$geno)
head(mir429exp)
```

```
##      geno        exp
## 1034   AA -1.6324392
## 1036   AB  0.7672951
## 1041   AA  1.8017162
## 1049   AB  0.3317057
## 1058   AA -1.5964887
## 1060   AA  1.9896869
```

## Visualize


```r
mir874genos<-ggplot(data=mir874exp, aes(x=geno,y=exp, col=geno)) +
geom_jitter(width=0.25) +
geom_segment(aes(x=0.75, y=median(mir874exp[mir874exp$geno=="AA","exp"]), xend=1.25, yend=median(mir874exp[mir874exp$geno=="AA","exp"])), color="black") +
geom_segment(aes(x=1.75, y=median(mir874exp[mir874exp$geno=="AB","exp"]), xend=2.25, yend=median(mir874exp[mir874exp$geno=="AB","exp"])), color="black") +
geom_segment(aes(x=2.75, y=median(mir874exp[mir874exp$geno=="BB","exp"]), xend=3.25, yend=median(mir874exp[mir874exp$geno=="BB","exp"])), color="black") +
labs(x="Genotype", y="Normalized miR-874 Expression", title="ALGA0016550 [T/C]", colour="Genotype") +
scale_color_manual(values=c("purple","orange","darkgreen"))+
theme_bw()
mir874genos
```

![plot of chunk geno_exp_mir874](figure/geno_exp_mir874-1.tiff)

```r
mir429genos<-ggplot(data=mir429exp, aes(x=geno,y=exp, col=geno)) +
geom_jitter(width=0.25) +
geom_segment(aes(x=0.75, y=median(mir429exp[mir429exp$geno=="AA","exp"]), xend=1.25, yend=median(mir429exp[mir429exp$geno=="AA","exp"])), color="black") +
geom_segment(aes(x=1.75, y=median(mir429exp[mir429exp$geno=="AB","exp"]), xend=2.25, yend=median(mir429exp[mir429exp$geno=="AB","exp"])), color="black") +
geom_segment(aes(x=2.75, y=median(mir429exp[mir429exp$geno=="BB","exp"]), xend=3.25, yend=median(mir429exp[mir429exp$geno=="BB","exp"])), color="black") +
labs(x="Genotype", y="Normalized miR-429 Expression", title="ALGA0118516 [T/C]", colour="Genotype") +
scale_color_manual(values=c("purple","orange","darkgreen"))+
theme_bw()
mir429genos
```

![plot of chunk geno_exp_mir429](figure/geno_exp_mir429-1.tiff)

## Save data
