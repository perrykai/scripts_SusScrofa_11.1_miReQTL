**Script:** `7_mirna_mrna_pheno_corr.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`

**Date:**  12/21/17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11/`

2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`, `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`, `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`

**Input File(s):** 

1. `voom.Rdata`, `MSUPRP_gpData_Ss11.Rdata`, `G_Z_matrix_Ss11.Rdata`

2. `5_filtered_targets_exp_rst.Rdata`, `13_target_mrna_coloc_pqtl.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`

3. `12_mrna_mirna_resid_exp.Rdata`, `9_mireqtl_pqtl_coloc_peaks.Rdata`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

**Output File(s):** ``

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to correlate residual miRNA expression and residual mRNA expression with the phenotypes of interest

(meaning pQTL phenotypes co-localizing with miR-eQTL target genes)

These results will be utilized for supporting miR-874 co-localized targets results.

## Install libraries


```r
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

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


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11/G_Z_matrix_Ss11.Rdata")
ls()
```

```
## [1] "dge"       "G"         "MSUPRP"    "MSUPRP168" "v"         "wcen"     
## [7] "Z"
```

Rename R object to differentiate between mRNA and microRNA


```r
Mdge <- dge
Mv <- v
Mwcen <- wcen
MG <- G
```

Remove R objects that will not be used in this analysis


```r
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168", "MSUPRP")))
ls()
```

```
## [1] "Mdge"      "MG"        "MSUPRP"    "MSUPRP168" "Mv"        "Mwcen"
```

Load microRNA targets


```r
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/5_filtered_targets_exp_rst.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/13_target_mrna_coloc_pqtl.Rdata")
```

Load microRNA expression data


```r
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()
```

```
##  [1] "annot"                "coloc"                "dge"                 
##  [4] "G"                    "Mdge"                 "MG"                  
##  [7] "MSUPRP"               "MSUPRP168"            "MSUPRP_miRNA"        
## [10] "Mv"                   "Mwcen"                "negcoloc"            
## [13] "summary_MSUPRP_miRNA" "targets.exp"          "targets.exp.sum"     
## [16] "v"                    "wtcen"                "Z"
```

Rename R object to differentiate microRNA data from mRNA data


```r
dge.mi <- dge
v.mi <- v
wcen.mi <- wtcen
G.mi <- G
```

Retain only the objects needed for the analysis


```r
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168", "MSUPRP",
	"dge.mi","v.mi","wcen.mi","G.mi","MSUPRP_miRNA", "targets.exp", "coloc", "negcoloc")))

load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/12_mrna_mirna_resid_exp.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/9_mireqtl_pqtl_coloc_peaks.Rdata")
ls()
```

```
##  [1] "coloc"        "dge.mi"       "G.mi"         "Mdge"        
##  [5] "MG"           "MSUPRP"       "MSUPRP168"    "MSUPRP_miRNA"
##  [9] "Mv"           "Mwcen"        "negcoloc"     "regul"       
## [13] "rfit.mi"      "rMfit"        "targets.exp"  "v.mi"        
## [17] "wcen.mi"
```

## Analysis

Extract the phenotype object and subset the 166 animals:


```r
pheno166<-MSUPRP$pheno[,,1]
pheno166<-pheno166[colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)],]
dim(pheno166)
```

```
## [1] 166  67
```

Subset the phenotypes colocalized:


```r
nm.pheno166<-as.character(unique(negcoloc$pheno))

colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)]
```

```
##   [1] "1034" "1036" "1041" "1049" "1058" "1060" "1080" "1082" "1085" "1096"
##  [11] "1100" "1107" "1110" "1111" "1113" "1116" "1123" "1134" "1136" "1145"
##  [21] "1147" "1154" "1158" "1170" "1177" "1179" "1192" "1194" "1197" "1199"
##  [31] "1205" "1207" "1237" "1240" "1242" "1265" "1267" "1278" "1282" "1291"
##  [41] "1295" "1300" "1304" "1321" "1323" "1424" "1425" "1426" "1431" "1434"
##  [51] "1435" "1444" "1445" "1449" "1458" "1482" "1484" "1491" "1493" "1502"
##  [61] "1504" "1510" "1512" "1517" "1523" "1529" "1532" "1534" "1537" "1543"
##  [71] "1578" "1580" "1589" "1592" "1593" "1594" "1625" "1627" "1638" "1640"
##  [81] "1644" "1652" "1662" "1669" "1677" "1685" "1687" "1695" "1697" "1746"
##  [91] "1758" "1776" "1778" "1782" "1784" "1785" "1789" "1793" "1798" "1800"
## [101] "1818" "1819" "1820" "1833" "1836" "1839" "1843" "1844" "1879" "1881"
## [111] "1884" "1886" "1889" "1891" "1903" "1904" "1907" "1910" "1914" "1916"
## [121] "1928" "1930" "1965" "1971" "1976" "1980" "1989" "1991" "1999" "2003"
## [131] "2018" "2020" "2022" "2024" "2026" "2027" "2029" "2030" "2064" "2071"
## [141] "2073" "2076" "2094" "2100" "2118" "2119" "2120" "2123" "2131" "2135"
## [151] "2141" "2143" "2152" "2154" "2164" "2168" "2195" "2197" "2229" "2231"
## [161] "2261" "2263" "2297" "2303" "2311" "2317"
```

```r
rfit.misub<-rfit.mi[colnames(v.mi$E)[colnames(v.mi$E) %in% colnames(Mv$E)],]
dim(rfit.misub)
```

```
## [1] 166  15
```

Build list of miRNAs and associated pQTL of interest (either miR-eQTL colocalizes, or target genes colocalize):


```r
mirqtl<-c("juiciness", "driploss", "ph_24h", "protein", "cook_yield", "WBS", "tenderness", "overtend")
pheno.mir.list<-list("last_lum", "num_ribs", mirqtl, mirqtl, mirqtl, mirqtl, mirqtl, nm.pheno166)
names(pheno.mir.list)<-c("ssc-miR-190b","ssc-miR-184","ssc-miR-345-3p","ssc-let-7d-5p","ssc-let-7g","ssc-miR-95","ssc-miR-1468","ssc-miR-874")
rm(mirqtl)
```

Build list of mRNAs and associated pQTL of interest (negatively correlated target genes that co-localize with pQTL phenotypes)

Target genes co-localized with pQTL:


```r
unique(negcoloc$ID)
```

```
##  [1] "XLOC_002346" "XLOC_002373" "XLOC_011766" "XLOC_011767" "XLOC_011794"
##  [6] "XLOC_013026" "XLOC_013029" "XLOC_013030" "XLOC_019568" "XLOC_021684"
## [11] "XLOC_021691" "XLOC_021695" "XLOC_022575" "XLOC_022612" "XLOC_022617"
## [16] "XLOC_022870" "XLOC_022875" "XLOC_023083" "XLOC_023143" "XLOC_023161"
## [21] "XLOC_003117" "XLOC_008680" "XLOC_008685" "XLOC_008850" "XLOC_008918"
## [26] "XLOC_009274" "XLOC_009276" "XLOC_009284" "XLOC_009349" "XLOC_009537"
## [31] "XLOC_009607"
```

Extract these genes from the rMfit object:


```r
sub.rMfit<-rMfit[,colnames(rMfit)%in%unique(negcoloc$ID)]
dim(sub.rMfit)
```

```
## [1] 166  31
```

Subset the miR-874 target genes:


```r
pheno.mir874<-pheno.mir.list$`ssc-miR-874`
```

Function to perform correlation analysis of mRNA and miRNA espression and export a summary of correlation analysis


```r
cor.exp <- function(target, transc, pheno, ...){
	x <- cor.test(transc, pheno, ...)
	rst <- data.frame(cor=x$estimate, 
		z=x$statistic, 
		pvalue=x$p.value)
	rownames(rst) <- target
	return(rst)
}
```

---
Run corrlation analysis

Correlation between residual expression of miRNA and phenotypes of interest:


```r
rst.cor<-lapply(names(pheno.mir.list), function(x) do.call(rbind, lapply(pheno.mir.list[[x]], 
	function(y) cor.exp(target=y, transc=rfit.misub[,x], pheno=pheno166[,y], alternative="two.sided", method="pearson"))))
names(rst.cor)<-names(pheno.mir.list)
rst.cor
```

```
## $`ssc-miR-190b`
##                 cor         z    pvalue
## last_lum -0.1164132 -1.501022 0.1352727
## 
## $`ssc-miR-184`
##                cor        z     pvalue
## num_ribs 0.2027452 2.258597 0.02572992
## 
## $`ssc-miR-345-3p`
##                    cor         z       pvalue
## juiciness  -0.10657496 -1.372643 1.717375e-01
## driploss   -0.26947966 -3.583595 4.465156e-04
## ph_24h      0.15411099  1.891088 6.058006e-02
## protein     0.29953851  3.983629 1.025856e-04
## cook_yield  0.35347443  4.824296 3.205696e-06
## WBS         0.07946751  1.017792 3.102851e-01
## tenderness -0.17338499 -2.254559 2.548458e-02
## overtend   -0.18569339 -2.420127 1.660749e-02
## 
## $`ssc-let-7d-5p`
##                    cor          z       pvalue
## juiciness   0.10918173  1.4066174 1.614330e-01
## driploss    0.21583473  2.8307543 5.225097e-03
## ph_24h     -0.12386780 -1.5134729 1.323067e-01
## protein    -0.34291206 -4.6319105 7.431531e-06
## cook_yield -0.15606024 -2.0171590 4.532003e-02
## WBS        -0.10119569 -1.2986467 1.958994e-01
## tenderness  0.05964547  0.7651971 4.452540e-01
## overtend    0.06110746  0.7840225 4.341578e-01
## 
## $`ssc-let-7g`
##                    cor          z       pvalue
## juiciness   0.08824985  1.1345762 2.582085e-01
## driploss    0.12866051  1.6614674 9.853003e-02
## ph_24h     -0.02701050 -0.3276044 7.436769e-01
## protein    -0.35028761 -4.7453035 4.565278e-06
## cook_yield -0.10200750 -1.3091737 1.923181e-01
## WBS        -0.05626429 -0.7194741 4.728791e-01
## tenderness  0.03693363  0.4733041 6.366256e-01
## overtend    0.03151038  0.4037302 6.869371e-01
## 
## $`ssc-miR-95`
##                     cor           z       pvalue
## juiciness   0.060228334  0.77270175 4.408111e-01
## driploss    0.097935795  1.26024847 2.093700e-01
## ph_24h     -0.026975174 -0.32717566 7.440005e-01
## protein    -0.314068304 -4.19747039 4.451765e-05
## cook_yield -0.084839197 -1.08707363 2.786092e-01
## WBS         0.010212214  0.13038762 8.964205e-01
## tenderness  0.007743023  0.09916205 9.211307e-01
## overtend   -0.017879477 -0.22900563 8.191499e-01
## 
## $`ssc-miR-1468`
##                   cor         z       pvalue
## juiciness   0.1149193  1.481500 1.403921e-01
## driploss    0.2137734  2.802419 5.683401e-03
## ph_24h     -0.1462412 -1.792349 7.513366e-02
## protein    -0.3333819 -4.486826 1.369823e-05
## cook_yield -0.1517696 -1.960374 5.165673e-02
## WBS        -0.2169999 -2.838096 5.115595e-03
## tenderness  0.1168444  1.506658 1.338220e-01
## overtend    0.1036929  1.335115 1.836891e-01
## 
## $`ssc-miR-874`
##                    cor          z      pvalue
## car_bf10    0.10145177  1.2979671 0.196143598
## WBS        -0.01766144 -0.2255214 0.821856199
## tenderness  0.04057355  0.5200232 0.603748628
## overtend    0.03277630  0.4199671 0.675059198
## lma_16wk    0.24084445  3.1778581 0.001773468
## dress_ptg   0.09489528  1.2207615 0.223928573
## num_ribs    0.06779136  0.7412216 0.460019283
## protein     0.11752885  1.5016813 0.135138155
## ph_24h      0.11591334  1.4149120 0.159209144
## driploss   -0.03396054 -0.4351582 0.664019735
## cook_yield  0.15109938  1.9515139 0.052710017
## juiciness   0.04417062  0.5662126 0.572023276
```

```r
rst.corfil<-lapply(rst.cor, function(x) x[x$pvalue<0.05,])
rst.corfil
```

```
## $`ssc-miR-190b`
## [1] cor    z      pvalue
## <0 rows> (or 0-length row.names)
## 
## $`ssc-miR-184`
##                cor        z     pvalue
## num_ribs 0.2027452 2.258597 0.02572992
## 
## $`ssc-miR-345-3p`
##                   cor         z       pvalue
## driploss   -0.2694797 -3.583595 4.465156e-04
## protein     0.2995385  3.983629 1.025856e-04
## cook_yield  0.3534744  4.824296 3.205696e-06
## tenderness -0.1733850 -2.254559 2.548458e-02
## overtend   -0.1856934 -2.420127 1.660749e-02
## 
## $`ssc-let-7d-5p`
##                   cor         z       pvalue
## driploss    0.2158347  2.830754 5.225097e-03
## protein    -0.3429121 -4.631911 7.431531e-06
## cook_yield -0.1560602 -2.017159 4.532003e-02
## 
## $`ssc-let-7g`
##                cor         z       pvalue
## protein -0.3502876 -4.745303 4.565278e-06
## 
## $`ssc-miR-95`
##                cor        z       pvalue
## protein -0.3140683 -4.19747 4.451765e-05
## 
## $`ssc-miR-1468`
##                 cor         z       pvalue
## driploss  0.2137734  2.802419 5.683401e-03
## protein  -0.3333819 -4.486826 1.369823e-05
## WBS      -0.2169999 -2.838096 5.115595e-03
## 
## $`ssc-miR-874`
##                cor        z      pvalue
## lma_16wk 0.2408444 3.177858 0.001773468
```

Correlation between residual expression of mRNA and phenotypes of interest:


```r
mrst.cor<-lapply(colnames(sub.rMfit), function(x) do.call(rbind, lapply(pheno.mir874, 
	function(y) cor.exp(target=y, transc=sub.rMfit[,x], pheno=pheno166[,y], alternative="two.sided", method="pearson"))))
names(mrst.cor)<-colnames(sub.rMfit)
mrst.cor
```

```
## $XLOC_019568
##                    cor          z       pvalue
## car_bf10    0.18239858  2.3611643 1.940697e-02
## WBS         0.07256556  0.9289040 3.543124e-01
## tenderness  0.07787184  1.0002836 3.186469e-01
## overtend    0.06332997  0.8126507 4.175965e-01
## lma_16wk   -0.14471802 -1.8730122 6.284632e-02
## dress_ptg  -0.18027401 -2.3470874 2.011518e-02
## num_ribs   -0.03367368 -0.3675449 7.138655e-01
## protein     0.10406029  1.3275845 1.861947e-01
## ph_24h      0.37361068  4.8834185 2.689911e-06
## driploss   -0.07060298 -0.9064212 3.660428e-01
## cook_yield  0.21025166  2.7456870 6.716838e-03
## juiciness  -0.02917921 -0.3738354 7.090096e-01
## 
## $XLOC_021684
##                    cor          z       pvalue
## car_bf10    0.14685628  1.8896633 6.058940e-02
## WBS         0.08106082  1.0383323 3.006532e-01
## tenderness  0.11101534  1.4305326 1.544668e-01
## overtend    0.09311102  1.1976056 2.327986e-01
## lma_16wk   -0.14617806 -1.8923193 6.021057e-02
## dress_ptg  -0.18174047 -2.3668295 1.910792e-02
## num_ribs   -0.10671408 -1.1707987 2.440189e-01
## protein     0.12983207  1.6614468 9.856975e-02
## ph_24h      0.34521550  4.4596813 1.619937e-05
## driploss   -0.05123095 -0.6569390 5.121410e-01
## cook_yield  0.23167012  3.0404843 2.752779e-03
## juiciness  -0.02234943 -0.2862839 7.750220e-01
## 
## $XLOC_022612
##                    cor           z     pvalue
## car_bf10    0.02802124  0.35679228 0.72171162
## WBS         0.05766443  0.73743720 0.46191674
## tenderness  0.01197977  0.15342698 0.87825019
## overtend    0.03880571  0.49733016 0.61962236
## lma_16wk   -0.12591741 -1.62546717 0.10598321
## dress_ptg  -0.08307652 -1.06758900 0.28727494
## num_ribs    0.04789099  0.52302919 0.60192645
## protein     0.16838650  2.16753533 0.03166322
## ph_24h      0.03760747  0.45628910 0.64885552
## driploss    0.00436737  0.05593016 0.95546552
## cook_yield  0.08303710  1.06382071 0.28898336
## juiciness  -0.11728875 -1.51246813 0.13233944
## 
## $XLOC_022617
##                     cor          z       pvalue
## car_bf10    0.154570433  1.9912922 4.812985e-02
## WBS         0.140212122  1.8079686 7.245491e-02
## tenderness  0.093757389  1.2059927 2.295573e-01
## overtend    0.086666312  1.1140621 2.668826e-01
## lma_16wk   -0.169041442 -2.1963951 2.946698e-02
## dress_ptg  -0.107112954 -1.3796524 1.695719e-01
## num_ribs    0.097905849  1.0731826 2.853606e-01
## protein     0.148936547  1.9111080 5.776758e-02
## ph_24h      0.404391586  5.3608821 3.139977e-07
## driploss   -0.020808962 -0.2665424 7.901563e-01
## cook_yield  0.232878988  3.0572570 2.611260e-03
## juiciness  -0.009751106 -0.1248810 9.007707e-01
## 
## $XLOC_003117
##                    cor          z       pvalue
## car_bf10    0.14764587  1.9000491 5.920326e-02
## WBS         0.09219514  1.1821034 2.388863e-01
## tenderness  0.11623004  1.4986280 1.358926e-01
## overtend    0.09785647  1.2592179 2.097409e-01
## lma_16wk   -0.13062976 -1.6873355 9.343988e-02
## dress_ptg  -0.15923970 -2.0656206 4.043749e-02
## num_ribs   -0.14924496 -1.6465108 1.022971e-01
## protein     0.08469177  1.0784929 2.824270e-01
## ph_24h      0.39320384  5.1849871 7.029873e-07
## driploss   -0.04522359 -0.5797376 5.628871e-01
## cook_yield  0.19652715  2.5589953 1.140730e-02
## juiciness   0.01299308  0.1664066 8.680420e-01
## 
## $XLOC_008680
##                     cor           z      pvalue
## car_bf10   -0.054610968 -0.69612297 0.487349345
## WBS         0.047908563  0.61235874 0.541153742
## tenderness -0.014332131 -0.18355969 0.854585640
## overtend   -0.019014334 -0.24354631 0.807886593
## lma_16wk    0.007036497  0.09011337 0.928307089
## dress_ptg  -0.192824331 -2.51658430 0.012810338
## num_ribs   -0.235449809 -2.64275114 0.009329981
## protein     0.031607376  0.40125312 0.688765833
## ph_24h     -0.167459707 -2.05942231 0.041218115
## driploss   -0.013342763 -0.17088595 0.864524091
## cook_yield -0.204837326 -2.67184156 0.008309573
## juiciness  -0.006201722 -0.07942232 0.936793545
## 
## $XLOC_008685
##                     cor           z     pvalue
## car_bf10    0.094746270  1.21137256 0.22751673
## WBS         0.047488945  0.60698309 0.54470669
## tenderness  0.060247977  0.77295469 0.44066180
## overtend    0.072431440  0.93001781 0.35372878
## lma_16wk   -0.154614951 -2.00413763 0.04670135
## dress_ptg  -0.094933930 -1.22126324 0.22373911
## num_ribs   -0.076500466 -0.83697427 0.40428472
## protein    -0.063783406 -0.81097202 0.41857895
## ph_24h      0.172058865  2.11768460 0.03588462
## driploss    0.065685605  0.84300677 0.40045247
## cook_yield  0.076696809  0.98209210 0.32750993
## juiciness  -0.006014295 -0.07702195 0.93869996
## 
## $XLOC_009274
##                    cor          z       pvalue
## car_bf10    0.19252375  2.4971431 1.351889e-02
## WBS         0.08485130  1.0872298 2.785403e-01
## tenderness  0.10178929  1.3103449 1.919114e-01
## overtend    0.08313671  1.0683679 2.869246e-01
## lma_16wk   -0.13979504 -1.8080038 7.243814e-02
## dress_ptg  -0.15400545 -1.9960448 4.758424e-02
## num_ribs   -0.13485089 -1.4846101 1.402913e-01
## protein     0.11449819  1.4624369 1.455708e-01
## ph_24h      0.41039725  5.4564808 2.012276e-07
## driploss   -0.01746421 -0.2236851 8.232807e-01
## cook_yield  0.24213429  3.1861756 1.728010e-03
## juiciness  -0.04479207 -0.5741947 5.666228e-01
## 
## $XLOC_009349
##                    cor          z     pvalue
## car_bf10    0.13777710  1.7705010 0.07852451
## WBS         0.02302001  0.2939777 0.76914896
## tenderness  0.13416549  1.7338323 0.08482710
## overtend    0.13976375  1.8075911 0.07250271
## lma_16wk   -0.04693223 -0.6016888 0.54821215
## dress_ptg  -0.03687177 -0.4725104 0.63719069
## num_ribs   -0.02389002 -0.2606838 0.79478694
## protein     0.02010601  0.2551683 0.79891879
## ph_24h      0.14577862  1.7865572 0.07607014
## driploss   -0.01432054 -0.1834113 0.85470190
## cook_yield  0.18691306  2.4291565 0.01622059
## juiciness   0.03317422  0.4250713 0.67134193
## 
## $XLOC_011766
##                    cor          z       pvalue
## car_bf10    0.14779193  1.9019707 5.894972e-02
## WBS         0.10001953  1.2833995 2.011737e-01
## tenderness  0.09554089  1.2291430 2.207788e-01
## overtend    0.08397731  1.0792466 2.820625e-01
## lma_16wk   -0.12353394 -1.5942175 1.128127e-01
## dress_ptg  -0.13242200 -1.7108962 8.899104e-02
## num_ribs   -0.08079213 -0.8842287 3.783558e-01
## protein     0.12528114  1.6022632 1.110576e-01
## ph_24h      0.36851527  4.8062671 3.760057e-06
## driploss   -0.02048652 -0.2624106 7.933342e-01
## cook_yield  0.20976145  2.7389903 6.848924e-03
## juiciness  -0.05243267 -0.6723907 5.022812e-01
## 
## $XLOC_021695
##                    cor          z       pvalue
## car_bf10    0.17587453  2.2739624 2.428097e-02
## WBS         0.11069626  1.4220145 1.569327e-01
## tenderness  0.09976434  1.2840127 2.009485e-01
## overtend    0.08758969  1.1260230 2.618008e-01
## lma_16wk   -0.15622566 -2.0255354 4.443373e-02
## dress_ptg  -0.15401181 -1.9961294 4.757494e-02
## num_ribs    0.02231666  0.2435067 8.080322e-01
## protein     0.13672906  1.7513451 8.179081e-02
## ph_24h      0.39623690  5.2323977 5.666549e-07
## driploss   -0.01089780 -0.1395682 8.891724e-01
## cook_yield  0.23068891  3.0268812 2.872710e-03
## juiciness  -0.02892981 -0.3706375 7.113858e-01
## 
## $XLOC_022870
##                     cor           z       pvalue
## car_bf10    0.154223373  1.98671209 4.864055e-02
## WBS         0.129103036  1.66218777 9.839702e-02
## tenderness  0.082507250  1.06022322 2.906020e-01
## overtend    0.077456840  0.99492057 3.212406e-01
## lma_16wk   -0.172850075 -2.24738837 2.594867e-02
## dress_ptg  -0.183508254 -2.39064992 1.795202e-02
## num_ribs    0.056942083  0.62217427 5.350176e-01
## protein     0.166335079  2.14037246 3.383044e-02
## ph_24h      0.337650950  4.34922408 2.540213e-05
## driploss   -0.004579278 -0.05864398 9.533071e-01
## cook_yield  0.193191934  2.51386831 1.291164e-02
## juiciness  -0.022918823 -0.29358125 7.694491e-01
## 
## $XLOC_022875
##                    cor          z       pvalue
## car_bf10    0.24738945  3.2497690 1.404666e-03
## WBS         0.10404130  1.3355585 1.835556e-01
## tenderness  0.08407754  1.0805438 2.814865e-01
## overtend    0.06725465  0.8632342 3.892694e-01
## lma_16wk   -0.19649215 -2.5663575 1.117178e-02
## dress_ptg  -0.18538109 -2.4159120 1.679411e-02
## num_ribs    0.14865395  1.6398430 1.036787e-01
## protein     0.13848519  1.7742761 7.790750e-02
## ph_24h      0.34552411  4.4642084 1.590083e-05
## driploss   -0.06720169 -0.8625515 3.896437e-01
## cook_yield  0.22847698  2.9962517 3.160580e-03
## juiciness  -0.02619177 -0.3355335 7.376515e-01
## 
## $XLOC_023161
##                     cor           z       pvalue
## car_bf10    0.199783912  2.59515238 1.032283e-02
## WBS         0.122705734  1.57853074 1.163824e-01
## tenderness  0.099279983  1.27771664 2.031550e-01
## overtend    0.089949090  1.15659882 2.491183e-01
## lma_16wk   -0.174060075 -2.26361050 2.490906e-02
## dress_ptg  -0.167216840 -2.17200187 3.129169e-02
## num_ribs    0.119381005  1.31167343 1.921553e-01
## protein     0.142968273  1.83289283 6.866561e-02
## ph_24h      0.390677653  5.14565242 8.398867e-07
## driploss   -0.001772876 -0.02270393 9.819140e-01
## cook_yield  0.224014262  2.93460303 3.821986e-03
## juiciness  -0.017835804 -0.22844608 8.195841e-01
## 
## $XLOC_008918
##                    cor          z       pvalue
## car_bf10    0.16789727  2.1677558 3.163701e-02
## WBS         0.08197376  1.0501050 2.952242e-01
## tenderness  0.11400117  1.4695076 1.436109e-01
## overtend    0.09305193  1.1968389 2.330966e-01
## lma_16wk   -0.13758869 -1.7789133 7.710639e-02
## dress_ptg  -0.13603781 -1.7584815 8.053088e-02
## num_ribs   -0.08900760 -0.9748274 3.316234e-01
## protein     0.11390342  1.4547400 1.476881e-01
## ph_24h      0.39659797  5.2380552 5.522179e-07
## driploss   -0.03087654 -0.3956012 6.929134e-01
## cook_yield  0.23566871  3.0960209 2.309495e-03
## juiciness  -0.03427187 -0.4391521 6.611294e-01
## 
## $XLOC_023143
##                    cor          z       pvalue
## car_bf10    0.16940157  2.1877492 3.012045e-02
## WBS         0.11261287  1.4469490 1.498316e-01
## tenderness  0.07702093  0.9892879 3.239797e-01
## overtend    0.06154314  0.7896336 4.308819e-01
## lma_16wk   -0.16472238 -2.1386903 3.394104e-02
## dress_ptg  -0.15880694 -2.0598615 4.099201e-02
## num_ribs   -0.10155062 -1.1135430 2.677195e-01
## protein     0.12661178  1.6195570 1.072846e-01
## ph_24h      0.36077912  4.6900856 6.185152e-06
## driploss   -0.06901370 -0.8859189 3.769584e-01
## cook_yield  0.23600303  3.1006718 2.275540e-03
## juiciness  -0.05567653 -0.7141152 4.761711e-01
## 
## $XLOC_021691
##                    cor          z       pvalue
## car_bf10    0.18337693  2.3742685 1.875451e-02
## WBS         0.14094826  1.8176526 7.095353e-02
## tenderness  0.07602096  0.9763688 3.303198e-01
## overtend    0.07221313  0.9272000 3.551852e-01
## lma_16wk   -0.14387828 -1.8619135 6.440459e-02
## dress_ptg  -0.15167689 -1.9651484 5.108621e-02
## num_ribs    0.11254280  1.2355465 2.190608e-01
## protein     0.16487790  2.1210954 3.544510e-02
## ph_24h      0.37722418  4.9384425 2.113872e-06
## driploss   -0.03745991 -0.4800579 6.318259e-01
## cook_yield  0.22771028  2.9856463 3.266285e-03
## juiciness  -0.04628735 -0.5934034 5.537288e-01
## 
## $XLOC_002346
##                    cor          z     pvalue
## car_bf10    0.06538272  0.8339707 0.40552552
## WBS         0.02529053  0.3229911 0.74711592
## tenderness  0.06969816  0.8947479 0.37223315
## overtend    0.06412543  0.8228998 0.41175989
## lma_16wk   -0.03370101 -0.4318288 0.66643303
## dress_ptg  -0.05960049 -0.7646179 0.44559791
## num_ribs   -0.02810254 -0.3066836 0.75962016
## protein     0.07219274  0.9184197 0.35977304
## ph_24h      0.16494734  2.0276543 0.04440224
## driploss   -0.05429969 -0.6964028 0.48716246
## cook_yield  0.17052333  2.2094566 0.02853700
## juiciness  -0.05738167 -0.7360567 0.46274765
## 
## $XLOC_011794
##                     cor           z      pvalue
## car_bf10    0.064597260  0.82390970 0.411202207
## WBS         0.122371071  1.57415997 0.117389440
## tenderness  0.005290738  0.06775545 0.946062842
## overtend    0.001091308  0.01397557 0.988866455
## lma_16wk   -0.186677507 -2.43341494 0.016031143
## dress_ptg  -0.185659352 -2.41966772 0.016627740
## num_ribs   -0.059024546 -0.64500633 0.520164660
## protein     0.025006107  0.31739117 0.751358041
## ph_24h      0.220412164  2.73973423 0.006910705
## driploss    0.037681364  0.48289986 0.629810836
## cook_yield  0.038431640  0.49102509 0.624069151
## juiciness  -0.010514827 -0.13466293 0.893043506
## 
## $XLOC_013030
##                     cor           z       pvalue
## car_bf10    0.061646117  0.78612212 4.329443e-01
## WBS        -0.003534103 -0.04512069 9.640663e-01
## tenderness  0.125329189  1.61775237 1.076377e-01
## overtend    0.091373734  1.17507045 2.416698e-01
## lma_16wk   -0.146002338 -1.88999495 6.052291e-02
## dress_ptg  -0.087438556 -1.12406514 2.626280e-01
## num_ribs   -0.116647215 -1.28121723 2.026085e-01
## protein     0.088148695  1.12285245 2.631718e-01
## ph_24h      0.330382006  4.24397986 3.871471e-05
## driploss   -0.036050002 -0.46196556 6.447182e-01
## cook_yield  0.159177461  2.05848752 4.113500e-02
## juiciness   0.057316184  0.73521392 4.632593e-01
## 
## $XLOC_009607
##                    cor          z       pvalue
## car_bf10    0.19685070  2.5555026 1.152360e-02
## WBS         0.11963465  1.5384421 1.258795e-01
## tenderness  0.06954990  0.8928353 3.732536e-01
## overtend    0.05799619  0.7439658 4.579617e-01
## lma_16wk   -0.11121838 -1.4331817 1.537095e-01
## dress_ptg  -0.14324456 -1.8535403 6.560138e-02
## num_ribs    0.03837238  0.4189018 6.760430e-01
## protein     0.08289601  1.0554651 2.927939e-01
## ph_24h      0.40425120  5.3586575 3.172476e-07
## driploss   -0.02829753 -0.3625303 7.174225e-01
## cook_yield  0.23218978  3.0476926 2.691114e-03
## juiciness  -0.06870355 -0.8819186 3.791116e-01
## 
## $XLOC_022575
##                     cor          z       pvalue
## car_bf10    0.207354911  2.6978324 7.718382e-03
## WBS         0.107871477  1.3852942 1.678553e-01
## tenderness  0.115667888  1.4912813 1.378086e-01
## overtend    0.107876061  1.3895968 1.665349e-01
## lma_16wk   -0.176841758 -2.3009440 2.265376e-02
## dress_ptg  -0.147835246 -1.9142487 5.732966e-02
## num_ribs    0.114538209  1.2577417 2.109477e-01
## protein     0.123871544  1.5839529 1.151668e-01
## ph_24h      0.348076567  4.5017138 1.362314e-05
## driploss    0.021992580  0.2817106 7.785206e-01
## cook_yield  0.212101216  2.7709730 6.238672e-03
## juiciness  -0.005363992 -0.0686936 9.453172e-01
## 
## $XLOC_023083
##                     cor          z       pvalue
## car_bf10    0.138678044  1.7823048 0.0765725633
## WBS         0.150976758  1.9498932 0.0529046198
## tenderness  0.039412077  0.5051133 0.6141573715
## overtend    0.038073019  0.4879263 0.6262537195
## lma_16wk   -0.179617458 -2.3382539 0.0205808231
## dress_ptg  -0.176579741 -2.2974250 0.0228584154
## num_ribs    0.157792542  1.7431511 0.0838901285
## protein     0.130310804  1.6676787 0.0973234044
## ph_24h      0.289575562  3.6680752 0.0003409604
## driploss    0.001145346  0.0146676 0.9883151992
## cook_yield  0.171433462  2.2216052 0.0276873952
## juiciness  -0.037009817 -0.4742818 0.6359298051
## 
## $XLOC_009276
##                      cor           z     pvalue
## car_bf10    0.1879782310  2.43599834 0.01593455
## WBS        -0.0009939919 -0.01269044 0.98989028
## tenderness -0.0747158453 -0.95951164 0.33871351
## overtend   -0.0678234878 -0.87056907 0.38526200
## lma_16wk   -0.1375546740 -1.77846511 0.07718021
## dress_ptg  -0.1320363426 -1.70582492 0.08993379
## num_ribs    0.0798773691  0.87415240 0.38379600
## protein    -0.0731203511 -0.93028350 0.35361709
## ph_24h      0.0468059472  0.56811460 0.57082392
## driploss    0.1125171543  1.45013128 0.14893169
## cook_yield -0.0443961227 -0.56737118 0.57124261
## juiciness  -0.0335973417 -0.43049894 0.66739795
## 
## $XLOC_002373
##                    cor          z       pvalue
## car_bf10    0.17223688  2.2254761 2.742988e-02
## WBS         0.09124641  1.1698364 2.437745e-01
## tenderness  0.10339674  1.3312597 1.849512e-01
## overtend    0.08303009  1.0669882 2.875453e-01
## lma_16wk   -0.13486610 -1.7430537 8.319846e-02
## dress_ptg  -0.15413197 -1.9977245 4.739983e-02
## num_ribs   -0.11809888 -1.2973860 1.970079e-01
## protein     0.13200585  1.6897536 9.301054e-02
## ph_24h      0.40005942  5.2924385 4.305126e-07
## driploss   -0.05237716 -0.6716768 5.027345e-01
## cook_yield  0.24825286  3.2719066 1.303963e-03
## juiciness  -0.04619424 -0.5922072 5.545275e-01
## 
## $XLOC_008850
##                    cor          z       pvalue
## car_bf10    0.18618357  2.4119021 1.698716e-02
## WBS         0.08477981  1.0863071 2.789470e-01
## tenderness  0.10793649  1.3903844 1.662962e-01
## overtend    0.09458491  1.2167327 2.254540e-01
## lma_16wk   -0.14727395 -1.9068192 5.829245e-02
## dress_ptg  -0.14453369 -1.8705756 6.318571e-02
## num_ribs   -0.11427326 -1.2547938 2.120124e-01
## protein     0.10179874  1.2984265 1.959976e-01
## ph_24h      0.39326537  5.1859468 6.999345e-07
## driploss   -0.02537062 -0.3250071 7.455899e-01
## cook_yield  0.24279947  3.1954760 1.676478e-03
## juiciness  -0.04419880 -0.5665745 5.717778e-01
## 
## $XLOC_009284
##                     cor           z       pvalue
## car_bf10    0.198617439  2.57937594 1.078645e-02
## WBS         0.065407517  0.83685930 4.038969e-01
## tenderness  0.092253081  1.18647549 2.371505e-01
## overtend    0.087799135  1.12873649 2.606574e-01
## lma_16wk   -0.190440396 -2.48429266 1.398505e-02
## dress_ptg  -0.143766019 -1.86043006 6.461529e-02
## num_ribs   -0.117720072 -1.29316599 1.984585e-01
## protein     0.072854088  0.92687781 3.553773e-01
## ph_24h      0.393442526  5.18871071 6.912150e-07
## driploss   -0.003616533 -0.04631453 9.631159e-01
## cook_yield  0.220082036  2.88044388 4.505351e-03
## juiciness  -0.027516384 -0.35251513 7.249046e-01
## 
## $XLOC_009537
##                    cor          z       pvalue
## car_bf10    0.01724843  0.2195693 8.264830e-01
## WBS        -0.04363238 -0.5575919 5.778880e-01
## tenderness  0.11754032  1.5157575 1.315057e-01
## overtend    0.11668923  1.5046301 1.343426e-01
## lma_16wk   -0.19606587 -2.5605670 1.135223e-02
## dress_ptg  -0.06464365 -0.8295777 4.079835e-01
## num_ribs   -0.12531213 -1.3778551 1.708336e-01
## protein     0.09806420  1.2503216 2.129963e-01
## ph_24h      0.21375889  2.6530093 8.854736e-03
## driploss   -0.05927842 -0.7604715 4.480647e-01
## cook_yield  0.30802918  4.1336439 5.699096e-05
## juiciness   0.09847188  1.2672142 2.068752e-01
## 
## $XLOC_011767
##                    cor          z       pvalue
## car_bf10    0.18943088  2.4555209 1.512461e-02
## WBS         0.07619639  0.9756467 3.306853e-01
## tenderness  0.12401978  1.6005850 1.113934e-01
## overtend    0.10817793  1.3935313 1.653448e-01
## lma_16wk   -0.12818045 -1.6551644 9.980355e-02
## dress_ptg  -0.15349895 -1.9893215 4.832845e-02
## num_ribs   -0.07160515 -0.7831302 4.351064e-01
## protein     0.12568377  1.6074951 1.099052e-01
## ph_24h      0.38800572  5.1041977 1.012177e-06
## driploss   -0.04255873 -0.5455119 5.861427e-01
## cook_yield  0.22752105  2.9830297 3.292860e-03
## juiciness  -0.02449085 -0.3137300 7.541248e-01
## 
## $XLOC_013026
##                    cor          z       pvalue
## car_bf10    0.15818164  2.0389945 0.0430758629
## WBS         0.04367661  0.5581582 0.5775021916
## tenderness  0.08332855  1.0708503 0.2858101417
## overtend    0.05874008  0.7535412 0.4522052308
## lma_16wk   -0.06517019 -0.8363636 0.4041673402
## dress_ptg  -0.17256289 -2.2435397 0.0262008045
## num_ribs    0.03057799  0.3337225 0.7391762751
## protein     0.10765775  1.3740095 0.1713486101
## ph_24h      0.29417774  3.7318471 0.0002709771
## driploss   -0.05545276 -0.7112362 0.4779482635
## cook_yield  0.20487917  2.6724113 0.0082960793
## juiciness  -0.04645407 -0.5955454 0.5522999391
## 
## $XLOC_013029
##                      cor            z     pvalue
## car_bf10    0.0307119687  0.391084027 0.69624899
## WBS         0.0607212555  0.776670235 0.43847896
## tenderness  0.0668589314  0.858132212 0.39207186
## overtend    0.0899180985  1.156197077 0.24928206
## lma_16wk   -0.1026845588 -1.321992065 0.18801174
## dress_ptg  -0.0030612240 -0.039202979 0.96877624
## num_ribs    0.1693568298  1.874543050 0.06330745
## protein     0.0002780801  0.003528441 0.99718909
## ph_24h      0.1690279310  2.079272864 0.03932934
## driploss    0.0174647211  0.223691675 0.82327563
## cook_yield  0.0550224950  0.703545982 0.48271932
## juiciness  -0.0020208726 -0.025879849 0.97938464
```

```r
mrst.corfil<-lapply(mrst.cor, function(x) x[x$pvalue<0.05,])
names(mrst.corfil)<-paste(names(mrst.corfil),as.character(negcoloc[match(names(mrst.corfil), negcoloc$ID),"genes"]), sep=".")
mrst.corfil
```

```
## $XLOC_019568.GPR157
##                   cor         z       pvalue
## car_bf10    0.1823986  2.361164 1.940697e-02
## dress_ptg  -0.1802740 -2.347087 2.011518e-02
## ph_24h      0.3736107  4.883419 2.689911e-06
## cook_yield  0.2102517  2.745687 6.716838e-03
## 
## $XLOC_021684.PFDN6
##                   cor         z       pvalue
## dress_ptg  -0.1817405 -2.366830 1.910792e-02
## ph_24h      0.3452155  4.459681 1.619937e-05
## cook_yield  0.2316701  3.040484 2.752779e-03
## 
## $XLOC_022612.ETV7
##               cor        z     pvalue
## protein 0.1683865 2.167535 0.03166322
## 
## $XLOC_022617.CDKN1A
##                   cor         z       pvalue
## car_bf10    0.1545704  1.991292 4.812985e-02
## lma_16wk   -0.1690414 -2.196395 2.946698e-02
## ph_24h      0.4043916  5.360882 3.139977e-07
## cook_yield  0.2328790  3.057257 2.611260e-03
## 
## $XLOC_003117.ATP8A2
##                   cor         z       pvalue
## dress_ptg  -0.1592397 -2.065621 4.043749e-02
## ph_24h      0.3932038  5.184987 7.029873e-07
## cook_yield  0.1965272  2.558995 1.140730e-02
## 
## $XLOC_008680.PLEKHB2
##                   cor         z      pvalue
## dress_ptg  -0.1928243 -2.516584 0.012810338
## num_ribs   -0.2354498 -2.642751 0.009329981
## ph_24h     -0.1674597 -2.059422 0.041218115
## cook_yield -0.2048373 -2.671842 0.008309573
## 
## $XLOC_008685.WDR33
##                 cor         z     pvalue
## lma_16wk -0.1546150 -2.004138 0.04670135
## ph_24h    0.1720589  2.117685 0.03588462
## 
## $XLOC_009274.PLEKHB2
##                   cor         z       pvalue
## car_bf10    0.1925238  2.497143 1.351889e-02
## dress_ptg  -0.1540054 -1.996045 4.758424e-02
## ph_24h      0.4103972  5.456481 2.012276e-07
## cook_yield  0.2421343  3.186176 1.728010e-03
## 
## $XLOC_009349.METTL8
##                  cor        z     pvalue
## cook_yield 0.1869131 2.429157 0.01622059
## 
## $XLOC_011766.CARNS1
##                  cor        z       pvalue
## ph_24h     0.3685153 4.806267 3.760057e-06
## cook_yield 0.2097614 2.738990 6.848924e-03
## 
## $XLOC_021695.NUDT3
##                   cor         z       pvalue
## car_bf10    0.1758745  2.273962 2.428097e-02
## lma_16wk   -0.1562257 -2.025535 4.443373e-02
## dress_ptg  -0.1540118 -1.996129 4.757494e-02
## ph_24h      0.3962369  5.232398 5.666549e-07
## cook_yield  0.2306889  3.026881 2.872710e-03
## 
## $XLOC_022870.KHNYN
##                   cor         z       pvalue
## car_bf10    0.1542234  1.986712 4.864055e-02
## lma_16wk   -0.1728501 -2.247388 2.594867e-02
## dress_ptg  -0.1835083 -2.390650 1.795202e-02
## protein     0.1663351  2.140372 3.383044e-02
## ph_24h      0.3376510  4.349224 2.540213e-05
## cook_yield  0.1931919  2.513868 1.291164e-02
## 
## $XLOC_022875.NOP9
##                   cor         z       pvalue
## car_bf10    0.2473894  3.249769 1.404666e-03
## lma_16wk   -0.1964922 -2.566357 1.117178e-02
## dress_ptg  -0.1853811 -2.415912 1.679411e-02
## ph_24h      0.3455241  4.464208 1.590083e-05
## cook_yield  0.2284770  2.996252 3.160580e-03
## 
## $XLOC_023161.POMT2
##                   cor         z       pvalue
## car_bf10    0.1997839  2.595152 1.032283e-02
## lma_16wk   -0.1740601 -2.263611 2.490906e-02
## dress_ptg  -0.1672168 -2.172002 3.129169e-02
## ph_24h      0.3906777  5.145652 8.398867e-07
## cook_yield  0.2240143  2.934603 3.821986e-03
## 
## $XLOC_008918.IGFBP5
##                  cor        z       pvalue
## car_bf10   0.1678973 2.167756 3.163701e-02
## ph_24h     0.3965980 5.238055 5.522179e-07
## cook_yield 0.2356687 3.096021 2.309495e-03
## 
## $XLOC_023143.PROX2
##                   cor         z       pvalue
## car_bf10    0.1694016  2.187749 3.012045e-02
## lma_16wk   -0.1647224 -2.138690 3.394104e-02
## dress_ptg  -0.1588069 -2.059862 4.099201e-02
## ph_24h      0.3607791  4.690086 6.185152e-06
## cook_yield  0.2360030  3.100672 2.275540e-03
## 
## $XLOC_021691.LEMD2
##                  cor        z       pvalue
## car_bf10   0.1833769 2.374268 1.875451e-02
## protein    0.1648779 2.121095 3.544510e-02
## ph_24h     0.3772242 4.938443 2.113872e-06
## cook_yield 0.2277103 2.985646 3.266285e-03
## 
## $XLOC_002346.UCK1
##                  cor        z     pvalue
## ph_24h     0.1649473 2.027654 0.04440224
## cook_yield 0.1705233 2.209457 0.02853700
## 
## $XLOC_011794.OVOL1
##                  cor         z      pvalue
## lma_16wk  -0.1866775 -2.433415 0.016031143
## dress_ptg -0.1856594 -2.419668 0.016627740
## ph_24h     0.2204122  2.739734 0.006910705
## 
## $XLOC_013030.OVOL1
##                  cor        z       pvalue
## ph_24h     0.3303820 4.243980 3.871471e-05
## cook_yield 0.1591775 2.058488 4.113500e-02
## 
## $XLOC_009607.KLHL30
##                  cor        z       pvalue
## car_bf10   0.1968507 2.555503 1.152360e-02
## ph_24h     0.4042512 5.358658 3.172476e-07
## cook_yield 0.2321898 3.047693 2.691114e-03
## 
## $XLOC_022575.ZNF451
##                   cor         z       pvalue
## car_bf10    0.2073549  2.697832 7.718382e-03
## lma_16wk   -0.1768418 -2.300944 2.265376e-02
## ph_24h      0.3480766  4.501714 1.362314e-05
## cook_yield  0.2121012  2.770973 6.238672e-03
## 
## $XLOC_023083.MAX
##                   cor         z       pvalue
## lma_16wk   -0.1796175 -2.338254 0.0205808231
## dress_ptg  -0.1765797 -2.297425 0.0228584154
## ph_24h      0.2895756  3.668075 0.0003409604
## cook_yield  0.1714335  2.221605 0.0276873952
## 
## $XLOC_009276.UGGT1
##                cor        z     pvalue
## car_bf10 0.1879782 2.435998 0.01593455
## 
## $XLOC_002373.RXRA
##                   cor         z       pvalue
## car_bf10    0.1722369  2.225476 2.742988e-02
## dress_ptg  -0.1541320 -1.997725 4.739983e-02
## ph_24h      0.4000594  5.292439 4.305126e-07
## cook_yield  0.2482529  3.271907 1.303963e-03
## 
## $XLOC_008850.NIF3L1
##                  cor        z       pvalue
## car_bf10   0.1861836 2.411902 1.698716e-02
## ph_24h     0.3932654 5.185947 6.999345e-07
## cook_yield 0.2427995 3.195476 1.676478e-03
## 
## $XLOC_009284.PRPF40A
##                   cor         z       pvalue
## car_bf10    0.1986174  2.579376 1.078645e-02
## lma_16wk   -0.1904404 -2.484293 1.398505e-02
## ph_24h      0.3934425  5.188711 6.912150e-07
## cook_yield  0.2200820  2.880444 4.505351e-03
## 
## $XLOC_009537.PTPRN
##                   cor         z       pvalue
## lma_16wk   -0.1960659 -2.560567 1.135223e-02
## ph_24h      0.2137589  2.653009 8.854736e-03
## cook_yield  0.3080292  4.133644 5.699096e-05
## 
## $XLOC_011767.PPP1CA
##                   cor         z       pvalue
## car_bf10    0.1894309  2.455521 1.512461e-02
## dress_ptg  -0.1534990 -1.989322 4.832845e-02
## ph_24h      0.3880057  5.104198 1.012177e-06
## cook_yield  0.2275210  2.983030 3.292860e-03
## 
## $XLOC_013026.CCDC85B
##                   cor         z       pvalue
## car_bf10    0.1581816  2.038995 0.0430758629
## dress_ptg  -0.1725629 -2.243540 0.0262008045
## ph_24h      0.2941777  3.731847 0.0002709771
## cook_yield  0.2048792  2.672411 0.0082960793
## 
## $XLOC_013029.SNX32
##              cor        z     pvalue
## ph_24h 0.1690279 2.079273 0.03932934
```

Summarize results; split negcoloc by XLOC ID in mrst.cor:


```r
splt.negcoloc<-lapply(names(mrst.cor), function(x) negcoloc[negcoloc$ID==x,])
names(splt.negcoloc)<-names(mrst.cor)
splt.negcoloc
```

```
## $XLOC_019568
##    chr          ID  genes       miRNA        cor    pheno cor.mag
## 35   6 XLOC_019568 GPR157 ssc-miR-874 -0.1602775 lma_16wk       -
## 
## $XLOC_021684
##    chr          ID genes       miRNA        cor     pheno cor.mag
## 39   7 XLOC_021684 PFDN6 ssc-miR-874 -0.1512231 dress_ptg       -
## 
## $XLOC_022612
##    chr          ID genes       miRNA        cor     pheno cor.mag
## 44   7 XLOC_022612  ETV7 ssc-miR-874 -0.1256663 dress_ptg       -
## 
## $XLOC_022617
##    chr          ID  genes       miRNA        cor     pheno cor.mag
## 45   7 XLOC_022617 CDKN1A ssc-miR-874 -0.0875502 dress_ptg       -
## 
## $XLOC_003117
##    chr          ID  genes       miRNA        cor     pheno cor.mag
## 58  11 XLOC_003117 ATP8A2 ssc-miR-874 -0.1738591 dress_ptg       -
## 
## $XLOC_008680
##    chr          ID   genes       miRNA        cor   pheno cor.mag
## 60  15 XLOC_008680 PLEKHB2 ssc-miR-874 -0.1192406 protein       -
## 
## $XLOC_008685
##    chr          ID genes       miRNA       cor   pheno cor.mag
## 62  15 XLOC_008685 WDR33 ssc-miR-874 -0.109456 protein       -
## 
## $XLOC_009274
##    chr          ID   genes       miRNA        cor   pheno cor.mag
## 68  15 XLOC_009274 PLEKHB2 ssc-miR-874 -0.1592552 protein       -
## 
## $XLOC_009349
##    chr          ID  genes       miRNA        cor   pheno cor.mag
## 71  15 XLOC_009349 METTL8 ssc-miR-874 -0.1030303 protein       -
## 
## $XLOC_011766
##    chr          ID  genes       miRNA        cor      pheno cor.mag
## 8    2 XLOC_011766 CARNS1 ssc-miR-874 -0.1198248        WBS       -
## 18   2 XLOC_011766 CARNS1 ssc-miR-874 -0.1198248 tenderness       -
## 
## $XLOC_021695
##    chr          ID genes       miRNA       cor     pheno cor.mag
## 41   7 XLOC_021695 NUDT3 ssc-miR-874 -0.132092 dress_ptg       -
## 
## $XLOC_022870
##    chr          ID genes       miRNA        cor    pheno cor.mag
## 50   7 XLOC_022870 KHNYN ssc-miR-874 -0.1599854 num_ribs       -
## 
## $XLOC_022875
##    chr          ID genes       miRNA        cor    pheno cor.mag
## 51   7 XLOC_022875  NOP9 ssc-miR-874 -0.1525374 num_ribs       -
## 
## $XLOC_023161
##    chr          ID genes       miRNA        cor    pheno cor.mag
## 56   7 XLOC_023161 POMT2 ssc-miR-874 -0.1209931 num_ribs       -
## 
## $XLOC_008918
##    chr          ID  genes       miRNA        cor      pheno cor.mag
## 65  15 XLOC_008918 IGFBP5 ssc-miR-874 -0.1161738    protein       -
## 77  15 XLOC_008918 IGFBP5 ssc-miR-874 -0.1161738     ph_24h       -
## 81  15 XLOC_008918 IGFBP5 ssc-miR-874 -0.1161738   driploss       -
## 88  15 XLOC_008918 IGFBP5 ssc-miR-874 -0.1161738 cook_yield       -
## 92  15 XLOC_008918 IGFBP5 ssc-miR-874 -0.1161738   overtend       -
## 
## $XLOC_023143
##    chr          ID genes       miRNA        cor    pheno cor.mag
## 54   7 XLOC_023143 PROX2 ssc-miR-874 -0.1107704 num_ribs       -
## 
## $XLOC_021691
##    chr          ID genes       miRNA        cor     pheno cor.mag
## 40   7 XLOC_021691 LEMD2 ssc-miR-874 -0.1481563 dress_ptg       -
## 
## $XLOC_002346
##   chr          ID genes       miRNA        cor    pheno cor.mag
## 4   1 XLOC_002346  UCK1 ssc-miR-874 -0.1034684 car_bf10       -
## 
## $XLOC_011794
##    chr          ID genes       miRNA       cor      pheno cor.mag
## 11   2 XLOC_011794 OVOL1 ssc-miR-874 -0.154582        WBS       -
## 21   2 XLOC_011794 OVOL1 ssc-miR-874 -0.154582 tenderness       -
## 28   2 XLOC_011794 OVOL1 ssc-miR-874 -0.154582   overtend       -
## 
## $XLOC_013030
##    chr          ID genes       miRNA        cor      pheno cor.mag
## 16   2 XLOC_013030 OVOL1 ssc-miR-874 -0.1202629        WBS       -
## 26   2 XLOC_013030 OVOL1 ssc-miR-874 -0.1202629 tenderness       -
## 33   2 XLOC_013030 OVOL1 ssc-miR-874 -0.1202629   overtend       -
## 
## $XLOC_009607
##    chr          ID  genes       miRNA        cor    pheno cor.mag
## 86  15 XLOC_009607 KLHL30 ssc-miR-874 -0.1280029 driploss       -
## 
## $XLOC_022575
##    chr          ID  genes       miRNA        cor     pheno cor.mag
## 42   7 XLOC_022575 ZNF451 ssc-miR-874 -0.1757576 dress_ptg       -
## 
## $XLOC_023083
##    chr          ID genes       miRNA        cor    pheno cor.mag
## 53   7 XLOC_023083   MAX ssc-miR-874 -0.0875502 num_ribs       -
## 
## $XLOC_009276
##    chr          ID genes       miRNA         cor   pheno cor.mag
## 69  15 XLOC_009276 UGGT1 ssc-miR-874 -0.09456006 protein       -
## 
## $XLOC_002373
##   chr          ID genes       miRNA       cor    pheno cor.mag
## 6   1 XLOC_002373  RXRA ssc-miR-874 -0.139686 car_bf10       -
## 
## $XLOC_008850
##    chr          ID  genes       miRNA        cor   pheno cor.mag
## 63  15 XLOC_008850 NIF3L1 ssc-miR-874 -0.1726908 protein       -
## 
## $XLOC_009284
##    chr          ID   genes       miRNA        cor   pheno cor.mag
## 70  15 XLOC_009284 PRPF40A ssc-miR-874 -0.1430449 protein       -
## 
## $XLOC_009537
##    chr          ID genes       miRNA        cor      pheno cor.mag
## 74  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376    protein       -
## 79  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376     ph_24h       -
## 84  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376   driploss       -
## 90  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376 cook_yield       -
## 93  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376   overtend       -
## 94  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376  juiciness       -
## 95  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376        WBS       -
## 96  15 XLOC_009537 PTPRN ssc-miR-874 -0.1230376 tenderness       -
## 
## $XLOC_011767
##    chr          ID  genes       miRNA        cor      pheno cor.mag
## 9    2 XLOC_011767 PPP1CA ssc-miR-874 -0.1373494        WBS       -
## 19   2 XLOC_011767 PPP1CA ssc-miR-874 -0.1373494 tenderness       -
## 
## $XLOC_013026
##    chr          ID   genes       miRNA        cor      pheno cor.mag
## 14   2 XLOC_013026 CCDC85B ssc-miR-874 -0.1220153        WBS       -
## 24   2 XLOC_013026 CCDC85B ssc-miR-874 -0.1220153 tenderness       -
## 31   2 XLOC_013026 CCDC85B ssc-miR-874 -0.1220153   overtend       -
## 
## $XLOC_013029
##    chr          ID genes       miRNA        cor      pheno cor.mag
## 15   2 XLOC_013029 SNX32 ssc-miR-874 -0.1296093        WBS       -
## 25   2 XLOC_013029 SNX32 ssc-miR-874 -0.1296093 tenderness       -
## 32   2 XLOC_013029 SNX32 ssc-miR-874 -0.1296093   overtend       -
```

## Save data


```r
save(rst.cor, rst.corfil, mrst.cor, mrst.corfil, splt.negcoloc, file="../16_mirna_mrna_pheno_corr.R")
```

