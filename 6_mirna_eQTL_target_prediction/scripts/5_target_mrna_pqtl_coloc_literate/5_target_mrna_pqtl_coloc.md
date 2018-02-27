**Script:** `5_target_mrna-pqtl-coloc.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`

**Date:**  12/13/17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/`, `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k`, `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z`

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

1. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`

**Input File(s):** 

1. `voom.Rdata`, `pQTL_60k.Rdata`, `inrange_function.Rdata`

1. `10_mrna_mirna_corr_rst.Rdata`

1. `MSUPRP_gpData_Ss11.Rdata`

**Output File Directory:** 

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

**Output File(s):** ``

1. `12_target_mrna_coloc_pqtl.Rdata`, `13_target_mrna_coloc_pqtl.txt`, `14_target_mrna_neg_coloc_pqtl.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to co-localize the significant target mRNAs with pQTL identified in the MSUPRP dataset

## Install libraries



```r
rm(list=ls())
library(methods)
library(limma)
library(edgeR)
library(gwaR)
library(parallel)
library(qvalue)
```

Session Information


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: CentOS release 6.8 (Final)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] qvalue_2.2.2 gwaR_1.0     edgeR_3.12.1 limma_3.26.9 knitr_1.17  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7      plyr_1.8.4       grid_3.2.0       gtable_0.2.0    
##  [5] magrittr_1.5     evaluate_0.10.1  scales_0.4.0     ggplot2_2.1.0   
##  [9] stringi_1.1.1    reshape2_1.4.2   splines_3.2.0    tools_3.2.0     
## [13] stringr_1.2.0    munsell_0.4.3    colorspace_1.3-2
```

## Load data

Load required R objects 


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/10_mrna_mirna_corr_rst.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")
ls()
```

```
##  [1] "dge"                 "GBLUP"               "GWAS"               
##  [4] "inrange"             "MSUPRP"              "MSUPRP168"          
##  [7] "QTL"                 "QTLpeaks"            "rst.corR"           
## [10] "sig.mrnaR"           "sig.mrnaR.names"     "sig.neg.mrnaR"      
## [13] "sig.neg.mrnaR.names" "summary.sigR"        "v"                  
## [16] "wcen"
```

## Analysis

### Annotation of corrleated target genes


```r
genes <- dge$genes
colnames(genes)[1] <- "chr"
annot <- do.call(rbind, lapply(names(sig.mrnaR), function(x) 
	data.frame(genes[sig.mrnaR[[x]],], miRNA=rep(x, length(sig.mrnaR[[x]])), 
		rst.corR[[x]][sig.mrnaR[[x]],c("cor", "pvalue", "qvalue")])))
annot$ID<-as.character(rownames(annot))
head(annot)
```

```
##             chr     start       end  width strand          ID    genes
## XLOC_001644   1 120807328 120946698 139371      - XLOC_001644    AP4E1
## XLOC_002556  10  24462354  24638203 175850      + XLOC_002556 PPP1R12B
## XLOC_004063  12  36434092  36496637  62546      + XLOC_004063    INTS2
## XLOC_014709   3  72240089  72281813  41725      + XLOC_014709     TIA1
## XLOC_007579  14 128962097 128966113   4017      + XLOC_007579   NANOS1
## XLOC_002971  10  46862232  46900082  37851      - XLOC_002971  SUV39H2
##                     locus           miRNA        cor       pvalue
## XLOC_001644 RLOC_00000903   ssc-let-7d-5p  0.1906535 2.662193e-04
## XLOC_002556 RLOC_00002598   ssc-let-7d-5p -0.2249726 1.689111e-05
## XLOC_004063 RLOC_00004295   ssc-let-7d-5p  0.1995619 1.353598e-04
## XLOC_014709 RLOC_00014880   ssc-let-7d-5p  0.1954728 1.852719e-04
## XLOC_007579 RLOC_00008098  ssc-miR-345-3p  0.2844104 5.352893e-08
## XLOC_002971 RLOC_00002822 ssc-miR-6782-3p -0.1575027 2.594113e-03
##                   qvalue
## XLOC_001644 3.894488e-02
## XLOC_002556 9.883916e-03
## XLOC_004063 3.613759e-02
## XLOC_014709 3.613759e-02
## XLOC_007579 5.830717e-05
## XLOC_002971 4.546732e-02
```

### Phenotypic QTL MSUPRP


```r
length(QTL)
```

```
## [1] 25
```

```r
names(QTL)
```

```
##  [1] "bf10_10wk"  "lrf_10wk"   "bf10_13wk"  "lrf_13wk"   "bf10_16wk" 
##  [6] "lma_16wk"   "lrf_16wk"   "bf10_19wk"  "lrf_19wk"   "bf10_22wk" 
## [11] "lrf_22wk"   "dress_ptg"  "cook_yield" "WBS"        "juiciness" 
## [16] "tenderness" "overtend"   "driploss"   "ph_24h"     "car_length"
## [21] "num_ribs"   "last_lum"   "car_bf10"   "loin"       "protein"
```

```r
QTL[[1]]
```

```
##             chr       pos         pval       qval
## ASGA0029650   6 144550432 3.070272e-06 0.04372783
## ASGA0029651   6 144640222 6.613699e-07 0.01562265
## ALGA0104402   6 147743951 7.312776e-07 0.01562265
```

pQTL peak information


```r
dim(QTLpeaks)
```

```
## [1] 46  9
```

```r
QTLpeaks$pheno.chr<-as.factor(make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE))
rownames(QTLpeaks)<-make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE)
head(QTLpeaks)
```

```
##                 pheno chr     start       end sig.snp         snp
## bf10_10wk.6 bf10_10wk   6 144550432 147743951       3 ASGA0029651
## lrf_10wk.6   lrf_10wk   6 114921710 153129074      51 ALGA0104402
## lrf_10wk.12  lrf_10wk  12  41734343  41734343       1 ASGA0054658
## bf10_13wk.3 bf10_13wk   3 112321459 112321459       1 H3GA0010564
## bf10_13wk.6 bf10_13wk   6 142218300 150042424      10 ALGA0104402
## lrf_13wk.6   lrf_13wk   6 108751759 151354003      17 ALGA0104402
##               pos.snp         pval         qval   pheno.chr
## bf10_10wk.6 144640222 6.613699e-07 1.562265e-02 bf10_10wk.6
## lrf_10wk.6  147743951 1.228888e-09 5.250669e-05  lrf_10wk.6
## lrf_10wk.12  41734343 4.619538e-05 4.028143e-02 lrf_10wk.12
## bf10_13wk.3 112321459 1.339311e-06 1.143386e-02 bf10_13wk.3
## bf10_13wk.6 147743951 8.468488e-08 3.049103e-03 bf10_13wk.6
## lrf_13wk.6  147743951 5.452684e-10 2.329768e-05  lrf_13wk.6
```

```r
table(QTLpeaks$pheno)
```

```
## 
##  bf10_10wk   lrf_10wk  bf10_13wk   lrf_13wk  bf10_16wk   lma_16wk 
##          1          2          2          1          1          1 
##   lrf_16wk  bf10_19wk   lrf_19wk  bf10_22wk   lrf_22wk  dress_ptg 
##          2          1          1          1          1          3 
## cook_yield        WBS  juiciness tenderness   overtend   driploss 
##          3          2          1          4          4          1 
##     ph_24h car_length   num_ribs   last_lum   car_bf10       loin 
##          1          2          1          5          2          2 
##    protein 
##          1
```

```r
table(QTLpeaks$chr)
```

```
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 15 
##  1  3  3  1  4 14  5  1  2  1  2  1  8
```

q-values pQTL GWAS


```r
qval <- GWAS$qvalue
dim(qval)
```

```
## [1] 42727    67
```

p-values pQTL GWAS


```r
pval <- GWAS$pvalue
dim(pval)
```

```
## [1] 42727    67
```

Standardized SNP effects pQTL GWA


```r
sdEff <- GWAS$Estimate
dim(sdEff)
```

```
## [1] 42727    67
```

```r
map<-MSUPRP$map
str(map)
```

```
## 'data.frame':	43130 obs. of  3 variables:
##  $ chr: Factor w/ 173 levels "0","1","10","11",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ pos: int  205163 261794 309120 289363 408640 381075 373626 239523 455096 477400 ...
##  $ snp: Factor w/ 8 levels "[A/C]","[A/G]",..: 2 1 2 1 7 7 8 1 2 7 ...
```

```r
map$chr<-gsub("X", "19", map$chr)
map$chr<-gsub("Y", "20", map$chr)
map$chr<-as.numeric(map$chr)
str(map)
```

```
## 'data.frame':	43130 obs. of  3 variables:
##  $ chr: num  1 1 1 1 1 1 1 1 1 1 ...
##  $ pos: int  205163 261794 309120 289363 408640 381075 373626 239523 455096 477400 ...
##  $ snp: Factor w/ 8 levels "[A/C]","[A/G]",..: 2 1 2 1 7 7 8 1 2 7 ...
```

Retain marker information for all pQTL


```r
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(map[names(sig[[x]]),], 
	std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]
length(sig)
```

```
## [1] 25
```

```r
names(sig)
```

```
##  [1] "bf10_10wk"  "lrf_10wk"   "bf10_13wk"  "lrf_13wk"   "bf10_16wk" 
##  [6] "lma_16wk"   "lrf_16wk"   "bf10_19wk"  "lrf_19wk"   "bf10_22wk" 
## [11] "lrf_22wk"   "dress_ptg"  "cook_yield" "WBS"        "juiciness" 
## [16] "tenderness" "overtend"   "driploss"   "ph_24h"     "car_length"
## [21] "num_ribs"   "last_lum"   "car_bf10"   "loin"       "protein"
```

If a pQTL contains more than one peak split each peak


```r
idx <- lapply(sig, function(x) as.numeric(names(table(x$chr))))
idx <- idx[unlist(lapply(idx, function(x) length(x) > 1))]

mp <- lapply(1:length(idx), function(x) lapply(1:length(idx[[x]]), function(y) sig[[names(idx)[x]]][sig[[names(idx)[x]]]$chr == idx[[x]][y],]))
names(mp) <- names(idx)

for (i in 1:length(mp)){
names(mp[[i]]) <- paste(names(mp)[[i]], 1:length(mp[[i]]), sep=".")
}

qtl <- sig[!names(sig) %in% names(mp)]
for (i in 1:length(mp)){
qtl <- c(qtl, mp[[i]])
}
```

pQTL genomic regions


```r
qtlP <- do.call(rbind, lapply(names(qtl), function(x) data.frame(pheno=strsplit(x, "[.]")[[1]][1], 
	chr=unique(qtl[[x]]$chr), start=min(qtl[[x]]$pos), end=max(qtl[[x]]$pos))))
rownames(qtlP) <- names(qtl)
qtlP <- qtlP[order(qtlP$chr),]
tmp <- list()
for(i in unique(qtlP$chr)){
	x <- qtlP[qtlP$chr == i,]
	tmp[[i]] <- x[order(x$start),]
}
qtlP <- do.call(rbind,tmp)
dim(qtlP)
```

```
## [1] 43  4
```

### Co-localization mRNA with pQTL 

Check all mRNA within a eQTL


```r
win <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single=NULL, range=c(start="start", end="end")))
names(win) <- rownames(qtlP)
```

mRNA overlaping left side of pQTL


```r
left <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="start", range=NULL))
names(left) <- rownames(qtlP)
```

mRNA overlaping right side of pQTL


```r
right <- lapply(1:nrow(qtlP), function(x) inrange(chr=qtlP[x,"chr"], start=qtlP[x,"start"], 
	end=qtlP[x,"end"], map=annot, single="end", range=NULL))
names(right) <- rownames(qtlP)
```

Merge all mRNA co-localizing with pQTL  


```r
# Merge within and left side
coloc <- lapply(names(win), function(x) rbind(win[[x]], 
	left[[x]][!as.character(left[[x]]$geneID) %in% as.character(win[[x]]$geneID) | 
	!as.character(left[[x]]$miRNA) == as.character(win[[x]]$miRNA),]))
names(coloc) <- names(win)

names(coloc[c(1,28)])
```

```
## [1] "car_bf10.1" "num_ribs"
```

For car_bf10.1 and num_ribs, there is one more right-overlapping peak than there is left-overlapping. 
So, when building the merged object, the two objects will have different lengths, triggering the warning 
and adding the last line of the longer list to the merged object.



```r
# Merge coloc and right side
coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]], 
	right[[x]][!as.character(rownames(right[[x]])) %in% as.character(rownames(coloc[[x]])) |
	!as.character(right[[x]]$miRNA) == as.character(coloc[[x]]$miRNA),]))
```

```
## Warning in as.character(right[[x]]$miRNA) == as.character(coloc[[x]]
## $miRNA): longer object length is not a multiple of shorter object length

## Warning in as.character(right[[x]]$miRNA) == as.character(coloc[[x]]
## $miRNA): longer object length is not a multiple of shorter object length
```

```r
names(coloc) <- names(win)

data.frame(left=sapply(1:length(left), function(x) nrow(left[[x]])), 
	right=sapply(1:length(right), function(x) nrow(right[[x]])),
	win=sapply(1:length(win), function(x) nrow(win[[x]])), 
	coloc=sapply(1:length(coloc), function(x) nrow(coloc[[x]])))
```

```
##    left right win coloc
## 1     6     7   6     7
## 2    10    10  10    10
## 3    10    10  10    10
## 4     7     7   7     7
## 5     0     0   0     0
## 6     0     0   0     0
## 7     0     0   0     0
## 8     0     0   0     0
## 9     0     0   0     0
## 10    0     0   0     0
## 11    0     0   0     0
## 12    0     0   0     0
## 13    1     1   1     1
## 14    1     1   1     1
## 15    1     1   1     1
## 16    1     1   1     1
## 17    0     0   0     0
## 18    0     0   0     0
## 19    0     0   0     0
## 20    0     0   0     0
## 21    0     0   0     0
## 22    0     0   0     0
## 23    0     0   0     0
## 24    0     0   0     0
## 25    0     0   0     0
## 26    0     0   0     0
## 27    7     7   7     7
## 28   11    12  11    12
## 29    0     0   0     0
## 30    0     0   0     0
## 31    0     0   0     0
## 32    0     0   0     0
## 33    2     2   2     2
## 34    0     0   0     0
## 35    0     0   0     0
## 36   17    17  17    17
## 37    4     4   4     4
## 38    7     7   7     7
## 39    4     4   4     4
## 40    2     2   2     2
## 41    1     1   1     1
## 42    1     1   1     1
## 43    1     1   1     1
```

Final list of mRNA targets significantly correlated with miRNA and co-localizing with a pQTL


```r
coloc <- do.call(rbind, lapply(names(coloc), function(x) data.frame(coloc[[x]], 
	pheno=rep(strsplit(x, "[.]")[[1]][1], nrow(coloc[[x]])))))
rownames(coloc) <- NULL
dim(coloc)
```

```
## [1] 96 13
```

```r
head(coloc)
```

```
##   chr     start       end  width strand          ID    genes         locus
## 1   1 273024154 273067252  43099      + XLOC_001198 ADAMTS13 RLOC_00002353
## 2   1 273703960 273797743  93784      + XLOC_001203     RXRA RLOC_00002368
## 3   1 273933365 274085480 152116      + XLOC_001204   COL5A1 RLOC_00002370
## 4   1 271446726 271453855   7130      - XLOC_002346     UCK1 RLOC_00002322
## 5   1 273303015 273476293 173279      - XLOC_002368     VAV2 RLOC_00002361
## 6   1 273796280 273796561    282      - XLOC_002373     RXRA RLOC_00002369
##         miRNA        cor      pvalue     qvalue    pheno
## 1 ssc-miR-874  0.1274188 0.014817681 0.01646063 car_bf10
## 2 ssc-miR-874  0.1256663 0.016248302 0.01714472 car_bf10
## 3 ssc-miR-874  0.1088719 0.037332748 0.02714232 car_bf10
## 4 ssc-miR-874 -0.1034684 0.047842001 0.03169112 car_bf10
## 5 ssc-miR-874  0.1428989 0.006278854 0.01100962 car_bf10
## 6 ssc-miR-874 -0.1396860 0.007553373 0.01216263 car_bf10
```

```r
table(coloc$miRNA)
```

```
## 
##   ssc-let-7d-5p  ssc-miR-345-3p ssc-miR-6782-3p     ssc-miR-874 
##               0               0               0              95 
##      ssc-miR-95 
##               1
```

```r
coloc[coloc$miRNA=="ssc-miR-874",c("genes", "miRNA", "qvalue", "pheno")]
```

```
##       genes       miRNA      qvalue      pheno
## 1  ADAMTS13 ssc-miR-874 0.016460629   car_bf10
## 2      RXRA ssc-miR-874 0.017144718   car_bf10
## 3    COL5A1 ssc-miR-874 0.027142317   car_bf10
## 4      UCK1 ssc-miR-874 0.031691123   car_bf10
## 5      VAV2 ssc-miR-874 0.011009619   car_bf10
## 6      RXRA ssc-miR-874 0.012162627   car_bf10
## 7     HMCN2 ssc-miR-874 0.016279821   car_bf10
## 8    CARNS1 ssc-miR-874 0.019611261        WBS
## 9    PPP1CA ssc-miR-874 0.012795924        WBS
## 10    BRMS1 ssc-miR-874 0.040503591        WBS
## 11    OVOL1 ssc-miR-874 0.009055118        WBS
## 12     RELA ssc-miR-874 0.004997449        WBS
## 13    ATG2A ssc-miR-874 0.008344526        WBS
## 14  CCDC85B ssc-miR-874 0.018768508        WBS
## 15    SNX32 ssc-miR-874 0.015775143        WBS
## 16    OVOL1 ssc-miR-874 0.019474731        WBS
## 17    SIPA1 ssc-miR-874 0.044040849        WBS
## 18   CARNS1 ssc-miR-874 0.019611261 tenderness
## 19   PPP1CA ssc-miR-874 0.012795924 tenderness
## 20    BRMS1 ssc-miR-874 0.040503591 tenderness
## 21    OVOL1 ssc-miR-874 0.009055118 tenderness
## 22     RELA ssc-miR-874 0.004997449 tenderness
## 23    ATG2A ssc-miR-874 0.008344526 tenderness
## 24  CCDC85B ssc-miR-874 0.018768508 tenderness
## 25    SNX32 ssc-miR-874 0.015775143 tenderness
## 26    OVOL1 ssc-miR-874 0.019474731 tenderness
## 27    SIPA1 ssc-miR-874 0.044040849 tenderness
## 28    OVOL1 ssc-miR-874 0.009055118   overtend
## 29     RELA ssc-miR-874 0.004997449   overtend
## 30    ATG2A ssc-miR-874 0.008344526   overtend
## 31  CCDC85B ssc-miR-874 0.018768508   overtend
## 32    SNX32 ssc-miR-874 0.015775143   overtend
## 33    OVOL1 ssc-miR-874 0.019474731   overtend
## 34    SIPA1 ssc-miR-874 0.044040849   overtend
## 35   GPR157 ssc-miR-874 0.008344526   lma_16wk
## 36    PQLC1 ssc-miR-874 0.009092896   lrf_13wk
## 37    PQLC1 ssc-miR-874 0.009092896   car_bf10
## 38    PQLC1 ssc-miR-874 0.009092896   lrf_10wk
## 39    PFDN6 ssc-miR-874 0.009441877  dress_ptg
## 40    LEMD2 ssc-miR-874 0.009874884  dress_ptg
## 41    NUDT3 ssc-miR-874 0.014626639  dress_ptg
## 42   ZNF451 ssc-miR-874 0.007715510  dress_ptg
## 43    LEMD2 ssc-miR-874 0.031659192  dress_ptg
## 44     ETV7 ssc-miR-874 0.017144718  dress_ptg
## 45   CDKN1A ssc-miR-874 0.048976188  dress_ptg
## 46   STOML1 ssc-miR-874 0.021253050   num_ribs
## 47     RGMA ssc-miR-874 0.016460629   num_ribs
## 48     JDP2 ssc-miR-874 0.035995446   num_ribs
## 49    CD276 ssc-miR-874 0.010435392   num_ribs
## 50    KHNYN ssc-miR-874 0.008344526   num_ribs
## 51     NOP9 ssc-miR-874 0.009358773   num_ribs
## 52     SPTB ssc-miR-874 0.034498485   num_ribs
## 53      MAX ssc-miR-874 0.048976188   num_ribs
## 54    PROX2 ssc-miR-874 0.025709409   num_ribs
## 55   ANGEL1 ssc-miR-874 0.024697145   num_ribs
## 56    POMT2 ssc-miR-874 0.019256573   num_ribs
## 57      PML ssc-miR-874 0.010970398   num_ribs
## 58   ATP8A2 ssc-miR-874 0.007715510  dress_ptg
## 59    USP12 ssc-miR-874 0.024170619  dress_ptg
## 60  PLEKHB2 ssc-miR-874 0.019897284    protein
## 61   HS6ST1 ssc-miR-874 0.009358773    protein
## 62    WDR33 ssc-miR-874 0.026845682    protein
## 63   NIF3L1 ssc-miR-874 0.007815367    protein
## 64    CFLAR ssc-miR-874 0.021356867    protein
## 65   IGFBP5 ssc-miR-874 0.021676537    protein
## 66    KCNE4 ssc-miR-874 0.022775615    protein
## 67  FAM168B ssc-miR-874 0.041786865    protein
## 68  PLEKHB2 ssc-miR-874 0.008599413    protein
## 69    UGGT1 ssc-miR-874 0.040503591    protein
## 70  PRPF40A ssc-miR-874 0.010970398    protein
## 71   METTL8 ssc-miR-874 0.031966317    protein
## 72    SF3B1 ssc-miR-874 0.009874884    protein
## 73  TMEM237 ssc-miR-874 0.008672540    protein
## 74    PTPRN ssc-miR-874 0.018249445    protein
## 75  FAM124B ssc-miR-874 0.016485794    protein
## 77   IGFBP5 ssc-miR-874 0.021676537     ph_24h
## 78    KCNE4 ssc-miR-874 0.022775615     ph_24h
## 79    PTPRN ssc-miR-874 0.018249445     ph_24h
## 80  FAM124B ssc-miR-874 0.016485794     ph_24h
## 81   IGFBP5 ssc-miR-874 0.021676537   driploss
## 82    KCNE4 ssc-miR-874 0.022775615   driploss
## 83  ATG16L1 ssc-miR-874 0.021253050   driploss
## 84    PTPRN ssc-miR-874 0.018249445   driploss
## 85  FAM124B ssc-miR-874 0.016485794   driploss
## 86   KLHL30 ssc-miR-874 0.016460629   driploss
## 87     HES6 ssc-miR-874 0.015117355   driploss
## 88   IGFBP5 ssc-miR-874 0.021676537 cook_yield
## 89    KCNE4 ssc-miR-874 0.022775615 cook_yield
## 90    PTPRN ssc-miR-874 0.018249445 cook_yield
## 91  FAM124B ssc-miR-874 0.016485794 cook_yield
## 92   IGFBP5 ssc-miR-874 0.021676537   overtend
## 93    PTPRN ssc-miR-874 0.018249445   overtend
## 94    PTPRN ssc-miR-874 0.018249445  juiciness
## 95    PTPRN ssc-miR-874 0.018249445        WBS
## 96    PTPRN ssc-miR-874 0.018249445 tenderness
```

Investigate the occurences of multiple rows:


```r
sort(table(as.character(coloc$genes)))
```

```
## 
## ADAMTS13   ANGEL1  ATG16L1   ATP8A2    CD276   CDKN1A    CFLAR   COL5A1 
##        1        1        1        1        1        1        1        1 
##     ETV7  FAM168B   GPR157     HES6    HMCN2   HS6ST1     JDP2    KHNYN 
##        1        1        1        1        1        1        1        1 
##   KLHL30      MAX   METTL8   NIF3L1     NOP9    NR4A2    NUDT3    PFDN6 
##        1        1        1        1        1        1        1        1 
##      PML    POMT2    PROX2  PRPF40A     RGMA    SF3B1     SPTB   STOML1 
##        1        1        1        1        1        1        1        1 
##  TMEM237     UCK1    UGGT1    USP12     VAV2    WDR33   ZNF451    BRMS1 
##        1        1        1        1        1        1        1        2 
##   CARNS1    LEMD2  PLEKHB2   PPP1CA     RXRA    ATG2A  CCDC85B    PQLC1 
##        2        2        2        2        2        3        3        3 
##     RELA    SIPA1    SNX32  FAM124B    KCNE4   IGFBP5    OVOL1    PTPRN 
##        3        3        3        4        4        5        6        8
```

```r
coloc$cor.mag<-ifelse(coloc$cor<0, "-", "+")
```

Convert miRNA names to human equivalent (based on seed sequence) for use in IPA:


```r
hsa.coloc<-coloc

miRid<-as.character(hsa.coloc$miRNA)
```

Convert ssc to hsa


```r
miRid<-gsub("ssc", "hsa", miRid)
```

ssc-miR-6782 is hsa-miR-6821


```r
miRid<-gsub("6782", "6821", miRid)
```

ssc-miR-874 is hsa-miR-874-3p


```r
miRid<-gsub("-874", "-874-3p", miRid)
```

ssc-miR-95 hsa-miR-95-3p


```r
miRid<-gsub("-95", "-95-3p", miRid)

hsa.coloc$hsa.miRNA<-miRid
```

Extract the pertinent columns for the significant negatively-associated mRNAs overlapping pQTLs:


```r
negcoloc<-coloc[coloc$cor < 0,c("chr", "ID", "genes", "miRNA", "cor", "pheno", "cor.mag")]
negcoloc
```

```
##    chr          ID   genes       miRNA         cor      pheno cor.mag
## 4    1 XLOC_002346    UCK1 ssc-miR-874 -0.10346842   car_bf10       -
## 6    1 XLOC_002373    RXRA ssc-miR-874 -0.13968602   car_bf10       -
## 8    2 XLOC_011766  CARNS1 ssc-miR-874 -0.11982475        WBS       -
## 9    2 XLOC_011767  PPP1CA ssc-miR-874 -0.13734940        WBS       -
## 11   2 XLOC_011794   OVOL1 ssc-miR-874 -0.15458196        WBS       -
## 14   2 XLOC_013026 CCDC85B ssc-miR-874 -0.12201533        WBS       -
## 15   2 XLOC_013029   SNX32 ssc-miR-874 -0.12960935        WBS       -
## 16   2 XLOC_013030   OVOL1 ssc-miR-874 -0.12026287        WBS       -
## 18   2 XLOC_011766  CARNS1 ssc-miR-874 -0.11982475 tenderness       -
## 19   2 XLOC_011767  PPP1CA ssc-miR-874 -0.13734940 tenderness       -
## 21   2 XLOC_011794   OVOL1 ssc-miR-874 -0.15458196 tenderness       -
## 24   2 XLOC_013026 CCDC85B ssc-miR-874 -0.12201533 tenderness       -
## 25   2 XLOC_013029   SNX32 ssc-miR-874 -0.12960935 tenderness       -
## 26   2 XLOC_013030   OVOL1 ssc-miR-874 -0.12026287 tenderness       -
## 28   2 XLOC_011794   OVOL1 ssc-miR-874 -0.15458196   overtend       -
## 31   2 XLOC_013026 CCDC85B ssc-miR-874 -0.12201533   overtend       -
## 32   2 XLOC_013029   SNX32 ssc-miR-874 -0.12960935   overtend       -
## 33   2 XLOC_013030   OVOL1 ssc-miR-874 -0.12026287   overtend       -
## 35   6 XLOC_019568  GPR157 ssc-miR-874 -0.16027747   lma_16wk       -
## 39   7 XLOC_021684   PFDN6 ssc-miR-874 -0.15122307  dress_ptg       -
## 40   7 XLOC_021691   LEMD2 ssc-miR-874 -0.14815626  dress_ptg       -
## 41   7 XLOC_021695   NUDT3 ssc-miR-874 -0.13209200  dress_ptg       -
## 42   7 XLOC_022575  ZNF451 ssc-miR-874 -0.17575758  dress_ptg       -
## 44   7 XLOC_022612    ETV7 ssc-miR-874 -0.12566630  dress_ptg       -
## 45   7 XLOC_022617  CDKN1A ssc-miR-874 -0.08755020  dress_ptg       -
## 50   7 XLOC_022870   KHNYN ssc-miR-874 -0.15998540   num_ribs       -
## 51   7 XLOC_022875    NOP9 ssc-miR-874 -0.15253742   num_ribs       -
## 53   7 XLOC_023083     MAX ssc-miR-874 -0.08755020   num_ribs       -
## 54   7 XLOC_023143   PROX2 ssc-miR-874 -0.11077035   num_ribs       -
## 56   7 XLOC_023161   POMT2 ssc-miR-874 -0.12099306   num_ribs       -
## 58  11 XLOC_003117  ATP8A2 ssc-miR-874 -0.17385907  dress_ptg       -
## 60  15 XLOC_008680 PLEKHB2 ssc-miR-874 -0.11924060    protein       -
## 62  15 XLOC_008685   WDR33 ssc-miR-874 -0.10945601    protein       -
## 63  15 XLOC_008850  NIF3L1 ssc-miR-874 -0.17269076    protein       -
## 65  15 XLOC_008918  IGFBP5 ssc-miR-874 -0.11617379    protein       -
## 68  15 XLOC_009274 PLEKHB2 ssc-miR-874 -0.15925520    protein       -
## 69  15 XLOC_009276   UGGT1 ssc-miR-874 -0.09456006    protein       -
## 70  15 XLOC_009284 PRPF40A ssc-miR-874 -0.14304491    protein       -
## 71  15 XLOC_009349  METTL8 ssc-miR-874 -0.10303030    protein       -
## 74  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760    protein       -
## 77  15 XLOC_008918  IGFBP5 ssc-miR-874 -0.11617379     ph_24h       -
## 79  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760     ph_24h       -
## 81  15 XLOC_008918  IGFBP5 ssc-miR-874 -0.11617379   driploss       -
## 84  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760   driploss       -
## 86  15 XLOC_009607  KLHL30 ssc-miR-874 -0.12800292   driploss       -
## 88  15 XLOC_008918  IGFBP5 ssc-miR-874 -0.11617379 cook_yield       -
## 90  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760 cook_yield       -
## 92  15 XLOC_008918  IGFBP5 ssc-miR-874 -0.11617379   overtend       -
## 93  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760   overtend       -
## 94  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760  juiciness       -
## 95  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760        WBS       -
## 96  15 XLOC_009537   PTPRN ssc-miR-874 -0.12303760 tenderness       -
```

```r
as.character(unique(negcoloc$genes))
```

```
##  [1] "UCK1"    "RXRA"    "CARNS1"  "PPP1CA"  "OVOL1"   "CCDC85B" "SNX32"  
##  [8] "GPR157"  "PFDN6"   "LEMD2"   "NUDT3"   "ZNF451"  "ETV7"    "CDKN1A" 
## [15] "KHNYN"   "NOP9"    "MAX"     "PROX2"   "POMT2"   "ATP8A2"  "PLEKHB2"
## [22] "WDR33"   "NIF3L1"  "IGFBP5"  "UGGT1"   "PRPF40A" "METTL8"  "PTPRN"  
## [29] "KLHL30"
```

```r
as.character(unique(negcoloc$pheno))
```

```
##  [1] "car_bf10"   "WBS"        "tenderness" "overtend"   "lma_16wk"  
##  [6] "dress_ptg"  "num_ribs"   "protein"    "ph_24h"     "driploss"  
## [11] "cook_yield" "juiciness"
```

```r
table(negcoloc$miRNA)
```

```
## 
##   ssc-let-7d-5p  ssc-miR-345-3p ssc-miR-6782-3p     ssc-miR-874 
##               0               0               0              52 
##      ssc-miR-95 
##               0
```

```r
table(as.character(negcoloc$pheno))
```

```
## 
##   car_bf10 cook_yield  dress_ptg   driploss  juiciness   lma_16wk 
##          2          2          7          3          1          1 
##   num_ribs   overtend     ph_24h    protein tenderness        WBS 
##          5          6          2          9          7          7
```

```r
table(as.character(negcoloc[negcoloc$miRNA=="ssc-miR-874" & negcoloc$chr=="2","pheno"]))
```

```
## 
##   overtend tenderness        WBS 
##          4          6          6
```

In my previous analysis, I had significant co-localization in 3 miR-6782-3p target genes in num_ribs traits

These three genes (STOML1, KHNYN, and RGMA) are present in this analysis, but their correlations did not reach 
statistical significance (q-val between 0.06 and 0.15).
## Save data


```r
save(annot, coloc, negcoloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/13_target_mrna_coloc_pqtl.Rdata")
write.table(hsa.coloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/14_target_mrna_coloc_pqtl.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(negcoloc, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/15_target_mrna_neg_coloc_pqtl.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

