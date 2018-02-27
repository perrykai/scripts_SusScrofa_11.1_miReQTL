**Script:** `3_filter_target_genes_exp.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`

**Date:**  12/13/17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

2. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom`

**Input File(s):** 

1. `4_filtered_targetscan_rst.Rdata`

2. `voom.Rdata`

**Output File Directory:** 

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

**Output File(s):** 

1. `5_filtered_targets_exp_rst.Rdata`

2. `6_DAVID_bkgd_gene_names_exp.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to filter the target prediction results based on genes expressed in the 168 LD samples DV obtained in her analysis.

The filtered list will dictate which genes will be included in the correlation analysis and pQTL co-localization.

## Install libraries


```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

rm(list=ls())

library(limma)
library(edgeR)
library(methods)
library(biomaRt)
```

## Load data

Load the list of filtered target prediction results:


```r
load("../4_filtered_targetscan_rst.Rdata")
```

Load DV's dge object to obtain her annotation file:


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
```

## Analysis



```r
genes<-dge$genes
genes$geneXID<-as.character(rownames(genes))
dim(genes)
```

```
## [1] 16152     9
```

```r
head(genes)
```

```
##             chr  start    end width strand            ID   genes
## XLOC_000001   1      7   5316  5310      + RLOC_00000001     TBP
## XLOC_000002   1  23368  40113 16746      + RLOC_00000003   PSMB1
## XLOC_000003   1 174035 175741  1707      + RLOC_00000006 FAM120B
## XLOC_000004   1 198992 211342 12351      + RLOC_00000007    DLL1
## XLOC_000005   1 227225 229005  1781      + RLOC_00000009    <NA>
## XLOC_000006   1 229894 235654  5761      + RLOC_00000010    <NA>
##                     locus     geneXID
## XLOC_000001 RLOC_00000001 XLOC_000001
## XLOC_000002 RLOC_00000003 XLOC_000002
## XLOC_000003 RLOC_00000006 XLOC_000003
## XLOC_000004 RLOC_00000007 XLOC_000004
## XLOC_000005 RLOC_00000009 XLOC_000005
## XLOC_000006 RLOC_00000010 XLOC_000006
```

How many genes have no name?


```r
sum(is.na(genes$genes))
```

```
## [1] 1552
```

How many genes have ensembl IDs, but no gene symbol?


```r
length(grep("ENSSSCG", as.character(genes$genes)))
```

```
## [1] 2274
```

Filter out the NA genes and the ensembl IDs:


```r
genes<-genes[!is.na(genes$genes),]
dim(genes)
```

```
## [1] 14600     9
```

```r
genes<-genes[!grepl("ENSSSCG", as.character(genes$genes)),]
```

Dimensions of retained genes object:


```r
dim(genes)
```

```
## [1] 12326     9
```

```r
head(genes)
```

```
##             chr  start    end  width strand            ID   genes
## XLOC_000001   1      7   5316   5310      + RLOC_00000001     TBP
## XLOC_000002   1  23368  40113  16746      + RLOC_00000003   PSMB1
## XLOC_000003   1 174035 175741   1707      + RLOC_00000006 FAM120B
## XLOC_000004   1 198992 211342  12351      + RLOC_00000007    DLL1
## XLOC_000013   1 555355 719928 164574      + RLOC_00000018   PHF10
## XLOC_000014   1 854559 883694  29136      + RLOC_00000020   THBS2
##                     locus     geneXID
## XLOC_000001 RLOC_00000001 XLOC_000001
## XLOC_000002 RLOC_00000003 XLOC_000002
## XLOC_000003 RLOC_00000006 XLOC_000003
## XLOC_000004 RLOC_00000007 XLOC_000004
## XLOC_000013 RLOC_00000018 XLOC_000013
## XLOC_000014 RLOC_00000020 XLOC_000014
```

How many unique genes:


```r
length(unique(genes$genes))
```

```
## [1] 11462
```

```r
tst<-targets.unique$"miR-128-3p"
head(tst)
```

```
##           a_Gene_ID miRNA_family_ID Site_type ensid_version
## 25  ENST00000002165      miR-128-3p   8mer-1a             6
## 62  ENST00000003302      miR-128-3p   7mer-m8             4
## 75  ENST00000004921      miR-128-3p   8mer-1a             3
## 164 ENST00000007722      miR-128-3p   7mer-m8             7
## 211 ENST00000012049      miR-128-3p   8mer-1a             5
## 230 ENST00000012443      miR-128-3p   7mer-m8             4
##     ensembl_gene_id ensembl_transc_id external_gene_name hgnc_symbol
## 25  ENSG00000001036   ENST00000002165              FUCA2       FUCA2
## 62  ENSG00000048028   ENST00000003302              USP28       USP28
## 75             <NA>              <NA>               <NA>        <NA>
## 164 ENSG00000005884   ENST00000007722              ITGA3       ITGA3
## 211 ENSG00000011478   ENST00000012049              QPCTL       QPCTL
## 230 ENSG00000011485   ENST00000012443              PPP5C       PPP5C
##     status
## 25   KNOWN
## 62   KNOWN
## 75    <NA>
## 164  KNOWN
## 211  KNOWN
## 230  KNOWN
```

```r
sum(tst$external_gene_name %in% genes$genes)
```

```
## [1] 1888
```

```r
dim(tst[tst$external_gene_name %in% genes$genes,])
```

```
## [1] 1888    9
```

Create a list of miRNA names to loop through:


```r
miRnm<-as.character(names(targets.unique))
miRnm
```

```
##  [1] "miR-200-3p/429"   "let-7-5p/98-5p"   "miR-128-3p"      
##  [4] "miR-140-5p"       "miR-6821-3p"      "miR-6888-3p"     
##  [7] "miR-874-3p"       "miR-345-3p"       "miR-6072/6891-3p"
## [10] "miR-1306-3p"      "miR-184"          "miR-190-5p"      
## [13] "miR-1468-5p"      "miR-95-3p"
```

Filter the genes targets based on those expressed in the dataset:


```r
targets.exp<-list()
targets.exp.sum<-list()
for(i in miRnm){
targets.exp[[i]]<-targets.unique[[i]][targets.unique[[i]]$external_gene_name %in% genes$genes,]
# Create summary file:
targets.exp.sum[[i]]$summary<-data.frame(
	# How many targets were input into biomart? 
	gene.input=nrow(targets.unique[[i]]),
	# How many targets were expressed?
	gene.output=nrow(targets.exp[[i]]),
	# What fraction of the targets input were expressed/retained?
	gene.prop=nrow(targets.exp[[i]])/nrow(targets.unique[[i]]))
}

str(targets.exp)
```

```
## List of 14
##  $ miR-200-3p/429  :'data.frame':	1541 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1541] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ miRNA_family_ID   : chr [1:1541] "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" ...
##   ..$ Site_type         : chr [1:1541] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1541] "3" "5" "8" "4" ...
##   ..$ ensembl_gene_id   : chr [1:1541] "ENSG00000003056" "ENSG00000002587" "ENSG00000001630" "ENSG00000075790" ...
##   ..$ ensembl_transc_id : chr [1:1541] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ external_gene_name: chr [1:1541] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ hgnc_symbol       : chr [1:1541] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ status            : chr [1:1541] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ let-7-5p/98-5p  :'data.frame':	891 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:891] "ENST00000001008" "ENST00000005386" "ENST00000007699" "ENST00000007722" ...
##   ..$ miRNA_family_ID   : chr [1:891] "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" ...
##   ..$ Site_type         : chr [1:891] "7mer-m8" "8mer-1a" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:891] "4" "3" "5" "7" ...
##   ..$ ensembl_gene_id   : chr [1:891] "ENSG00000004478" "ENSG00000005175" "ENSG00000006047" "ENSG00000005884" ...
##   ..$ ensembl_transc_id : chr [1:891] "ENST00000001008" "ENST00000005386" "ENST00000007699" "ENST00000007722" ...
##   ..$ external_gene_name: chr [1:891] "FKBP4" "RPAP3" "YBX2" "ITGA3" ...
##   ..$ hgnc_symbol       : chr [1:891] "FKBP4" "RPAP3" "YBX2" "ITGA3" ...
##   ..$ status            : chr [1:891] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-128-3p      :'data.frame':	1888 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1888] "ENST00000002165" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ miRNA_family_ID   : chr [1:1888] "miR-128-3p" "miR-128-3p" "miR-128-3p" "miR-128-3p" ...
##   ..$ Site_type         : chr [1:1888] "8mer-1a" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1888] "6" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1888] "ENSG00000001036" "ENSG00000048028" "ENSG00000005884" "ENSG00000011478" ...
##   ..$ ensembl_transc_id : chr [1:1888] "ENST00000002165" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ external_gene_name: chr [1:1888] "FUCA2" "USP28" "ITGA3" "QPCTL" ...
##   ..$ hgnc_symbol       : chr [1:1888] "FUCA2" "USP28" "ITGA3" "QPCTL" ...
##   ..$ status            : chr [1:1888] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-140-5p      :'data.frame':	957 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:957] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000053468" ...
##   ..$ miRNA_family_ID   : chr [1:957] "miR-140-5p" "miR-140-5p" "miR-140-5p" "miR-140-5p" ...
##   ..$ Site_type         : chr [1:957] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:957] "5" "3" "3" "3" ...
##   ..$ ensembl_gene_id   : chr [1:957] "ENSG00000002587" "ENSG00000001461" "ENSG00000005175" "ENSG00000048544" ...
##   ..$ ensembl_transc_id : chr [1:957] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000053468" ...
##   ..$ external_gene_name: chr [1:957] "HS3ST1" "NIPAL3" "RPAP3" "MRPS10" ...
##   ..$ hgnc_symbol       : chr [1:957] "HS3ST1" "NIPAL3" "RPAP3" "MRPS10" ...
##   ..$ status            : chr [1:957] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6821-3p     :'data.frame':	674 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:674] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000013222" ...
##   ..$ miRNA_family_ID   : chr [1:674] "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" ...
##   ..$ Site_type         : chr [1:674] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:674] "5" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:674] "ENSG00000002587" "ENSG00000075790" "ENSG00000010270" "ENSG00000241644" ...
##   ..$ ensembl_transc_id : chr [1:674] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000013222" ...
##   ..$ external_gene_name: chr [1:674] "HS3ST1" "BCAP29" "STARD3NL" "INMT" ...
##   ..$ hgnc_symbol       : chr [1:674] "HS3ST1" "BCAP29" "STARD3NL" "INMT" ...
##   ..$ status            : chr [1:674] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6888-3p     :'data.frame':	1711 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1711] "ENST00000002596" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ miRNA_family_ID   : chr [1:1711] "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" ...
##   ..$ Site_type         : chr [1:1711] "7mer-m8" "8mer-1a" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1711] "5" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1711] "ENSG00000002587" "ENSG00000048028" "ENSG00000005884" "ENSG00000011478" ...
##   ..$ ensembl_transc_id : chr [1:1711] "ENST00000002596" "ENST00000003302" "ENST00000007722" "ENST00000012049" ...
##   ..$ external_gene_name: chr [1:1711] "HS3ST1" "USP28" "ITGA3" "QPCTL" ...
##   ..$ hgnc_symbol       : chr [1:1711] "HS3ST1" "USP28" "ITGA3" "QPCTL" ...
##   ..$ status            : chr [1:1711] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-874-3p      :'data.frame':	1358 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1358] "ENST00000002596" "ENST00000046640" "ENST00000062104" "ENST00000078429" ...
##   ..$ miRNA_family_ID   : chr [1:1358] "miR-874-3p" "miR-874-3p" "miR-874-3p" "miR-874-3p" ...
##   ..$ Site_type         : chr [1:1358] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1358] "5" "3" "2" "4" ...
##   ..$ ensembl_gene_id   : chr [1:1358] "ENSG00000002587" "ENSG00000040531" "ENSG00000053438" "ENSG00000088256" ...
##   ..$ ensembl_transc_id : chr [1:1358] "ENST00000002596" "ENST00000046640" "ENST00000062104" "ENST00000078429" ...
##   ..$ external_gene_name: chr [1:1358] "HS3ST1" "CTNS" "NNAT" "GNA11" ...
##   ..$ hgnc_symbol       : chr [1:1358] "HS3ST1" "CTNS" "NNAT" "GNA11" ...
##   ..$ status            : chr [1:1358] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-345-3p      :'data.frame':	1237 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1237] "ENST00000005259" "ENST00000011653" "ENST00000033079" "ENST00000072869" ...
##   ..$ miRNA_family_ID   : chr [1:1237] "miR-345-3p" "miR-345-3p" "miR-345-3p" "miR-345-3p" ...
##   ..$ Site_type         : chr [1:1237] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1237] "4" "4" "3" "4" ...
##   ..$ ensembl_gene_id   : chr [1:1237] "ENSG00000075790" "ENSG00000010610" "ENSG00000031003" "ENSG00000133597" ...
##   ..$ ensembl_transc_id : chr [1:1237] "ENST00000005259" "ENST00000011653" "ENST00000033079" "ENST00000072869" ...
##   ..$ external_gene_name: chr [1:1237] "BCAP29" "CD4" "FAM13B" "ADCK2" ...
##   ..$ hgnc_symbol       : chr [1:1237] "BCAP29" "CD4" "FAM13B" "ADCK2" ...
##   ..$ status            : chr [1:1237] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6072/6891-3p:'data.frame':	1064 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1064] "ENST00000003912" "ENST00000007722" "ENST00000013125" "ENST00000014914" ...
##   ..$ miRNA_family_ID   : chr [1:1064] "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" ...
##   ..$ Site_type         : chr [1:1064] "7mer-m8" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1064] "3" "7" "4" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1064] "ENSG00000001461" "ENSG00000005884" "ENSG00000012983" "ENSG00000013588" ...
##   ..$ ensembl_transc_id : chr [1:1064] "ENST00000003912" "ENST00000007722" "ENST00000013125" "ENST00000014914" ...
##   ..$ external_gene_name: chr [1:1064] "NIPAL3" "ITGA3" "MAP4K5" "GPRC5A" ...
##   ..$ hgnc_symbol       : chr [1:1064] "NIPAL3" "ITGA3" "MAP4K5" "GPRC5A" ...
##   ..$ status            : chr [1:1064] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1306-3p     :'data.frame':	175 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:175] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000216487" ...
##   ..$ miRNA_family_ID   : chr [1:175] "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" ...
##   ..$ Site_type         : chr [1:175] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:175] "8" "4" "3" "7" ...
##   ..$ ensembl_gene_id   : chr [1:175] "ENSG00000006453" "ENSG00000088256" "ENSG00000100344" "ENSG00000100599" ...
##   ..$ ensembl_transc_id : chr [1:175] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000216487" ...
##   ..$ external_gene_name: chr [1:175] "BAIAP2L1" "GNA11" "PNPLA3" "RIN3" ...
##   ..$ hgnc_symbol       : chr [1:175] "BAIAP2L1" "GNA11" "PNPLA3" "RIN3" ...
##   ..$ status            : chr [1:175] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-184         :'data.frame':	219 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:219] "ENST00000200181" "ENST00000214869" "ENST00000216075" "ENST00000220592" ...
##   ..$ miRNA_family_ID   : chr [1:219] "miR-184" "miR-184" "miR-184" "miR-184" ...
##   ..$ Site_type         : chr [1:219] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:219] "3" "2" "6" "5" ...
##   ..$ ensembl_gene_id   : chr [1:219] "ENSG00000132470" "ENSG00000099203" "ENSG00000100253" "ENSG00000123908" ...
##   ..$ ensembl_transc_id : chr [1:219] "ENST00000200181" "ENST00000214869" "ENST00000216075" "ENST00000220592" ...
##   ..$ external_gene_name: chr [1:219] "ITGB4" "TMED1" "MIOX" "AGO2" ...
##   ..$ hgnc_symbol       : chr [1:219] "ITGB4" "TMED1" "MIOX" "AGO2" ...
##   ..$ status            : chr [1:219] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-190-5p      :'data.frame':	620 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:620] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000179259" ...
##   ..$ miRNA_family_ID   : chr [1:620] "miR-190-5p" "miR-190-5p" "miR-190-5p" "miR-190-5p" ...
##   ..$ Site_type         : chr [1:620] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:620] "5" "4" "3" "4" ...
##   ..$ ensembl_gene_id   : chr [1:620] "ENSG00000014919" "ENSG00000054965" "ENSG00000005001" "ENSG00000078237" ...
##   ..$ ensembl_transc_id : chr [1:620] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000179259" ...
##   ..$ external_gene_name: chr [1:620] "COX15" "FAM168A" "PRSS22" "TIGAR" ...
##   ..$ hgnc_symbol       : chr [1:620] "COX15" "FAM168A" "PRSS22" "TIGAR" ...
##   ..$ status            : chr [1:620] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1468-5p     :'data.frame':	248 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:248] "ENST00000173229" "ENST00000215838" "ENST00000219097" "ENST00000219252" ...
##   ..$ miRNA_family_ID   : chr [1:248] "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" ...
##   ..$ Site_type         : chr [1:248] "8mer-1a" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:248] "2" "3" "2" "5" ...
##   ..$ ensembl_gene_id   : chr [1:248] "ENSG00000065320" "ENSG00000185339" "ENSG00000091651" "ENSG00000102978" ...
##   ..$ ensembl_transc_id : chr [1:248] "ENST00000173229" "ENST00000215838" "ENST00000219097" "ENST00000219252" ...
##   ..$ external_gene_name: chr [1:248] "NTN1" "TCN2" "ORC6" "POLR2C" ...
##   ..$ hgnc_symbol       : chr [1:248] "NTN1" "TCN2" "ORC6" "POLR2C" ...
##   ..$ status            : chr [1:248] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-95-3p       :'data.frame':	126 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:126] "ENST00000216181" "ENST00000216797" "ENST00000236671" "ENST00000243964" ...
##   ..$ miRNA_family_ID   : chr [1:126] "miR-95-3p" "miR-95-3p" "miR-95-3p" "miR-95-3p" ...
##   ..$ Site_type         : chr [1:126] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:126] "5" "5" "2" "3" ...
##   ..$ ensembl_gene_id   : chr [1:126] "ENSG00000100345" "ENSG00000100906" "ENSG00000117984" "ENSG00000124140" ...
##   ..$ ensembl_transc_id : chr [1:126] "ENST00000216181" "ENST00000216797" "ENST00000236671" "ENST00000243964" ...
##   ..$ external_gene_name: chr [1:126] "MYH9" "NFKBIA" "CTSD" "SLC12A5" ...
##   ..$ hgnc_symbol       : chr [1:126] "MYH9" "NFKBIA" "CTSD" "SLC12A5" ...
##   ..$ status            : chr [1:126] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
```

```r
targets.exp.sum
```

```
## $`miR-200-3p/429`
## $`miR-200-3p/429`$summary
##   gene.input gene.output gene.prop
## 1       2301        1541 0.6697088
## 
## 
## $`let-7-5p/98-5p`
## $`let-7-5p/98-5p`$summary
##   gene.input gene.output gene.prop
## 1       1398         891 0.6373391
## 
## 
## $`miR-128-3p`
## $`miR-128-3p`$summary
##   gene.input gene.output gene.prop
## 1       2982        1888 0.6331321
## 
## 
## $`miR-140-5p`
## $`miR-140-5p`$summary
##   gene.input gene.output gene.prop
## 1       1518         957 0.6304348
## 
## 
## $`miR-6821-3p`
## $`miR-6821-3p`$summary
##   gene.input gene.output gene.prop
## 1       1118         674 0.6028623
## 
## 
## $`miR-6888-3p`
## $`miR-6888-3p`$summary
##   gene.input gene.output gene.prop
## 1       2816        1711 0.6075994
## 
## 
## $`miR-874-3p`
## $`miR-874-3p`$summary
##   gene.input gene.output gene.prop
## 1       2175        1358 0.6243678
## 
## 
## $`miR-345-3p`
## $`miR-345-3p`$summary
##   gene.input gene.output gene.prop
## 1       1950        1237  0.634359
## 
## 
## $`miR-6072/6891-3p`
## $`miR-6072/6891-3p`$summary
##   gene.input gene.output gene.prop
## 1       1680        1064 0.6333333
## 
## 
## $`miR-1306-3p`
## $`miR-1306-3p`$summary
##   gene.input gene.output gene.prop
## 1        264         175 0.6628788
## 
## 
## $`miR-184`
## $`miR-184`$summary
##   gene.input gene.output gene.prop
## 1        365         219       0.6
## 
## 
## $`miR-190-5p`
## $`miR-190-5p`$summary
##   gene.input gene.output gene.prop
## 1        979         620 0.6332993
## 
## 
## $`miR-1468-5p`
## $`miR-1468-5p`$summary
##   gene.input gene.output gene.prop
## 1        394         248 0.6294416
## 
## 
## $`miR-95-3p`
## $`miR-95-3p`$summary
##   gene.input gene.output gene.prop
## 1        196         126 0.6428571
```

---

Extract the gene symbols from the full expressed dataset to convert them to their ensembl gene IDs for use in DAVID:



```r
length(unique(genes$genes))
```

```
## [1] 11462
```

## Visualize
## Save data


```r
save(targets.exp, targets.exp.sum, file="../5_filtered_targets_exp_rst.Rdata")
write.table(unique(genes$genes), file="../6_DAVID_bkgd_gene_names_exp.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

