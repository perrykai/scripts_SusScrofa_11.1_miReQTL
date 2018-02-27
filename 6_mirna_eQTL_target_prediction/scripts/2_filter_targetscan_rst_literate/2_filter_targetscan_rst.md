**Script:** `2_filter_targetscan_rest.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`

**Date:**  `12/13/17`

**Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

**Input File(s):** `3_targetscan_mirna_eqtl_rst.txt`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`

**Output File(s):** `4_filtered_targetscan_rst.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to filter the miRNA target prediction results 
based on target site type, retaining 8mer-1a and 7mer-m8, the two most preferentially conserved and effective classes of target sites.  
The results will then be filtered to one miRNA-gene target pair to create a list of putative target genes for correlation analysis and co-localization analysis with pQTL.

Will need to translate transcript IDs back to gene names once filtering is complete. 
Also, final data will be organized as nested list by miRNA.

THIS ANALYSIS CONDUCTED IN R/3.1.0

## Install libraries


```r
library(biomaRt)
```

```
## Loading required package: methods
```

## Load data


```r
rm(list=ls())
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

tprst<-read.table("../3_targetscan_mirna_eqtl_rst.txt", 
	sep="\t",
	header=TRUE)
head(tprst)
```

```
##           a_Gene_ID miRNA_family_ID species_ID MSA_start MSA_end UTR_start
## 1            CDR1as     miR-1306-3p       9606      3654    3659      1424
## 2 ENST00000000412.3      miR-128-3p       9606      1207    1214       684
## 3 ENST00000000412.3     miR-1468-5p       9606      1732    1738       955
## 4 ENST00000000412.3  miR-200-3p/429       9606      1563    1569       881
## 5 ENST00000000412.3      miR-345-3p       9606      1824    1830      1040
## 6 ENST00000001008.4  let-7-5p/98-5p       9606      2007    2058      1059
##   UTR_end Group_num Site_type miRNA.in.this.species Group_type
## 1    1429         1      6mer                     x       6mer
## 2     689         2      6mer                     x       6mer
## 3     961         3   7mer-1a                     x    7mer-1a
## 4     887         4   7mer-m8                     x    7mer-m8
## 5    1046         5   7mer-1a                     x    7mer-1a
## 6    1065         6   7mer-m8                     x    7mer-m8
##   Species_in_this_group Species_in_this_group_with_this_site_type
## 1                  9606                                        NA
## 2                  9606                                        NA
## 3                  9606                                        NA
## 4                  9606                                        NA
## 5                  9606                                        NA
## 6                  9606                                        NA
##   ORF_overlap
## 1           0
## 2           0
## 3           1
## 4           0
## 5           1
## 6           0
```

```r
dim(tprst)
```

```
## [1] 146956     14
```

## Analysis

First, filter the target prediction results based on site type: retain 8mer and 7mer-m8:


```r
table(tprst$Site_type)
```

```
## 
##    6mer 7mer-1a 7mer-m8 8mer-1a 
##   73188   28974   33066   11728
```

```r
tp8m<-tprst[tprst$Site_type==c("8mer-1a", "7mer-m8"),c(1,2,9)]
tp8m$miRNA_family_ID<-as.character(tp8m$miRNA_family_ID)
tp8m$Site_type<-as.character(tp8m$Site_type)
```


The ensembl id's have a decimal at the end, indicating which version the ID comes from. In this case, the vesion number will be removed.

Remove the version numbers (".1", ".2" or ".3") after the ensembl ids to another column for each gene to make it compatible with biomaRt

See website http://useast.ensembl.org/Help/View?id=181 for details on the stability of ensembl IDs


```r
tp8m$ensid_version <- as.character(lapply(strsplit(as.character(tp8m$a_Gene_ID), '.', fixed = TRUE), "[", 2))
tp8m$a_Gene_ID<-as.character(gsub("\\..*","",tp8m$a_Gene_ID))

head(tp8m)
```

```
##          a_Gene_ID miRNA_family_ID Site_type ensid_version
## 4  ENST00000000412  miR-200-3p/429   7mer-m8             3
## 6  ENST00000001008  let-7-5p/98-5p   7mer-m8             4
## 25 ENST00000002165      miR-128-3p   8mer-1a             6
## 32 ENST00000002596      miR-140-5p   7mer-m8             5
## 35 ENST00000002596  miR-200-3p/429   8mer-1a             5
## 40 ENST00000002596     miR-6821-3p   7mer-m8             5
```

```r
table(tp8m$ensid_version)
```

```
## 
##    1   10   11   12    2    3    4    5    6    7    8    9 
## 5315   17   11    2 3987 4635 3637 2261 1259  632  288  116
```

How many putative targets are retained?


```r
dim(tp8m)
```

```
## [1] 22160     4
```

What fraction of results remain?


```r
nrow(tp8m)/nrow(tprst)
```

```
## [1] 0.1507934
```

Excellent, much easier to work with!

Do I still have putative targets for all my miRNA families?


```r
length(unique(tp8m$miRNA_family_ID))
```

```
## [1] 14
```

```r
miRnm<-unique(as.character(tp8m$miRNA_family_ID))
```

How many of each miRNA do I have targets for?


```r
for (i in miRnm){
	print(sum(tp8m$miRNA_family_ID == i))
}
```

```
## [1] 2633
## [1] 1479
## [1] 3374
## [1] 1605
## [1] 1188
## [1] 3248
## [1] 2450
## [1] 2111
## [1] 1802
## [1] 264
## [1] 377
## [1] 1030
## [1] 401
## [1] 198
```

If our concern in a gene list of targets for each miRNA, all that needs to be done is to create a list 
containing the transcript IDs for each putative target.


```r
targets.unique<-list()
for(i in miRnm){
	targets.unique[[i]]<-tp8m[tp8m$miRNA_family_ID==i,]
	targets.unique[[i]]<-subset(targets.unique[[i]], !duplicated(a_Gene_ID))
}
str(targets.unique)
```

```
## List of 14
##  $ miR-200-3p/429  :'data.frame':	2301 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:2301] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ miRNA_family_ID: chr [1:2301] "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" ...
##   ..$ Site_type      : chr [1:2301] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:2301] "3" "5" "8" "4" ...
##  $ let-7-5p/98-5p  :'data.frame':	1398 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:1398] "ENST00000001008" "ENST00000002829" "ENST00000005386" "ENST00000007699" ...
##   ..$ miRNA_family_ID: chr [1:1398] "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" ...
##   ..$ Site_type      : chr [1:1398] "7mer-m8" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:1398] "4" "3" "3" "5" ...
##  $ miR-128-3p      :'data.frame':	2982 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:2982] "ENST00000002165" "ENST00000003302" "ENST00000004921" "ENST00000007722" ...
##   ..$ miRNA_family_ID: chr [1:2982] "miR-128-3p" "miR-128-3p" "miR-128-3p" "miR-128-3p" ...
##   ..$ Site_type      : chr [1:2982] "8mer-1a" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:2982] "6" "4" "3" "7" ...
##  $ miR-140-5p      :'data.frame':	1518 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:1518] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000025301" ...
##   ..$ miRNA_family_ID: chr [1:1518] "miR-140-5p" "miR-140-5p" "miR-140-5p" "miR-140-5p" ...
##   ..$ Site_type      : chr [1:1518] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version  : chr [1:1518] "5" "3" "3" "2" ...
##  $ miR-6821-3p     :'data.frame':	1118 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:1118] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000013222" ...
##   ..$ miRNA_family_ID: chr [1:1118] "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" ...
##   ..$ Site_type      : chr [1:1118] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:1118] "5" "4" "7" "5" ...
##  $ miR-6888-3p     :'data.frame':	2816 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:2816] "ENST00000002596" "ENST00000003302" "ENST00000004980" "ENST00000007722" ...
##   ..$ miRNA_family_ID: chr [1:2816] "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" ...
##   ..$ Site_type      : chr [1:2816] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:2816] "5" "4" "5" "7" ...
##  $ miR-874-3p      :'data.frame':	2175 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:2175] "ENST00000002596" "ENST00000004921" "ENST00000008180" "ENST00000020926" ...
##   ..$ miRNA_family_ID: chr [1:2175] "miR-874-3p" "miR-874-3p" "miR-874-3p" "miR-874-3p" ...
##   ..$ Site_type      : chr [1:2175] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:2175] "5" "3" "9" "3" ...
##  $ miR-345-3p      :'data.frame':	1950 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:1950] "ENST00000002829" "ENST00000005259" "ENST00000011653" "ENST00000020926" ...
##   ..$ miRNA_family_ID: chr [1:1950] "miR-345-3p" "miR-345-3p" "miR-345-3p" "miR-345-3p" ...
##   ..$ Site_type      : chr [1:1950] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:1950] "3" "4" "4" "3" ...
##  $ miR-6072/6891-3p:'data.frame':	1680 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:1680] "ENST00000003912" "ENST00000007722" "ENST00000012443" "ENST00000013125" ...
##   ..$ miRNA_family_ID: chr [1:1680] "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" ...
##   ..$ Site_type      : chr [1:1680] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version  : chr [1:1680] "3" "7" "4" "4" ...
##  $ miR-1306-3p     :'data.frame':	264 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:264] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000216487" ...
##   ..$ miRNA_family_ID: chr [1:264] "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" ...
##   ..$ Site_type      : chr [1:264] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:264] "8" "4" "3" "7" ...
##  $ miR-184         :'data.frame':	365 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:365] "ENST00000013070" "ENST00000200181" "ENST00000214869" "ENST00000216075" ...
##   ..$ miRNA_family_ID: chr [1:365] "miR-184" "miR-184" "miR-184" "miR-184" ...
##   ..$ Site_type      : chr [1:365] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:365] "6" "3" "2" "6" ...
##  $ miR-190-5p      :'data.frame':	979 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:979] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000179259" ...
##   ..$ miRNA_family_ID: chr [1:979] "miR-190-5p" "miR-190-5p" "miR-190-5p" "miR-190-5p" ...
##   ..$ Site_type      : chr [1:979] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:979] "5" "4" "3" "4" ...
##  $ miR-1468-5p     :'data.frame':	394 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:394] "ENST00000173229" "ENST00000215838" "ENST00000217260" "ENST00000219097" ...
##   ..$ miRNA_family_ID: chr [1:394] "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" ...
##   ..$ Site_type      : chr [1:394] "8mer-1a" "7mer-m8" "8mer-1a" "8mer-1a" ...
##   ..$ ensid_version  : chr [1:394] "2" "3" "4" "2" ...
##  $ miR-95-3p       :'data.frame':	196 obs. of  4 variables:
##   ..$ a_Gene_ID      : chr [1:196] "ENST00000216181" "ENST00000216797" "ENST00000236671" "ENST00000242109" ...
##   ..$ miRNA_family_ID: chr [1:196] "miR-95-3p" "miR-95-3p" "miR-95-3p" "miR-95-3p" ...
##   ..$ Site_type      : chr [1:196] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version  : chr [1:196] "5" "5" "2" "3" ...
```

---

Convert the ensembl ID to human gene ID:

Use biomaRt package (v/2.20.0) to obtain annotation information for Ensembl dataset BLAST hits:


```r
mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="jul2016.archive.ensembl.org", dataset="hsapiens_gene_ensembl")

filters<-listFilters(mart)
dim(filters)
```

```
## [1] 310   2
```

```r
filters[60:70,]
```

```
##                     name                                     description
## 60          with_uniparc                              with UniParc ID(s)
## 61 with_uniprot_genename                       with UniProt Gene Name(s)
## 62 with_uniprotswissprot             with UniProt/SwissProt Accession(s)
## 63  with_uniprotsptrembl                with UniProt/TrEMBL Accession(s)
## 64         with_wikigene                             with WikiGene ID(s)
## 65       ensembl_gene_id       Ensembl Gene ID(s) [e.g. ENSG00000139618]
## 66 ensembl_transcript_id Ensembl Transcript ID(s) [e.g. ENST00000380152]
## 67    ensembl_peptide_id    Ensembl protein ID(s) [e.g. ENSP00000369497]
## 68       ensembl_exon_id       Ensembl exon ID(s) [e.g. ENSE00001508081]
## 69               hgnc_id                     HGNC ID(s) [e.g. HGNC:8030]
## 70           hgnc_symbol                      HGNC symbol(s) [e.g. NTN3]
```

```r
attributes<-listAttributes(mart)
dim(attributes)
```

```
## [1] 1464    2
```

```r
targets.rst <- list()
system.time(
for (i in miRnm){
targets.rst[[i]]<- getBM(attributes=c("ensembl_transcript_id", 
		"ensembl_gene_id", 
		"external_gene_name",
		"transcript_status", 
		"transcript_version", 
		"ens_hs_transcript",
		"hgnc_symbol",
		"hgnc_transcript_name",
		"gene_biotype"),
		filters="ensembl_transcript_id", 
		values=targets.unique[[i]]$a_Gene_ID, 
		mart=mart)
# Use match function to obtain the gene id, transcript id, gene name, hgnc gene symbol(check for equality), and gene status for each ensembl match and add it as a column to targets.unique:
targets.unique[[i]]$ensembl_gene_id<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "ensembl_gene_id"])
targets.unique[[i]]$ensembl_transc_id<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "ensembl_transcript_id"])
targets.unique[[i]]$external_gene_name<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "external_gene_name"])
targets.unique[[i]]$hgnc_symbol<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "hgnc_symbol"])
targets.unique[[i]]$status<-as.character(targets.rst[[i]][match(targets.unique[[i]]$a_Gene_ID, targets.rst[[i]]$ensembl_transcript_id), "transcript_status"])
})
```

```
##    user  system elapsed 
##   0.281   0.014  26.275
```

```r
str(targets.unique)
```

```
## List of 14
##  $ miR-200-3p/429  :'data.frame':	2301 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:2301] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ miRNA_family_ID   : chr [1:2301] "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" "miR-200-3p/429" ...
##   ..$ Site_type         : chr [1:2301] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:2301] "3" "5" "8" "4" ...
##   ..$ ensembl_gene_id   : chr [1:2301] "ENSG00000003056" "ENSG00000002587" "ENSG00000001630" "ENSG00000075790" ...
##   ..$ ensembl_transc_id : chr [1:2301] "ENST00000000412" "ENST00000002596" "ENST00000003100" "ENST00000005259" ...
##   ..$ external_gene_name: chr [1:2301] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ hgnc_symbol       : chr [1:2301] "M6PR" "HS3ST1" "CYP51A1" "BCAP29" ...
##   ..$ status            : chr [1:2301] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ let-7-5p/98-5p  :'data.frame':	1398 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1398] "ENST00000001008" "ENST00000002829" "ENST00000005386" "ENST00000007699" ...
##   ..$ miRNA_family_ID   : chr [1:1398] "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" "let-7-5p/98-5p" ...
##   ..$ Site_type         : chr [1:1398] "7mer-m8" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1398] "4" "3" "3" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1398] "ENSG00000004478" "ENSG00000001617" "ENSG00000005175" "ENSG00000006047" ...
##   ..$ ensembl_transc_id : chr [1:1398] "ENST00000001008" "ENST00000002829" "ENST00000005386" "ENST00000007699" ...
##   ..$ external_gene_name: chr [1:1398] "FKBP4" "SEMA3F" "RPAP3" "YBX2" ...
##   ..$ hgnc_symbol       : chr [1:1398] "FKBP4" "SEMA3F" "RPAP3" "YBX2" ...
##   ..$ status            : chr [1:1398] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-128-3p      :'data.frame':	2982 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:2982] "ENST00000002165" "ENST00000003302" "ENST00000004921" "ENST00000007722" ...
##   ..$ miRNA_family_ID   : chr [1:2982] "miR-128-3p" "miR-128-3p" "miR-128-3p" "miR-128-3p" ...
##   ..$ Site_type         : chr [1:2982] "8mer-1a" "7mer-m8" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:2982] "6" "4" "3" "7" ...
##   ..$ ensembl_gene_id   : chr [1:2982] "ENSG00000001036" "ENSG00000048028" NA "ENSG00000005884" ...
##   ..$ ensembl_transc_id : chr [1:2982] "ENST00000002165" "ENST00000003302" NA "ENST00000007722" ...
##   ..$ external_gene_name: chr [1:2982] "FUCA2" "USP28" NA "ITGA3" ...
##   ..$ hgnc_symbol       : chr [1:2982] "FUCA2" "USP28" NA "ITGA3" ...
##   ..$ status            : chr [1:2982] "KNOWN" "KNOWN" NA "KNOWN" ...
##  $ miR-140-5p      :'data.frame':	1518 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1518] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000025301" ...
##   ..$ miRNA_family_ID   : chr [1:1518] "miR-140-5p" "miR-140-5p" "miR-140-5p" "miR-140-5p" ...
##   ..$ Site_type         : chr [1:1518] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1518] "5" "3" "3" "2" ...
##   ..$ ensembl_gene_id   : chr [1:1518] "ENSG00000002587" "ENSG00000001461" "ENSG00000005175" "ENSG00000023516" ...
##   ..$ ensembl_transc_id : chr [1:1518] "ENST00000002596" "ENST00000003912" "ENST00000005386" "ENST00000025301" ...
##   ..$ external_gene_name: chr [1:1518] "HS3ST1" "NIPAL3" "RPAP3" "AKAP11" ...
##   ..$ hgnc_symbol       : chr [1:1518] "HS3ST1" "NIPAL3" "RPAP3" "AKAP11" ...
##   ..$ status            : chr [1:1518] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6821-3p     :'data.frame':	1118 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1118] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000013222" ...
##   ..$ miRNA_family_ID   : chr [1:1118] "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" "miR-6821-3p" ...
##   ..$ Site_type         : chr [1:1118] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1118] "5" "4" "7" "5" ...
##   ..$ ensembl_gene_id   : chr [1:1118] "ENSG00000002587" "ENSG00000075790" "ENSG00000010270" "ENSG00000241644" ...
##   ..$ ensembl_transc_id : chr [1:1118] "ENST00000002596" "ENST00000005259" "ENST00000009041" "ENST00000013222" ...
##   ..$ external_gene_name: chr [1:1118] "HS3ST1" "BCAP29" "STARD3NL" "INMT" ...
##   ..$ hgnc_symbol       : chr [1:1118] "HS3ST1" "BCAP29" "STARD3NL" "INMT" ...
##   ..$ status            : chr [1:1118] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6888-3p     :'data.frame':	2816 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:2816] "ENST00000002596" "ENST00000003302" "ENST00000004980" "ENST00000007722" ...
##   ..$ miRNA_family_ID   : chr [1:2816] "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" "miR-6888-3p" ...
##   ..$ Site_type         : chr [1:2816] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:2816] "5" "4" "5" "7" ...
##   ..$ ensembl_gene_id   : chr [1:2816] "ENSG00000002587" "ENSG00000048028" NA "ENSG00000005884" ...
##   ..$ ensembl_transc_id : chr [1:2816] "ENST00000002596" "ENST00000003302" NA "ENST00000007722" ...
##   ..$ external_gene_name: chr [1:2816] "HS3ST1" "USP28" NA "ITGA3" ...
##   ..$ hgnc_symbol       : chr [1:2816] "HS3ST1" "USP28" NA "ITGA3" ...
##   ..$ status            : chr [1:2816] "KNOWN" "KNOWN" NA "KNOWN" ...
##  $ miR-874-3p      :'data.frame':	2175 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:2175] "ENST00000002596" "ENST00000004921" "ENST00000008180" "ENST00000020926" ...
##   ..$ miRNA_family_ID   : chr [1:2175] "miR-874-3p" "miR-874-3p" "miR-874-3p" "miR-874-3p" ...
##   ..$ Site_type         : chr [1:2175] "7mer-m8" "8mer-1a" "8mer-1a" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:2175] "5" "3" "9" "3" ...
##   ..$ ensembl_gene_id   : chr [1:2175] "ENSG00000002587" NA "ENSG00000008517" "ENSG00000019505" ...
##   ..$ ensembl_transc_id : chr [1:2175] "ENST00000002596" NA "ENST00000008180" "ENST00000020926" ...
##   ..$ external_gene_name: chr [1:2175] "HS3ST1" NA "IL32" "SYT13" ...
##   ..$ hgnc_symbol       : chr [1:2175] "HS3ST1" NA "IL32" "SYT13" ...
##   ..$ status            : chr [1:2175] "KNOWN" NA "KNOWN" "KNOWN" ...
##  $ miR-345-3p      :'data.frame':	1950 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1950] "ENST00000002829" "ENST00000005259" "ENST00000011653" "ENST00000020926" ...
##   ..$ miRNA_family_ID   : chr [1:1950] "miR-345-3p" "miR-345-3p" "miR-345-3p" "miR-345-3p" ...
##   ..$ Site_type         : chr [1:1950] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:1950] "3" "4" "4" "3" ...
##   ..$ ensembl_gene_id   : chr [1:1950] "ENSG00000001617" "ENSG00000075790" "ENSG00000010610" "ENSG00000019505" ...
##   ..$ ensembl_transc_id : chr [1:1950] "ENST00000002829" "ENST00000005259" "ENST00000011653" "ENST00000020926" ...
##   ..$ external_gene_name: chr [1:1950] "SEMA3F" "BCAP29" "CD4" "SYT13" ...
##   ..$ hgnc_symbol       : chr [1:1950] "SEMA3F" "BCAP29" "CD4" "SYT13" ...
##   ..$ status            : chr [1:1950] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-6072/6891-3p:'data.frame':	1680 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:1680] "ENST00000003912" "ENST00000007722" "ENST00000012443" "ENST00000013125" ...
##   ..$ miRNA_family_ID   : chr [1:1680] "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" "miR-6072/6891-3p" ...
##   ..$ Site_type         : chr [1:1680] "7mer-m8" "7mer-m8" "7mer-m8" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:1680] "3" "7" "4" "4" ...
##   ..$ ensembl_gene_id   : chr [1:1680] "ENSG00000001461" "ENSG00000005884" "ENSG00000011485" "ENSG00000012983" ...
##   ..$ ensembl_transc_id : chr [1:1680] "ENST00000003912" "ENST00000007722" "ENST00000012443" "ENST00000013125" ...
##   ..$ external_gene_name: chr [1:1680] "NIPAL3" "ITGA3" "PPP5C" "MAP4K5" ...
##   ..$ hgnc_symbol       : chr [1:1680] "NIPAL3" "ITGA3" "PPP5C" "MAP4K5" ...
##   ..$ status            : chr [1:1680] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1306-3p     :'data.frame':	264 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:264] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000216487" ...
##   ..$ miRNA_family_ID   : chr [1:264] "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" "miR-1306-3p" ...
##   ..$ Site_type         : chr [1:264] "7mer-m8" "8mer-1a" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:264] "8" "4" "3" "7" ...
##   ..$ ensembl_gene_id   : chr [1:264] "ENSG00000006453" "ENSG00000088256" "ENSG00000100344" "ENSG00000100599" ...
##   ..$ ensembl_transc_id : chr [1:264] "ENST00000005260" "ENST00000078429" "ENST00000216180" "ENST00000216487" ...
##   ..$ external_gene_name: chr [1:264] "BAIAP2L1" "GNA11" "PNPLA3" "RIN3" ...
##   ..$ hgnc_symbol       : chr [1:264] "BAIAP2L1" "GNA11" "PNPLA3" "RIN3" ...
##   ..$ status            : chr [1:264] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-184         :'data.frame':	365 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:365] "ENST00000013070" "ENST00000200181" "ENST00000214869" "ENST00000216075" ...
##   ..$ miRNA_family_ID   : chr [1:365] "miR-184" "miR-184" "miR-184" "miR-184" ...
##   ..$ Site_type         : chr [1:365] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:365] "6" "3" "2" "6" ...
##   ..$ ensembl_gene_id   : chr [1:365] "ENSG00000012963" "ENSG00000132470" "ENSG00000099203" "ENSG00000100253" ...
##   ..$ ensembl_transc_id : chr [1:365] "ENST00000013070" "ENST00000200181" "ENST00000214869" "ENST00000216075" ...
##   ..$ external_gene_name: chr [1:365] "UBR7" "ITGB4" "TMED1" "MIOX" ...
##   ..$ hgnc_symbol       : chr [1:365] "UBR7" "ITGB4" "TMED1" "MIOX" ...
##   ..$ status            : chr [1:365] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-190-5p      :'data.frame':	979 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:979] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000179259" ...
##   ..$ miRNA_family_ID   : chr [1:979] "miR-190-5p" "miR-190-5p" "miR-190-5p" "miR-190-5p" ...
##   ..$ Site_type         : chr [1:979] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:979] "5" "4" "3" "4" ...
##   ..$ ensembl_gene_id   : chr [1:979] "ENSG00000014919" "ENSG00000054965" "ENSG00000005001" "ENSG00000078237" ...
##   ..$ ensembl_transc_id : chr [1:979] "ENST00000016171" "ENST00000064778" "ENST00000161006" "ENST00000179259" ...
##   ..$ external_gene_name: chr [1:979] "COX15" "FAM168A" "PRSS22" "TIGAR" ...
##   ..$ hgnc_symbol       : chr [1:979] "COX15" "FAM168A" "PRSS22" "TIGAR" ...
##   ..$ status            : chr [1:979] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-1468-5p     :'data.frame':	394 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:394] "ENST00000173229" "ENST00000215838" "ENST00000217260" "ENST00000219097" ...
##   ..$ miRNA_family_ID   : chr [1:394] "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" "miR-1468-5p" ...
##   ..$ Site_type         : chr [1:394] "8mer-1a" "7mer-m8" "8mer-1a" "8mer-1a" ...
##   ..$ ensid_version     : chr [1:394] "2" "3" "4" "2" ...
##   ..$ ensembl_gene_id   : chr [1:394] "ENSG00000065320" "ENSG00000185339" "ENSG00000101282" "ENSG00000091651" ...
##   ..$ ensembl_transc_id : chr [1:394] "ENST00000173229" "ENST00000215838" "ENST00000217260" "ENST00000219097" ...
##   ..$ external_gene_name: chr [1:394] "NTN1" "TCN2" "RSPO4" "ORC6" ...
##   ..$ hgnc_symbol       : chr [1:394] "NTN1" "TCN2" "RSPO4" "ORC6" ...
##   ..$ status            : chr [1:394] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
##  $ miR-95-3p       :'data.frame':	196 obs. of  9 variables:
##   ..$ a_Gene_ID         : chr [1:196] "ENST00000216181" "ENST00000216797" "ENST00000236671" "ENST00000242109" ...
##   ..$ miRNA_family_ID   : chr [1:196] "miR-95-3p" "miR-95-3p" "miR-95-3p" "miR-95-3p" ...
##   ..$ Site_type         : chr [1:196] "7mer-m8" "7mer-m8" "7mer-m8" "7mer-m8" ...
##   ..$ ensid_version     : chr [1:196] "5" "5" "2" "3" ...
##   ..$ ensembl_gene_id   : chr [1:196] "ENSG00000100345" "ENSG00000100906" "ENSG00000117984" "ENSG00000122548" ...
##   ..$ ensembl_transc_id : chr [1:196] "ENST00000216181" "ENST00000216797" "ENST00000236671" "ENST00000242109" ...
##   ..$ external_gene_name: chr [1:196] "MYH9" "NFKBIA" "CTSD" "KIAA0087" ...
##   ..$ hgnc_symbol       : chr [1:196] "MYH9" "NFKBIA" "CTSD" "KIAA0087" ...
##   ..$ status            : chr [1:196] "KNOWN" "KNOWN" "KNOWN" "KNOWN" ...
```

```r
targets.rst.sum<-list()
for (i in miRnm){
targets.rst.sum[[i]]$summary<-data.frame(
	# How many transcript IDs were input into biomart? 
	transc.id.input=nrow(targets.unique[[i]]),
	# How many transcript IDs are NA (possibly retired)?
	transc.id.na=sum(is.na(targets.unique[[i]]$ensembl_transc_id)),
	# What fraction of the input transcript IDs were NA?
	transc.id.na.prop=sum(is.na(targets.unique[[i]]$ensembl_transc_id)) / nrow(targets.unique[[i]]),
	# Of the transcripts remaining (not NA), check for any differences between the input transcript ID and the output ensembl_transc_id:
	transc.id.diff=sum(targets.unique[[i]]$a_Gene_ID != targets.unique[[i]]$ensembl_transc_id, na.rm=TRUE),
	# Do any of the gene names differ between the "external gene name" and the "HGNC symbol" columns?
	gene.name.diff=sum(targets.unique[[i]]$external_gene_name != targets.unique[[i]]$hgnc_symbol, na.rm=TRUE))
}

targets.rst.sum
```

```
## $`miR-200-3p/429`
## $`miR-200-3p/429`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            2301          143        0.06214689              0
##   gene.name.diff
## 1             11
## 
## 
## $`let-7-5p/98-5p`
## $`let-7-5p/98-5p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            1398          106         0.0758226              0
##   gene.name.diff
## 1             16
## 
## 
## $`miR-128-3p`
## $`miR-128-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            2982          216        0.07243461              0
##   gene.name.diff
## 1             23
## 
## 
## $`miR-140-5p`
## $`miR-140-5p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            1518           87        0.05731225              0
##   gene.name.diff
## 1             15
## 
## 
## $`miR-6821-3p`
## $`miR-6821-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            1118           93        0.08318426              0
##   gene.name.diff
## 1             12
## 
## 
## $`miR-6888-3p`
## $`miR-6888-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            2816          221        0.07848011              0
##   gene.name.diff
## 1             31
## 
## 
## $`miR-874-3p`
## $`miR-874-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            2175          180        0.08275862              0
##   gene.name.diff
## 1             15
## 
## 
## $`miR-345-3p`
## $`miR-345-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            1950          160        0.08205128              0
##   gene.name.diff
## 1             15
## 
## 
## $`miR-6072/6891-3p`
## $`miR-6072/6891-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1            1680          135        0.08035714              0
##   gene.name.diff
## 1             14
## 
## 
## $`miR-1306-3p`
## $`miR-1306-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1             264           20        0.07575758              0
##   gene.name.diff
## 1              2
## 
## 
## $`miR-184`
## $`miR-184`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1             365           25        0.06849315              0
##   gene.name.diff
## 1              5
## 
## 
## $`miR-190-5p`
## $`miR-190-5p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1             979           62        0.06332993              0
##   gene.name.diff
## 1              4
## 
## 
## $`miR-1468-5p`
## $`miR-1468-5p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1             394           20        0.05076142              0
##   gene.name.diff
## 1              2
## 
## 
## $`miR-95-3p`
## $`miR-95-3p`$summary
##   transc.id.input transc.id.na transc.id.na.prop transc.id.diff
## 1             196           14        0.07142857              0
##   gene.name.diff
## 1              3
```

## Visualize
## Save data

Save the list of unique targets for each miRNA, and the summary file 


```r
save(targets.unique, targets.rst, targets.rst.sum, file="../4_filtered_targetscan_rst.Rdata")
```

