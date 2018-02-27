**Script:** `4_filter_novel_mirna_hsa_blast_results.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts/`

**Date:**  `11/29/17`

**Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/`

**Input File(s):** 

1. `1_novel_mirna_precursor_blast_hsa_e5.txt`
2. `1_novel_mirna_mature_blast_hsa_e5.txt`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/1_blast_novel_mirna_output/`

**Output File(s):** 

1. `2_filtered_novel_mirna_precursor_abundance_e5.txt`
2. `2_novel_mirna_mature_abundance_e5.txt`
3. `3_filtered_novel_mirna_precursor_ids_e5.txt`
4. `3_filtered_novel_mirna_mature_ids_e5.txt`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to summarize the output from the BLAST query in order to characterize the candidate novel miRNAs present in the small RNA seq data.

So, need to load the precursor BLAST dataset and the full dataset of candidate novel miRNA and compare the names of sequences to see if the most abundant miRNA had BLAST results.

## Install libraries


```r
rm(list=ls())
```

## Load data


```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts")

hsa.blast.precursor<-read.table("../1_blast_novel_mirna_output/1_novel_mirna_precursor_blast_hsa_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
mature.hsa.blast<-read.table("../1_blast_novel_mirna_output/1_novel_mirna_mature_blast_hsa_e5.txt", header=FALSE, sep="", col.names=c("query_id","dbseq_id","perc_identical","length","mismatch","gapopen","query_start","query_end","dbseq_start","dbseq_end","evalue","bitscore"))
```

--------------------


```r
dim(hsa.blast.precursor)
```

```
## [1] 64 12
```

```r
head(hsa.blast.precursor)
```

```
##                                  query_id     dbseq_id perc_identical
## 1 seq|15_35462|candidatenovelpercursormiR  hsa-mir-26b         100.00
## 2  seq|X_41696|candidatenovelpercursormiR  hsa-mir-660          97.67
## 3   seq|3_8963|candidatenovelpercursormiR   hsa-mir-93          96.00
## 4   seq|3_8965|candidatenovelpercursormiR hsa-mir-106b         100.00
## 5   seq|3_9200|candidatenovelpercursormiR hsa-mir-193b         100.00
## 6   seq|3_7677|candidatenovelpercursormiR  hsa-mir-590          96.00
##   length mismatch gapopen query_start query_end dbseq_start dbseq_end
## 1     49        0       0           1        49          12        60
## 2     43        1       0           1        43          16        58
## 3     50        2       0           1        50          11        60
## 4     49        0       0           1        49          12        60
## 5     46        0       0           1        46          14        59
## 6     25        1       0           1        25          16        40
##   evalue bitscore
## 1  2e-23     97.6
## 2  2e-17     77.8
## 3  3e-19     83.8
## 4  2e-23     97.6
## 5  1e-21     91.7
## 6  1e-06     42.1
```

```r
str(hsa.blast.precursor)
```

```
## 'data.frame':	64 obs. of  12 variables:
##  $ query_id      : Factor w/ 44 levels "seq|1_1142|candidatenovelpercursormiR",..: 10 41 14 15 16 13 20 7 7 19 ...
##  $ dbseq_id      : Factor w/ 56 levels "hsa-mir-106b",..: 18 52 56 1 8 48 34 16 17 33 ...
##  $ perc_identical: num  100 97.7 96 100 100 ...
##  $ length        : int  49 43 50 49 46 25 53 51 24 55 ...
##  $ mismatch      : int  0 1 2 0 0 1 2 2 1 4 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  49 43 50 49 46 25 53 51 24 55 ...
##  $ dbseq_start   : int  12 16 11 12 14 16 4 10 14 6 ...
##  $ dbseq_end     : int  60 58 60 60 59 40 56 60 37 60 ...
##  $ evalue        : num  2e-23 2e-17 3e-19 2e-23 1e-21 ...
##  $ bitscore      : num  97.6 77.8 83.8 97.6 91.7 42.1 89.7 85.7 40.1 77.8 ...
```

```r
hsa.blast.precursor$dbseq_id<-as.character(hsa.blast.precursor$dbseq_id)
hsa.blast.precursor$query_id<-as.character(hsa.blast.precursor$query_id)
```

--------------------


```r
dim(mature.hsa.blast)
```

```
## [1] 43 12
```

```r
head(mature.hsa.blast)
```

```
##                               query_id        dbseq_id perc_identical
## 1 seq|15_35462|candidatenovelmaturemiR  hsa-miR-26b-5p         100.00
## 2  seq|X_41696|candidatenovelmaturemiR  hsa-miR-660-5p         100.00
## 3   seq|3_8963|candidatenovelmaturemiR   hsa-miR-93-5p         100.00
## 4   seq|3_8965|candidatenovelmaturemiR hsa-miR-106b-5p         100.00
## 5   seq|3_9200|candidatenovelmaturemiR hsa-miR-193b-3p          95.45
## 6   seq|3_7677|candidatenovelmaturemiR  hsa-miR-590-3p         100.00
##   length mismatch gapopen query_start query_end dbseq_start dbseq_end
## 1     21        0       0           1        21           1        21
## 2     22        0       0           1        22           1        22
## 3     23        0       0           1        23           1        23
## 4     21        0       0           1        21           1        21
## 5     22        1       0           1        22           1        22
## 6     21        0       0           1        21           1        21
##   evalue bitscore
## 1  9e-08     42.1
## 2  2e-08     44.1
## 3  6e-09     46.1
## 4  9e-08     42.1
## 5  5e-06     36.2
## 6  8e-08     42.1
```

```r
str(mature.hsa.blast)
```

```
## 'data.frame':	43 obs. of  12 variables:
##  $ query_id      : Factor w/ 38 levels "seq|1_1142|candidatenovelmaturemiR",..: 8 36 12 13 14 11 6 32 32 31 ...
##  $ dbseq_id      : Factor w/ 42 levels "hsa-miR-106b-5p",..: 14 40 42 1 6 36 13 8 9 39 ...
##  $ perc_identical: num  100 100 100 100 95.5 ...
##  $ length        : int  21 22 23 21 22 21 22 22 21 21 ...
##  $ mismatch      : int  0 0 0 0 1 0 1 0 0 0 ...
##  $ gapopen       : int  0 0 0 0 0 0 0 0 0 0 ...
##  $ query_start   : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ query_end     : int  21 22 23 21 22 21 22 22 21 21 ...
##  $ dbseq_start   : int  1 1 1 1 1 1 1 1 21 1 ...
##  $ dbseq_end     : int  21 22 23 21 22 21 22 22 1 21 ...
##  $ evalue        : num  9e-08 2e-08 6e-09 9e-08 5e-06 8e-08 5e-06 2e-08 9e-08 8e-08 ...
##  $ bitscore      : num  42.1 44.1 46.1 42.1 36.2 42.1 36.2 44.1 42.1 42.1 ...
```

```r
mature.hsa.blast$dbseq_id<-as.character(mature.hsa.blast$dbseq_id)
mature.hsa.blast$query_id<-as.character(mature.hsa.blast$query_id)
```

--------------------


```r
load("../5_novel_miRNA_filtered.Rdata")

dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 153  17
```

```r
names(novelmir10sigrandfoldmincounts)
```

```
##  [1] "provisional.id"                                                   
##  [2] "miRDeep2.score"                                                   
##  [3] "estimated.probability.that.the.miRNA.candidate.is.a.true.positive"
##  [4] "rfam.alert"                                                       
##  [5] "total.read.count"                                                 
##  [6] "mature.read.count"                                                
##  [7] "loop.read.count"                                                  
##  [8] "star.read.count"                                                  
##  [9] "significant.randfold.p.value"                                     
## [10] "miRBase.miRNA"                                                    
## [11] "example.miRBase.miRNA.with.the.same.seed"                         
## [12] "UCSC.browser"                                                     
## [13] "NCBI.blastn"                                                      
## [14] "consensus.mature.sequence"                                        
## [15] "consensus.star.sequence"                                          
## [16] "consensus.precursor.sequence"                                     
## [17] "precursor.coordinate"
```

Make the consensus.mature.sequence column into a character vector to count the length of the strings


```r
novelmir10sigrandfoldmincounts$consensus.mature.sequence<-as.character(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
rownames(novelmir10sigrandfoldmincounts)<- novelmir10sigrandfoldmincounts$provisional.id
head(novelmir10sigrandfoldmincounts)
```

```
##          provisional.id miRDeep2.score
## 15_35462       15_35462       572747.8
## X_41696         X_41696       136777.2
## 2_5946           2_5946       128796.8
## 3_8963           3_8963       117024.8
## 3_8965           3_8965       116410.3
## 3_9200           3_9200        89149.5
##          estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 15_35462                                                         91 +/- 2%
## X_41696                                                          91 +/- 2%
## 2_5946                                                           91 +/- 2%
## 3_8963                                                           91 +/- 2%
## 3_8965                                                           91 +/- 2%
## 3_9200                                                           91 +/- 2%
##          rfam.alert total.read.count mature.read.count loop.read.count
## 15_35462          -          1123411           1119830               0
## X_41696           -           268275            268073               0
## 2_5946            -           252627            249026               0
## 3_8963            -           229530            207261               0
## 3_8965            -           228324            203407               5
## 3_9200            -           174854            163337               0
##          star.read.count significant.randfold.p.value miRBase.miRNA
## 15_35462            3581                          yes             -
## X_41696              202                          yes             -
## 2_5946              3601                          yes             -
## 3_8963             22269                          yes             -
## 3_8965             24912                          yes             -
## 3_9200             11517                          yes             -
##          example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 15_35462                           hsa-miR-26a-5p            -           -
## X_41696                            hsa-miR-660-5p            -           -
## 2_5946                                          -            -           -
## 3_8963                              hsa-miR-17-5p            -           -
## 3_8965                              hsa-miR-17-5p            -           -
## 3_9200                            hsa-miR-193a-3p            -           -
##          consensus.mature.sequence consensus.star.sequence
## 15_35462    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## X_41696     uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 2_5946      uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 3_8963     caaagugcuguucgugcagguag acugcugagcuagcacuucccga
## 3_8965      uaaagugcugacagugcagaua  ccgcacuguggguacuugcugc
## 3_9200      aacuggcccacaaagucccgcu  cgggguuuugagggcgagauga
##                                            consensus.precursor.sequence
## 15_35462      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## X_41696      uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 2_5946    agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 3_8963   caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 3_8965   uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 3_9200      cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
##               precursor.coordinate
## 15_35462 15:120453419..120453476:+
## X_41696     X:43714081..43714139:+
## 2_5946        2:1474435..1474496:-
## 3_8963        3:8006929..8006991:-
## 3_8965        3:8007144..8007206:-
## 3_9200      3:28913101..28913160:-
```

```r
novelmir10sigrandfoldmincounts$consensus.mature.sequence.length<-nchar(novelmir10sigrandfoldmincounts$consensus.mature.sequence)
head(novelmir10sigrandfoldmincounts)
```

```
##          provisional.id miRDeep2.score
## 15_35462       15_35462       572747.8
## X_41696         X_41696       136777.2
## 2_5946           2_5946       128796.8
## 3_8963           3_8963       117024.8
## 3_8965           3_8965       116410.3
## 3_9200           3_9200        89149.5
##          estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 15_35462                                                         91 +/- 2%
## X_41696                                                          91 +/- 2%
## 2_5946                                                           91 +/- 2%
## 3_8963                                                           91 +/- 2%
## 3_8965                                                           91 +/- 2%
## 3_9200                                                           91 +/- 2%
##          rfam.alert total.read.count mature.read.count loop.read.count
## 15_35462          -          1123411           1119830               0
## X_41696           -           268275            268073               0
## 2_5946            -           252627            249026               0
## 3_8963            -           229530            207261               0
## 3_8965            -           228324            203407               5
## 3_9200            -           174854            163337               0
##          star.read.count significant.randfold.p.value miRBase.miRNA
## 15_35462            3581                          yes             -
## X_41696              202                          yes             -
## 2_5946              3601                          yes             -
## 3_8963             22269                          yes             -
## 3_8965             24912                          yes             -
## 3_9200             11517                          yes             -
##          example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 15_35462                           hsa-miR-26a-5p            -           -
## X_41696                            hsa-miR-660-5p            -           -
## 2_5946                                          -            -           -
## 3_8963                              hsa-miR-17-5p            -           -
## 3_8965                              hsa-miR-17-5p            -           -
## 3_9200                            hsa-miR-193a-3p            -           -
##          consensus.mature.sequence consensus.star.sequence
## 15_35462    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## X_41696     uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 2_5946      uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 3_8963     caaagugcuguucgugcagguag acugcugagcuagcacuucccga
## 3_8965      uaaagugcugacagugcagaua  ccgcacuguggguacuugcugc
## 3_9200      aacuggcccacaaagucccgcu  cgggguuuugagggcgagauga
##                                            consensus.precursor.sequence
## 15_35462      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## X_41696      uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 2_5946    agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 3_8963   caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 3_8965   uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 3_9200      cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
##               precursor.coordinate consensus.mature.sequence.length
## 15_35462 15:120453419..120453476:+                               22
## X_41696     X:43714081..43714139:+                               22
## 2_5946        2:1474435..1474496:-                               22
## 3_8963        3:8006929..8006991:-                               23
## 3_8965        3:8007144..8007206:-                               22
## 3_9200      3:28913101..28913160:-                               22
```

## Analysis

### 1. The goal is to compare the candidate novel miRNAs with BLAST results to the most abundant candidate novel miRNAs.

First, return the name of the sequences to the way they were before BLASTing at e-value = 1x10^-5.


```r
hsa.blast.precursor$seqname<-sapply(strsplit(hsa.blast.precursor$query_id, "|", fixed=TRUE),'[',2)
head(hsa.blast.precursor)
```

```
##                                  query_id     dbseq_id perc_identical
## 1 seq|15_35462|candidatenovelpercursormiR  hsa-mir-26b         100.00
## 2  seq|X_41696|candidatenovelpercursormiR  hsa-mir-660          97.67
## 3   seq|3_8963|candidatenovelpercursormiR   hsa-mir-93          96.00
## 4   seq|3_8965|candidatenovelpercursormiR hsa-mir-106b         100.00
## 5   seq|3_9200|candidatenovelpercursormiR hsa-mir-193b         100.00
## 6   seq|3_7677|candidatenovelpercursormiR  hsa-mir-590          96.00
##   length mismatch gapopen query_start query_end dbseq_start dbseq_end
## 1     49        0       0           1        49          12        60
## 2     43        1       0           1        43          16        58
## 3     50        2       0           1        50          11        60
## 4     49        0       0           1        49          12        60
## 5     46        0       0           1        46          14        59
## 6     25        1       0           1        25          16        40
##   evalue bitscore  seqname
## 1  2e-23     97.6 15_35462
## 2  2e-17     77.8  X_41696
## 3  3e-19     83.8   3_8963
## 4  2e-23     97.6   3_8965
## 5  1e-21     91.7   3_9200
## 6  1e-06     42.1   3_7677
```

Identify the unique number of sequences with BLAST hits at eval = 1x10^-5


```r
seqids<-unique(hsa.blast.precursor$seqname)
length(seqids)
```

```
## [1] 44
```

```r
seqids
```

```
##  [1] "15_35462"             "X_41696"              "3_8963"              
##  [4] "3_8965"               "3_9200"               "3_7677"              
##  [7] "AEMK02000452.1_43623" "13_28704"             "AEMK02000452.1_43616"
## [10] "AEMK02000452.1_43759" "AEMK02000452.1_43701" "AEMK02000452.1_43667"
## [13] "AEMK02000452.1_43691" "AEMK02000452.1_43697" "AEMK02000452.1_43641"
## [16] "AEMK02000452.1_43631" "AEMK02000452.1_43677" "AEMK02000452.1_43685"
## [19] "12_26305"             "AEMK02000452.1_43651" "AEMK02000452.1_43635"
## [22] "1_1142"               "X_41688"              "AEMK02000452.1_43639"
## [25] "12_27078"             "X_41792"              "7_18685"             
## [28] "13_29363"             "6_15075"              "AEMK02000452.1_43633"
## [31] "AEMK02000452.1_43688" "AEMK02000452.1_43665" "17_39198"            
## [34] "AEMK02000452.1_43794" "AEMK02000682.1_43991" "14_32585"            
## [37] "X_42552"              "1_2701"               "X_42661"             
## [40] "12_27143"             "AEMK02000452.1_43671" "AEMK02000452.1_43638"
## [43] "1_1896"               "2_7087"
```

Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset


```r
match(seqids, rownames(novelmir10sigrandfoldmincounts))
```

```
##  [1]   1   2   4   5   6   7   8   9  10  11  13  14  17  18  20  21  22
## [18]  24  25  27  28  29  30  31  33  36  39  43  48  52  60  64  68  70
## [35]  73  75  86  88  93  94 101 112 115 141
```

```r
rownames(novelmir10sigrandfoldmincounts)%in%seqids
```

```
##   [1]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [12] FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
##  [23] FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE
##  [34] FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
##  [45] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
##  [56] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
##  [67] FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE
##  [78] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE
##  [89] FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
## [100] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [111] FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
## [122] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [133] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
## [144] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```

Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.


```r
novelmir.abundance<-novelmir10sigrandfoldmincounts[seqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count", "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]
```

Is this object ordered by miRDeep2 score or by total.read.count?


```r
sum(rownames(novelmir.abundance[order(novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(novelmir.abundance))
```

```
## [1] 0
```

```r
sum(rownames(novelmir.abundance[order(novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(novelmir.abundance))
```

```
## [1] 0
```

Combine the information; are the predicted miRNA precursors with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?

This object (candidate novel precursors with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.


```r
novelmir.abundance<-cbind(novelmir.abundance[match(hsa.blast.precursor$seqname, rownames(novelmir.abundance)),], hsa.blast.precursor)
sum(novelmir.abundance$seqname != novelmir.abundance$provisional.id)
```

```
## [1] 0
```

```r
head(novelmir.abundance)
```

```
##          provisional.id miRDeep2.score total.read.count mature.read.count
## 15_35462       15_35462       572747.8          1123411           1119830
## X_41696         X_41696       136777.2           268275            268073
## 3_8963           3_8963       117024.8           229530            207261
## 3_8965           3_8965       116410.3           228324            203407
## 3_9200           3_9200        89149.5           174854            163337
## 3_7677           3_7677        57547.6           112871            104821
##          star.read.count consensus.mature.sequence
## 15_35462            3581    uucaaguaauucaggauagguu
## X_41696              202    uacccauugcauaucggaguug
## 3_8963             22269   caaagugcuguucgugcagguag
## 3_8965             24912    uaaagugcugacagugcagaua
## 3_9200             11517    aacuggcccacaaagucccgcu
## 3_7677              8050     uaauuuuauguauaagcuagu
##          consensus.mature.sequence.length
## 15_35462                               22
## X_41696                                22
## 3_8963                                 23
## 3_8965                                 22
## 3_9200                                 22
## 3_7677                                 21
##          example.miRBase.miRNA.with.the.same.seed
## 15_35462                           hsa-miR-26a-5p
## X_41696                            hsa-miR-660-5p
## 3_8963                              hsa-miR-17-5p
## 3_8965                              hsa-miR-17-5p
## 3_9200                            hsa-miR-193a-3p
## 3_7677                             hsa-miR-590-3p
##               precursor.coordinate                                query_id
## 15_35462 15:120453419..120453476:+ seq|15_35462|candidatenovelpercursormiR
## X_41696     X:43714081..43714139:+  seq|X_41696|candidatenovelpercursormiR
## 3_8963        3:8006929..8006991:-   seq|3_8963|candidatenovelpercursormiR
## 3_8965        3:8007144..8007206:-   seq|3_8965|candidatenovelpercursormiR
## 3_9200      3:28913101..28913160:-   seq|3_9200|candidatenovelpercursormiR
## 3_7677      3:11341773..11341834:+   seq|3_7677|candidatenovelpercursormiR
##              dbseq_id perc_identical length mismatch gapopen query_start
## 15_35462  hsa-mir-26b         100.00     49        0       0           1
## X_41696   hsa-mir-660          97.67     43        1       0           1
## 3_8963     hsa-mir-93          96.00     50        2       0           1
## 3_8965   hsa-mir-106b         100.00     49        0       0           1
## 3_9200   hsa-mir-193b         100.00     46        0       0           1
## 3_7677    hsa-mir-590          96.00     25        1       0           1
##          query_end dbseq_start dbseq_end evalue bitscore  seqname
## 15_35462        49          12        60  2e-23     97.6 15_35462
## X_41696         43          16        58  2e-17     77.8  X_41696
## 3_8963          50          11        60  3e-19     83.8   3_8963
## 3_8965          49          12        60  2e-23     97.6   3_8965
## 3_9200          46          14        59  1e-21     91.7   3_9200
## 3_7677          25          16        40  1e-06     42.1   3_7677
```

```r
novelmir.abundance<-novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(novelmir.abundance)
```

```
## [1] 64 10
```

```r
head(novelmir.abundance)
```

```
##          provisional.id miRDeep2.score total.read.count
## 15_35462       15_35462       572747.8          1123411
## X_41696         X_41696       136777.2           268275
## 3_8963           3_8963       117024.8           229530
## 3_8965           3_8965       116410.3           228324
## 3_9200           3_9200        89149.5           174854
## 3_7677           3_7677        57547.6           112871
##          consensus.mature.sequence consensus.mature.sequence.length
## 15_35462    uucaaguaauucaggauagguu                               22
## X_41696     uacccauugcauaucggaguug                               22
## 3_8963     caaagugcuguucgugcagguag                               23
## 3_8965      uaaagugcugacagugcagaua                               22
## 3_9200      aacuggcccacaaagucccgcu                               22
## 3_7677       uaauuuuauguauaagcuagu                               21
##          example.miRBase.miRNA.with.the.same.seed     dbseq_id
## 15_35462                           hsa-miR-26a-5p  hsa-mir-26b
## X_41696                            hsa-miR-660-5p  hsa-mir-660
## 3_8963                              hsa-miR-17-5p   hsa-mir-93
## 3_8965                              hsa-miR-17-5p hsa-mir-106b
## 3_9200                            hsa-miR-193a-3p hsa-mir-193b
## 3_7677                             hsa-miR-590-3p  hsa-mir-590
##          perc_identical evalue      precursor.coordinate
## 15_35462         100.00  2e-23 15:120453419..120453476:+
## X_41696           97.67  2e-17    X:43714081..43714139:+
## 3_8963            96.00  3e-19      3:8006929..8006991:-
## 3_8965           100.00  2e-23      3:8007144..8007206:-
## 3_9200           100.00  1e-21    3:28913101..28913160:-
## 3_7677            96.00  1e-06    3:11341773..11341834:+
```

```r
novelmir.provisional.ids<-unique(novelmir.abundance$provisional.id)
```

---------------------------------------

Now to repeat the same analysis on the candidate novel mature miRNAs with BLAST results at an e-value of 1x10-5

First, return the name of the mature sequences to the way they were before BLASTing at e-value = 1x10^-5.


```r
mature.hsa.blast$seqname<-sapply(strsplit(mature.hsa.blast$query_id, "|", fixed=TRUE),'[',2)
head(mature.hsa.blast)
```

```
##                               query_id        dbseq_id perc_identical
## 1 seq|15_35462|candidatenovelmaturemiR  hsa-miR-26b-5p         100.00
## 2  seq|X_41696|candidatenovelmaturemiR  hsa-miR-660-5p         100.00
## 3   seq|3_8963|candidatenovelmaturemiR   hsa-miR-93-5p         100.00
## 4   seq|3_8965|candidatenovelmaturemiR hsa-miR-106b-5p         100.00
## 5   seq|3_9200|candidatenovelmaturemiR hsa-miR-193b-3p          95.45
## 6   seq|3_7677|candidatenovelmaturemiR  hsa-miR-590-3p         100.00
##   length mismatch gapopen query_start query_end dbseq_start dbseq_end
## 1     21        0       0           1        21           1        21
## 2     22        0       0           1        22           1        22
## 3     23        0       0           1        23           1        23
## 4     21        0       0           1        21           1        21
## 5     22        1       0           1        22           1        22
## 6     21        0       0           1        21           1        21
##   evalue bitscore  seqname
## 1  9e-08     42.1 15_35462
## 2  2e-08     44.1  X_41696
## 3  6e-09     46.1   3_8963
## 4  9e-08     42.1   3_8965
## 5  5e-06     36.2   3_9200
## 6  8e-08     42.1   3_7677
```

Identify the unique number of sequences with BLAST hits at eval = 1x10^-5


```r
matureseqids<-unique(mature.hsa.blast$seqname)
length(matureseqids)
```

```
## [1] 38
```

```r
matureseqids
```

```
##  [1] "15_35462"             "X_41696"              "3_8963"              
##  [4] "3_8965"               "3_9200"               "3_7677"              
##  [7] "13_28704"             "AEMK02000452.1_43759" "AEMK02000452.1_43701"
## [10] "AEMK02000452.1_43667" "AEMK02000452.1_43691" "AEMK02000452.1_43697"
## [13] "2_5944"               "AEMK02000452.1_43641" "AEMK02000452.1_43631"
## [16] "AEMK02000452.1_43677" "AEMK02000452.1_43685" "12_26305"            
## [19] "AEMK02000452.1_43651" "AEMK02000452.1_43635" "1_1142"              
## [22] "X_41688"              "AEMK02000452.1_43639" "12_27078"            
## [25] "X_41792"              "7_18685"              "13_29363"            
## [28] "6_15075"              "AEMK02000452.1_43633" "AEMK02000452.1_43688"
## [31] "AEMK02000452.1_43794" "AEMK02000682.1_43991" "X_42552"             
## [34] "12_27143"             "AEMK02000452.1_43671" "AEMK02000452.1_43638"
## [37] "1_1896"               "2_7087"
```

Determine where the sequences with BLAST results at eval = 1x10^-5 match in the candidate novel miRNAs dataset


```r
match(matureseqids, rownames(novelmir10sigrandfoldmincounts))
```

```
##  [1]   1   2   4   5   6   7   9  11  13  14  17  18  19  20  21  22  24
## [18]  25  27  28  29  30  31  33  36  39  43  48  52  60  70  73  86  94
## [35] 101 112 115 141
```

```r
rownames(novelmir10sigrandfoldmincounts)%in%matureseqids
```

```
##   [1]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE
##  [12] FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [23] FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE
##  [34] FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
##  [45] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
##  [56] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
##  [67] FALSE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
##  [78] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
##  [89] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
## [100] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [111] FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
## [122] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [133] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
## [144] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```

Subset the information on the candidate miRNAs based on those with BLAST results at eval = 1x10^-5.


```r
mature.novelmir.abundance<-novelmir10sigrandfoldmincounts[matureseqids, c("provisional.id", "miRDeep2.score", "total.read.count", "mature.read.count", "star.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "precursor.coordinate")]
```

Is this object ordered by miRDeep2 score or by total.read.count?


```r
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$miRDeep2.score, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))
```

```
## [1] 0
```

```r
sum(rownames(mature.novelmir.abundance[order(mature.novelmir.abundance$total.read.count, decreasing = TRUE),])!=rownames(mature.novelmir.abundance))
```

```
## [1] 0
```

Combine the information; are the predicted miRNAs with BLAST results at eval = 1x10^-5 the most abundant predicted miRNAs?

This object (candidate novels with BLAST hits at eval = 1x10^-5) will be saved to visually inspect and summarize.


```r
mature.novelmir.abundance<-cbind(mature.novelmir.abundance[match(mature.hsa.blast$seqname, rownames(mature.novelmir.abundance)),], mature.hsa.blast)
sum(mature.novelmir.abundance$seqname != mature.novelmir.abundance$provisional.id)
```

```
## [1] 0
```

```r
head(mature.novelmir.abundance)
```

```
##          provisional.id miRDeep2.score total.read.count mature.read.count
## 15_35462       15_35462       572747.8          1123411           1119830
## X_41696         X_41696       136777.2           268275            268073
## 3_8963           3_8963       117024.8           229530            207261
## 3_8965           3_8965       116410.3           228324            203407
## 3_9200           3_9200        89149.5           174854            163337
## 3_7677           3_7677        57547.6           112871            104821
##          star.read.count consensus.mature.sequence
## 15_35462            3581    uucaaguaauucaggauagguu
## X_41696              202    uacccauugcauaucggaguug
## 3_8963             22269   caaagugcuguucgugcagguag
## 3_8965             24912    uaaagugcugacagugcagaua
## 3_9200             11517    aacuggcccacaaagucccgcu
## 3_7677              8050     uaauuuuauguauaagcuagu
##          consensus.mature.sequence.length
## 15_35462                               22
## X_41696                                22
## 3_8963                                 23
## 3_8965                                 22
## 3_9200                                 22
## 3_7677                                 21
##          example.miRBase.miRNA.with.the.same.seed
## 15_35462                           hsa-miR-26a-5p
## X_41696                            hsa-miR-660-5p
## 3_8963                              hsa-miR-17-5p
## 3_8965                              hsa-miR-17-5p
## 3_9200                            hsa-miR-193a-3p
## 3_7677                             hsa-miR-590-3p
##               precursor.coordinate                             query_id
## 15_35462 15:120453419..120453476:+ seq|15_35462|candidatenovelmaturemiR
## X_41696     X:43714081..43714139:+  seq|X_41696|candidatenovelmaturemiR
## 3_8963        3:8006929..8006991:-   seq|3_8963|candidatenovelmaturemiR
## 3_8965        3:8007144..8007206:-   seq|3_8965|candidatenovelmaturemiR
## 3_9200      3:28913101..28913160:-   seq|3_9200|candidatenovelmaturemiR
## 3_7677      3:11341773..11341834:+   seq|3_7677|candidatenovelmaturemiR
##                 dbseq_id perc_identical length mismatch gapopen
## 15_35462  hsa-miR-26b-5p         100.00     21        0       0
## X_41696   hsa-miR-660-5p         100.00     22        0       0
## 3_8963     hsa-miR-93-5p         100.00     23        0       0
## 3_8965   hsa-miR-106b-5p         100.00     21        0       0
## 3_9200   hsa-miR-193b-3p          95.45     22        1       0
## 3_7677    hsa-miR-590-3p         100.00     21        0       0
##          query_start query_end dbseq_start dbseq_end evalue bitscore
## 15_35462           1        21           1        21  9e-08     42.1
## X_41696            1        22           1        22  2e-08     44.1
## 3_8963             1        23           1        23  6e-09     46.1
## 3_8965             1        21           1        21  9e-08     42.1
## 3_9200             1        22           1        22  5e-06     36.2
## 3_7677             1        21           1        21  8e-08     42.1
##           seqname
## 15_35462 15_35462
## X_41696   X_41696
## 3_8963     3_8963
## 3_8965     3_8965
## 3_9200     3_9200
## 3_7677     3_7677
```

```r
mature.novelmir.abundance<-mature.novelmir.abundance[,c("provisional.id", "miRDeep2.score", "total.read.count",  "consensus.mature.sequence", "consensus.mature.sequence.length", "example.miRBase.miRNA.with.the.same.seed", "dbseq_id", "perc_identical", "evalue", "precursor.coordinate")]
dim(mature.novelmir.abundance)
```

```
## [1] 43 10
```

```r
head(mature.novelmir.abundance)
```

```
##          provisional.id miRDeep2.score total.read.count
## 15_35462       15_35462       572747.8          1123411
## X_41696         X_41696       136777.2           268275
## 3_8963           3_8963       117024.8           229530
## 3_8965           3_8965       116410.3           228324
## 3_9200           3_9200        89149.5           174854
## 3_7677           3_7677        57547.6           112871
##          consensus.mature.sequence consensus.mature.sequence.length
## 15_35462    uucaaguaauucaggauagguu                               22
## X_41696     uacccauugcauaucggaguug                               22
## 3_8963     caaagugcuguucgugcagguag                               23
## 3_8965      uaaagugcugacagugcagaua                               22
## 3_9200      aacuggcccacaaagucccgcu                               22
## 3_7677       uaauuuuauguauaagcuagu                               21
##          example.miRBase.miRNA.with.the.same.seed        dbseq_id
## 15_35462                           hsa-miR-26a-5p  hsa-miR-26b-5p
## X_41696                            hsa-miR-660-5p  hsa-miR-660-5p
## 3_8963                              hsa-miR-17-5p   hsa-miR-93-5p
## 3_8965                              hsa-miR-17-5p hsa-miR-106b-5p
## 3_9200                            hsa-miR-193a-3p hsa-miR-193b-3p
## 3_7677                             hsa-miR-590-3p  hsa-miR-590-3p
##          perc_identical evalue      precursor.coordinate
## 15_35462         100.00  9e-08 15:120453419..120453476:+
## X_41696          100.00  2e-08    X:43714081..43714139:+
## 3_8963           100.00  6e-09      3:8006929..8006991:-
## 3_8965           100.00  9e-08      3:8007144..8007206:-
## 3_9200            95.45  5e-06    3:28913101..28913160:-
## 3_7677           100.00  8e-08    3:11341773..11341834:+
```

```r
mature.novelmir.provisional.ids<-unique(mature.novelmir.abundance$provisional.id)
```

## Save data


```r
write.table(novelmir.abundance, file="../1_blast_novel_mirna_output/2_filtered_novel_mirna_precursor_abundance_e5.txt")
write.table(mature.novelmir.abundance, file="../1_blast_novel_mirna_output/2_novel_mirna_mature_abundance_e5.txt")
```

Create a list of the pertinent precursor provisional.ids to extract the correct pdfs for a supplemental figure


```r
write.table(novelmir.provisional.ids, file="../1_blast_novel_mirna_output/3_filtered_novel_mirna_precursor_ids_e5.txt", row.names=FALSE, col.names=FALSE)
write.table(mature.novelmir.provisional.ids, file="../1_blast_novel_mirna_output/3_filtered_novel_mirna_mature_ids_e5.txt", row.names=FALSE, col.names=FALSE)
```

