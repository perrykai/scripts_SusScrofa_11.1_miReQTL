**Script:** `2_novel_mirna_seq_analysis.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`

**Date:**  11/29/17

**Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`

**Input File(s):** `2_predicted_novel_mirna.csv`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`

**Output File(s):** 

1. `5_novel_miRNA_filtered.Rdata`
2. `6_novel_mature_mir.fa`
3. `7_novel_precursor_mir.fa`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

The objective of this script is to filter the extracted putative novel miRNAs for the following characteristics:
 
1. Novel miRNA candidate must have > 90% probability of being a true positive (filter by miRDeep2 score and then estimated probability that the miRNA candidate is a true positive)
2. Hairpins must have significant Randfold p-values
3. Minimum read counts for putative mature and star strand sequences

Additionally, the novel candidate miRNAs will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)

## Install libraries


```r
library("ShortRead")
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: methods
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     grep, grepl, intersect, is.unsorted, lapply, lengths, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unlist, unsplit
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: Biostrings
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: XVector
```

```
## Loading required package: Rsamtools
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomicAlignments
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts")
```

## Load data



```r
novelmir<-read.table("../2_predicted_novel_mirna.csv", sep="\t", header=TRUE)
dim(novelmir)
```

```
## [1] 1201   17
```

```r
colnames(novelmir)
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

```r
head(novelmir)
```

```
##   provisional.id miRDeep2.score
## 1       15_35462       572747.8
## 2        X_41696       136777.2
## 3         2_5946       128796.8
## 4         3_8963       117024.8
## 5         3_8965       116410.3
## 6         3_9200        89149.5
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 4                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123411           1119830               0
## 2          -           268275            268073               0
## 3          -           252627            249026               0
## 4          -           229530            207261               0
## 5          -           228324            203407               5
## 6          -           174854            163337               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 4           22269                          yes             -
## 5           24912                          yes             -
## 6           11517                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 4                            hsa-miR-17-5p            -           -
## 5                            hsa-miR-17-5p            -           -
## 6                          hsa-miR-193a-3p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 4   caaagugcuguucgugcagguag acugcugagcuagcacuucccga
## 5    uaaagugcugacagugcagaua  ccgcacuguggguacuugcugc
## 6    aacuggcccacaaagucccgcu  cgggguuuugagggcgagauga
##                                     consensus.precursor.sequence
## 1      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2     uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3  agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 4 caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 5 uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 6    cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
##        precursor.coordinate
## 1 15:120453419..120453476:+
## 2    X:43714081..43714139:+
## 3      2:1474435..1474496:-
## 4      3:8006929..8006991:-
## 5      3:8007144..8007206:-
## 6    3:28913101..28913160:-
```

## Analysis

### 1. Filter by the miRDeep2 score and the estimated probability of the novel miRNA being true positives:


```r
sum(novelmir$miRDeep2.score >= 10)
```

```
## [1] 369
```

```r
table(novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)
```

```
## 
## 19 +/- 3% 50 +/- 2% 52 +/- 2% 56 +/- 2% 79 +/- 2% 91 +/- 1% 91 +/- 2% 
##       155       125        33       311        46       112       419
```

```r
novelmir[novelmir$estimated.probability.that.the.miRNA.candidate.is.a.true.positive == "91 +/- 1%", "miRDeep2.score"]
```

```
##   [1] 5.9 5.9 5.9 5.9 5.9 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8
##  [18] 5.7 5.7 5.7 5.7 5.7 5.7 5.7 5.7 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6 5.6
##  [35] 5.6 5.6 5.6 5.6 5.6 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.5 5.4
##  [52] 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.3 5.3 5.3 5.3 5.3 5.3
##  [69] 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.3 5.2 5.2 5.2 5.2 5.2
##  [86] 5.2 5.2 5.2 5.2 5.2 5.2 5.2 5.2 5.1 5.1 5.1 5.1 5.1 5.1 5.1 5.1 5.1
## [103] 5.1 5.1 5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0
```

First filter by miRDeep2 score, then by estimated probability:


```r
novelmir10<-novelmir[novelmir$miRDeep2.score >= 10,]
dim(novelmir10)
```

```
## [1] 369  17
```

```r
table(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive)
```

```
## 
## 19 +/- 3% 50 +/- 2% 52 +/- 2% 56 +/- 2% 79 +/- 2% 91 +/- 1% 91 +/- 2% 
##         0         0         0         0         0         0       369
```

```r
sum(novelmir10$estimated.probability.that.the.miRNA.candidate.is.a.true.positive != "91 +/- 2%")
```

```
## [1] 0
```

Filtering by miRDeep2 score removed any miRNA with estimated probability < "91 +/- 2%"

### 2. Hairpins must have a significant Randfold p-value


```r
sum(novelmir10$significant.randfold.p.value == "yes")
```

```
## [1] 283
```

271 potential miRNA candidates have significant Randfold p-value


```r
novelmir10sigrandfold<-novelmir10[novelmir10$significant.randfold.p.value == "yes",]
dim(novelmir10sigrandfold)
```

```
## [1] 283  17
```

```r
head(novelmir10sigrandfold)
```

```
##   provisional.id miRDeep2.score
## 1       15_35462       572747.8
## 2        X_41696       136777.2
## 3         2_5946       128796.8
## 4         3_8963       117024.8
## 5         3_8965       116410.3
## 6         3_9200        89149.5
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 4                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123411           1119830               0
## 2          -           268275            268073               0
## 3          -           252627            249026               0
## 4          -           229530            207261               0
## 5          -           228324            203407               5
## 6          -           174854            163337               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 4           22269                          yes             -
## 5           24912                          yes             -
## 6           11517                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 4                            hsa-miR-17-5p            -           -
## 5                            hsa-miR-17-5p            -           -
## 6                          hsa-miR-193a-3p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 4   caaagugcuguucgugcagguag acugcugagcuagcacuucccga
## 5    uaaagugcugacagugcagaua  ccgcacuguggguacuugcugc
## 6    aacuggcccacaaagucccgcu  cgggguuuugagggcgagauga
##                                     consensus.precursor.sequence
## 1      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2     uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3  agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 4 caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 5 uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 6    cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
##        precursor.coordinate
## 1 15:120453419..120453476:+
## 2    X:43714081..43714139:+
## 3      2:1474435..1474496:-
## 4      3:8006929..8006991:-
## 5      3:8007144..8007206:-
## 6    3:28913101..28913160:-
```

Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences


```r
table(novelmir10sigrandfold$rfam.alert)
```

```
## 
##         -      rRNA rRNA/tRNA 
##       283         0         0
```

```r
sum(novelmir10sigrandfold$rfam.alert != "-")
```

```
## [1] 0
```

### 3. Minimum read counts for putative mature and star strand sequences: require at least 10 counts for each


```r
summary(novelmir10sigrandfold$mature.read.count)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##       1.0      26.5     102.0   10420.0     563.5 1120000.0
```

```r
summary(novelmir10sigrandfold$star.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     1.0     4.0    15.0   413.8    52.5 24910.0
```

```r
sum(novelmir10sigrandfold$mature.read.count <= 10)
```

```
## [1] 17
```

```r
novelmir10sigrandfoldmincounts<-novelmir10sigrandfold[novelmir10sigrandfold$mature.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 266  17
```

```r
sum(novelmir10sigrandfoldmincounts$star.read.count <= 10)
```

```
## [1] 113
```

```r
novelmir10sigrandfoldmincounts<-novelmir10sigrandfoldmincounts[novelmir10sigrandfoldmincounts$star.read.count > 10,]
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 153  17
```

```r
sum(novelmir10sigrandfoldmincounts$star.read.count <=10)
```

```
## [1] 0
```

```r
sum(novelmir10sigrandfoldmincounts$mature.read.count <=10)
```

```
## [1] 0
```

Investigate the final output:


```r
dim(novelmir10sigrandfoldmincounts)
```

```
## [1] 153  17
```

```r
head(novelmir10sigrandfoldmincounts)
```

```
##   provisional.id miRDeep2.score
## 1       15_35462       572747.8
## 2        X_41696       136777.2
## 3         2_5946       128796.8
## 4         3_8963       117024.8
## 5         3_8965       116410.3
## 6         3_9200        89149.5
##   estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1                                                         91 +/- 2%
## 2                                                         91 +/- 2%
## 3                                                         91 +/- 2%
## 4                                                         91 +/- 2%
## 5                                                         91 +/- 2%
## 6                                                         91 +/- 2%
##   rfam.alert total.read.count mature.read.count loop.read.count
## 1          -          1123411           1119830               0
## 2          -           268275            268073               0
## 3          -           252627            249026               0
## 4          -           229530            207261               0
## 5          -           228324            203407               5
## 6          -           174854            163337               0
##   star.read.count significant.randfold.p.value miRBase.miRNA
## 1            3581                          yes             -
## 2             202                          yes             -
## 3            3601                          yes             -
## 4           22269                          yes             -
## 5           24912                          yes             -
## 6           11517                          yes             -
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                           hsa-miR-26a-5p            -           -
## 2                           hsa-miR-660-5p            -           -
## 3                                        -            -           -
## 4                            hsa-miR-17-5p            -           -
## 5                            hsa-miR-17-5p            -           -
## 6                          hsa-miR-193a-3p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uucaaguaauucaggauagguu  ccuguucuccauuacuuggcuc
## 2    uacccauugcauaucggaguug  accuccuaugugcaugguuuac
## 3    uggugccugacgucuuggcagu agccagggcugcaggcacugaca
## 4   caaagugcuguucgugcagguag acugcugagcuagcacuucccga
## 5    uaaagugcugacagugcagaua  ccgcacuguggguacuugcugc
## 6    aacuggcccacaaagucccgcu  cgggguuuugagggcgagauga
##                                     consensus.precursor.sequence
## 1      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2     uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3  agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 4 caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 5 uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 6    cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
##        precursor.coordinate
## 1 15:120453419..120453476:+
## 2    X:43714081..43714139:+
## 3      2:1474435..1474496:-
## 4      3:8006929..8006991:-
## 5      3:8007144..8007206:-
## 6    3:28913101..28913160:-
```

How many of the novel miRNA candidates have a human miRBase miRNA with the same seed sequence


```r
sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed != "-")
```

```
## [1] 62
```

```r
sum(novelmir10sigrandfoldmincounts$example.miRBase.miRNA.with.the.same.seed == "-")
```

```
## [1] 91
```

What is the summary of the mature and star read counts now?


```r
summary(novelmir10sigrandfoldmincounts$mature.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      11      81     300   18980    2615 1120000
```

```r
summary(novelmir10sigrandfoldmincounts$star.read.count)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    11.0    21.0    49.0   760.7   150.0 24910.0
```

Subset the provisional.id, consensus.mature.sequence and consensus.precursor.sequence for use in BLASTN against known human and mouse miRBase sequences:


```r
novelmircandidateBLAST<-novelmir10sigrandfoldmincounts[,c("provisional.id", "consensus.mature.sequence", "consensus.precursor.sequence")]
dim(novelmircandidateBLAST)
```

```
## [1] 153   3
```

```r
head(novelmircandidateBLAST)
```

```
##   provisional.id consensus.mature.sequence
## 1       15_35462    uucaaguaauucaggauagguu
## 2        X_41696    uacccauugcauaucggaguug
## 3         2_5946    uggugccugacgucuuggcagu
## 4         3_8963   caaagugcuguucgugcagguag
## 5         3_8965    uaaagugcugacagugcagaua
## 6         3_9200    aacuggcccacaaagucccgcu
##                                     consensus.precursor.sequence
## 1      uucaaguaauucaggauagguugugugcuguccagccuguucuccauuacuuggcuc
## 2     uacccauugcauaucggaguugugaauucucaaagcaccuccuaugugcaugguuuac
## 3  agccagggcugcaggcacugacauucacccaugguauuguggugccugacgucuuggcagu
## 4 caaagugcuguucgugcagguagugugauaaccuaaccuacugcugagcuagcacuucccga
## 5 uaaagugcugacagugcagauagugguccucuccgugcuaccgcacuguggguacuugcugc
## 6    cgggguuuugagggcgagaugaguuuauguuuuauccaacuggcccacaaagucccgcu
```

### 4. The novel candidate miRNAs (both precursor and mature) will be written into fasta file format for BLAST against human and mouse miRBase databases (mature seq and precursor seq separately)

#### First, prepare the sequence names for the candidate novel mature miR


```r
seqnamesmature<-paste("seq|",novelmircandidateBLAST$provisional.id, "|candidate novel mature miR")
seqnamesmature<-gsub(" ", "", seqnamesmature, fixed = TRUE)
head(seqnamesmature)
```

```
## [1] "seq|15_35462|candidatenovelmaturemiR"
## [2] "seq|X_41696|candidatenovelmaturemiR" 
## [3] "seq|2_5946|candidatenovelmaturemiR"  
## [4] "seq|3_8963|candidatenovelmaturemiR"  
## [5] "seq|3_8965|candidatenovelmaturemiR"  
## [6] "seq|3_9200|candidatenovelmaturemiR"
```

Create the BStringSet object with the mature sequence names


```r
matureids<-BStringSet(seqnamesmature)
head(matureids)
```

```
##   A BStringSet instance of length 6
##     width seq
## [1]    36 seq|15_35462|candidatenovelmaturemiR
## [2]    35 seq|X_41696|candidatenovelmaturemiR
## [3]    34 seq|2_5946|candidatenovelmaturemiR
## [4]    34 seq|3_8963|candidatenovelmaturemiR
## [5]    34 seq|3_8965|candidatenovelmaturemiR
## [6]    34 seq|3_9200|candidatenovelmaturemiR
```

Prepare the candidate novel mature sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet


```r
novelmircandidateBLAST$consensus.mature.seq.tadjust <- gsub("u","t", novelmircandidateBLAST$consensus.mature.sequence)
```

Create the DNAStringSet object with the candidate novel mature sequence reads


```r
novelmatureseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.mature.seq.tadjust)
head(novelmatureseqreads)
```

```
##   A DNAStringSet instance of length 6
##     width seq
## [1]    22 TTCAAGTAATTCAGGATAGGTT
## [2]    22 TACCCATTGCATATCGGAGTTG
## [3]    22 TGGTGCCTGACGTCTTGGCAGT
## [4]    23 CAAAGTGCTGTTCGTGCAGGTAG
## [5]    22 TAAAGTGCTGACAGTGCAGATA
## [6]    22 AACTGGCCCACAAAGTCCCGCT
```

Use the ShortRead command to combine the candidate novel mature sequence IDs and the sequences


```r
candidatenovelmaturefasta<-ShortRead(sread=novelmatureseqreads ,id=matureids)
```

#### Then, prepare the sequence names for the candidate novel precursor miR:


```r
seqnamesprecursor<-paste("seq|", novelmircandidateBLAST$provisional.id, "|candidate novel percursor miR")
seqnamesprecursor<-gsub(" ", "", seqnamesprecursor, fixed = TRUE)

head(seqnamesprecursor)
```

```
## [1] "seq|15_35462|candidatenovelpercursormiR"
## [2] "seq|X_41696|candidatenovelpercursormiR" 
## [3] "seq|2_5946|candidatenovelpercursormiR"  
## [4] "seq|3_8963|candidatenovelpercursormiR"  
## [5] "seq|3_8965|candidatenovelpercursormiR"  
## [6] "seq|3_9200|candidatenovelpercursormiR"
```

Create the BStringSet object with the precursor sequence names


```r
precursorids<-BStringSet(seqnamesprecursor)
head(precursorids)
```

```
##   A BStringSet instance of length 6
##     width seq
## [1]    39 seq|15_35462|candidatenovelpercursormiR
## [2]    38 seq|X_41696|candidatenovelpercursormiR
## [3]    37 seq|2_5946|candidatenovelpercursormiR
## [4]    37 seq|3_8963|candidatenovelpercursormiR
## [5]    37 seq|3_8965|candidatenovelpercursormiR
## [6]    37 seq|3_9200|candidatenovelpercursormiR
```

Prepare the candidate precursor sequences for compatability with DNAStringSet: Substitute "t" for "u" as DNAStringSet only uses DNA alphabet


```r
novelmircandidateBLAST$consensus.precursor.seq.tadjust<- gsub("u","t", novelmircandidateBLAST$consensus.precursor.sequence)
```

Create the DNAStringSet object with the novel candidate precursor sequence reads


```r
novelprecursorseqreads<-DNAStringSet(novelmircandidateBLAST$consensus.precursor.seq.tadjust)
head(novelprecursorseqreads)
```

```
##   A DNAStringSet instance of length 6
##     width seq
## [1]    57 TTCAAGTAATTCAGGATAGGTTGTGTGCTGTCCAGCCTGTTCTCCATTACTTGGCTC
## [2]    58 TACCCATTGCATATCGGAGTTGTGAATTCTCAAAGCACCTCCTATGTGCATGGTTTAC
## [3]    61 AGCCAGGGCTGCAGGCACTGACATTCACCCATGGTATTGTGGTGCCTGACGTCTTGGCAGT
## [4]    62 CAAAGTGCTGTTCGTGCAGGTAGTGTGATAACCTAACCTACTGCTGAGCTAGCACTTCCCGA
## [5]    62 TAAAGTGCTGACAGTGCAGATAGTGGTCCTCTCCGTGCTACCGCACTGTGGGTACTTGCTGC
## [6]    59 CGGGGTTTTGAGGGCGAGATGAGTTTATGTTTTATCCAACTGGCCCACAAAGTCCCGCT
```

Use the ShortRead command to combine the candidate novel precursor sequence IDs and the sequences


```r
candidatenovelprecursorfasta<-ShortRead(sread=novelprecursorseqreads ,id=precursorids)
```

## Visualize

## Save data



```r
save(novelmir10sigrandfoldmincounts, file="../5_novel_miRNA_filtered.Rdata")
writeFasta(candidatenovelmaturefasta, file="../6_novel_mature_mir.fa")
writeFasta(candidatenovelprecursorfasta, file="../7_novel_precursor_mir.fa")
```

