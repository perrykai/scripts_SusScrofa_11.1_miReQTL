**Script:** `1_extract_mirna_from_core.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts`

**Date:**  11/29/17

**Input File Directory:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/1_miRNA_map_quantify/2_quant_predict_output`

**Input File(s):**  `result_27_11_2017_t_18_08_27.csv`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp`

**Output File(s):** 

1. `1_core_output_stats.csv`

2. `2_predicted_novel_mirna.csv`

3. `3_mature_mirna_detected.csv`


**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

This code extracts the miRDeep2 core module output from the `result_27_11_2017_t_18_08_27.csv` file.
This .csv file contains the miRDeep2 score distribution for the novel miRNA prediction step, 
the novel miRNAs predicted by miRDeep2, the mature miRBase miRNAs detected by miRDeep2, and the 
miRBase miRNAs not detected by miRDeep2. The first three of these items are extracted from this .csv file using this script.
The objective here is to characterize the known and novel miRNAs present in this dataset, isolate the sequences meeting miRDeep2 score of 10 or greater,
to isolate the sequences at that score cutoff having homologous seed sequence with a human miRBase miRNA,
and to estimate the false discovery rate of the miRDeep2 prediction step at miRDeep2 score 10. 

## Install libraries



```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/2_mirna_char_exp/scripts/")
```

## Load Data

To isolate the first section of the csv containing the miRDeep2 distribution scores:


```r
sts<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", nrows=21, header = TRUE)

head(sts)
```

```
##   miRDeep2.score novel.miRNAs.reported.by.miRDeep2
## 1             10                               369
## 2              9                               376
## 3              8                               391
## 4              7                               406
## 5              6                               419
## 6              5                               531
##   novel.miRNAs..estimated.false.positives
## 1                                34 +/- 6
## 2                                35 +/- 6
## 3                                36 +/- 6
## 4                                37 +/- 6
## 5                                39 +/- 6
## 6                                49 +/- 7
##   novel.miRNAs..estimated.true.positives known.miRNAs.in.species
## 1                  335 +/- 6 (91 +/- 2%)                     411
## 2                  341 +/- 6 (91 +/- 2%)                     411
## 3                  355 +/- 6 (91 +/- 2%)                     411
## 4                  369 +/- 6 (91 +/- 2%)                     411
## 5                  380 +/- 6 (91 +/- 2%)                     411
## 6                  482 +/- 7 (91 +/- 1%)                     411
##   known.miRNAs.in.data known.miRNAs.detected.by.miRDeep2
## 1                  333                         274 (82%)
## 2                  333                         274 (82%)
## 3                  333                         274 (82%)
## 4                  333                         274 (82%)
## 5                  333                         274 (82%)
## 6                  333                         300 (90%)
##   estimated.signal.to.noise excision.gearing
## 1                      15.3                4
## 2                      15.0                4
## 3                      15.0                4
## 4                      15.0                4
## 5                      14.7                4
## 6                      13.8                4
```

```r
tail(sts)
```

```
##    miRDeep2.score novel.miRNAs.reported.by.miRDeep2
## 16             -5                              2092
## 17             -6                              2205
## 18             -7                              2306
## 19             -8                              2412
## 20             -9                              2499
## 21            -10                              2572
##    novel.miRNAs..estimated.false.positives
## 16                             2136 +/- 45
## 17                             2280 +/- 46
## 18                             2377 +/- 47
## 19                             2459 +/- 48
## 20                             2519 +/- 47
## 21                             2570 +/- 47
##    novel.miRNAs..estimated.true.positives known.miRNAs.in.species
## 16                    4 +/- 12 (0 +/- 1%)                     411
## 17                     1 +/- 6 (0 +/- 0%)                     411
## 18                     2 +/- 9 (0 +/- 0%)                     411
## 19                    4 +/- 15 (0 +/- 1%)                     411
## 20                   11 +/- 23 (0 +/- 1%)                     411
## 21                   19 +/- 30 (1 +/- 1%)                     411
##    known.miRNAs.in.data known.miRNAs.detected.by.miRDeep2
## 16                  333                         312 (94%)
## 17                  333                         312 (94%)
## 18                  333                         312 (94%)
## 19                  333                         312 (94%)
## 20                  333                         312 (94%)
## 21                  333                         312 (94%)
##    estimated.signal.to.noise excision.gearing
## 16                       1.1                4
## 17                       1.1                4
## 18                       1.1                4
## 19                       1.1                4
## 20                       1.1                4
## 21                       1.1                4
```

To isolate the second section of the csv containing the novel predicted miRNAs:


```r
novelmirna<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", skip=26, nrow=1201,  header = TRUE, fill=TRUE)

head(novelmirna)
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

```r
tail(novelmirna)
```

```
##      provisional.id miRDeep2.score
## 1196       17_39309              0
## 1197        8_20276              0
## 1198       18_40535              0
## 1199         2_7487              0
## 1200       14_33505              0
## 1201        4_11082              0
##      estimated.probability.that.the.miRNA.candidate.is.a.true.positive
## 1196                                                         19 +/- 3%
## 1197                                                         19 +/- 3%
## 1198                                                         19 +/- 3%
## 1199                                                         19 +/- 3%
## 1200                                                         19 +/- 3%
## 1201                                                         19 +/- 3%
##      rfam.alert total.read.count mature.read.count loop.read.count
## 1196          -               17                17               0
## 1197          -               15                15               0
## 1198          -               47                47               0
## 1199          -               41                41               0
## 1200          -                8                 8               0
## 1201          -             1201              1201               0
##      star.read.count significant.randfold.p.value miRBase.miRNA
## 1196               0                           no             -
## 1197               0                           no             -
## 1198               0                           no             -
## 1199               0                           no             -
## 1200               0                           no             -
## 1201               0                          yes             -
##      example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1196                          hsa-miR-5009-5p            -           -
## 1197                              hsa-miR-484            -           -
## 1198                           hsa-miR-23a-5p            -           -
## 1199                              hsa-miR-608            -           -
## 1200                              hsa-miR-346            -           -
## 1201                                        -            -           -
##      consensus.mature.sequence consensus.star.sequence
## 1196        cuggacuuugaguuagaa      cuggccucauguucacuc
## 1197        ucaggcucuugggaccuc     ggaaaggagcagagugacg
## 1198        uggguuccugggggagua      ccucucagugucuugagc
## 1199        gggggugggguuggggca         gccuucugcguuugc
## 1200   ugucugcccgcaugccugccucu  aggcaggggcugggccugcagc
## 1201       uauguuaagaaguauguau     auguauuuaaauacauugg
##                                                   consensus.precursor.sequence
## 1196                                    cuggccucauguucacuccccuggacuuugaguuagaa
## 1197 ucaggcucuugggaccucccagugcaguggaggcaucagaaacgguaaacggggggaaaggagcagagugacg
## 1198                    ccucucagugucuugagcuucuccuagcaugguggcuggguuccugggggagua
## 1199                               gggggugggguuggggcagggauccuuggccuucugcguuugc
## 1200              ugucugcccgcaugccugccucucuguggcucugaaggaggcaggggcugggccugcagc
## 1201                            uauguuaagaaguauguaugugcauguauguauuuaaauacauugg
##          precursor.coordinate
## 1196  17:21714192..21714230:+
## 1197     8:9124509..9124582:+
## 1198  18:29611958..29612012:+
## 1199 2:149820991..149821034:-
## 1200  14:87144360..87144420:-
## 1201 4:124551930..124551976:+
```

To isolate the third section of the csv containing the miRBase miRNAs detected by miRDeep2:


```r
md<-read.table("../../1_miRNA_map_quantify/2_quant_predict_output/result_27_11_2017_t_18_08_27.csv", sep = "\t", skip = 1232, nrow = 289, header = TRUE, fill = TRUE)

head(md)
```

```
##     tag.id miRDeep2.score
## 1  6_17172       21446006
## 2 17_39671       21434705
## 3  6_17170        4844823
## 4 17_39683        4841198
## 5   2_5843        4659925
## 6  7_18211        4053377
##   estimated.probability.that.the.miRNA.is.a.true.positive rfam.alert
## 1                                               91 +/- 2%          -
## 2                                               91 +/- 2%          -
## 3                                               91 +/- 2%          -
## 4                                               91 +/- 2%          -
## 5                                               91 +/- 2%          -
## 6                                               91 +/- 2%          -
##   total.read.count mature.read.count loop.read.count star.read.count
## 1         42065408          42041498              96           23814
## 2         42043241          42031779            1981            9481
## 3          9502910           9103050             448          399412
## 4          9495799           9103022             282          392495
## 5          9140232           9095069               0           45163
## 6          7950516           7948388               2            2126
##   significant.randfold.p.value mature.miRBase.miRNA
## 1                          yes            ssc-miR-1
## 2                          yes            ssc-miR-1
## 3                          yes      ssc-miR-133a-5p
## 4                          yes      ssc-miR-133a-5p
## 5                          yes          ssc-miR-378
## 6                          yes          ssc-miR-206
##   example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 1                             hsa-miR-1-3p            -           -
## 2                             hsa-miR-1-3p            -           -
## 3                                        -            -           -
## 4                                        -            -           -
## 5                          hsa-miR-378a-3p            -           -
## 6                             hsa-miR-1-3p            -           -
##   consensus.mature.sequence consensus.star.sequence
## 1    uggaauguaaagaaguauguau acauacuucuuuauguacccaua
## 2    uggaauguaaagaaguauguau acauacuucuuuaugugcccaua
## 3    uugguccccuucaaccagcugu  agcugguaaaauggaaccaaau
## 4    uugguccccuucaaccagcugu  agcugguaaaauggaaccaaau
## 5    acuggacuuggagucagaaggc  cuccugacuccagguccugugu
## 6    uggaauguaaggaaguguguga acaugcuucuuuauauccccaua
##                                    consensus.precursor.sequence
## 1 acauacuucuuuauguacccauaugaacauacaaugcuauggaauguaaagaaguauguau
## 2 acauacuucuuuaugugcccauauggaccugcuaagcuauggaauguaaagaaguauguau
## 3  agcugguaaaauggaaccaaaucgccucuucaauggauuugguccccuucaaccagcugu
## 4  agcugguaaaauggaaccaaaucaacuguugaauggauuugguccccuucaaccagcugu
## 5  cuccugacuccagguccuguguguugccuggaaauagcacuggacuuggagucagaaggc
## 6  acaugcuucuuuauauccccauacggauuacuugacuauggaauguaaggaaguguguga
##       precursor.coordinate
## 1 6:107002176..107002237:-
## 2  17:61894059..61894120:+
## 3 6:106998908..106998968:-
## 4  17:61904802..61904862:+
## 5 2:150825874..150825934:+
## 6   7:45976303..45976363:+
```

```r
tail(md)
```

```
##       tag.id miRDeep2.score
## 284  4_11543           -3.7
## 285  X_42869           -4.0
## 286  8_21909           -5.4
## 287   3_8157           -6.0
## 288 10_24814           -7.3
## 289  7_18212           -8.3
##     estimated.probability.that.the.miRNA.is.a.true.positive rfam.alert
## 284                                                2 +/- 2%          -
## 285                                                2 +/- 2%          -
## 286                                                0 +/- 0%          -
## 287                                                0 +/- 0%          -
## 288                                                0 +/- 1%          -
## 289                                                0 +/- 1%          -
##     total.read.count mature.read.count loop.read.count star.read.count
## 284               12                12               0               0
## 285             1616              1616               0               0
## 286            33484             33484               0               0
## 287          3618856           3618854               2               0
## 288               75                75               0               0
## 289          7948388           7948388               0               0
##     significant.randfold.p.value mature.miRBase.miRNA
## 284                           no         ssc-miR-124a
## 285                           no         ssc-miR-2483
## 286                           no      ssc-miR-9843-3p
## 287                           no           ssc-let-7f
## 288                           no          ssc-miR-215
## 289                           no          ssc-miR-206
##     example.miRBase.miRNA.with.the.same.seed UCSC.browser NCBI.blastn
## 284                                        -            -           -
## 285                                        -            -           -
## 286                                        -            -           -
## 287                            hsa-let-7a-5p            -           -
## 288                                        -            -           -
## 289                             hsa-miR-1-3p            -           -
##     consensus.mature.sequence    consensus.star.sequence
## 284    uuaaggcacgcggugaaugcca       gacacuagagcgucuuagaa
## 285    aaacaucugguugguugagaga     uacuuagucucuggauguuuug
## 286    ucugugaacuagaaaccucugg    agucaaguuucagcucacuuaac
## 287    ugagguaguagauuguauaguu         uugugcucuauucccgaa
## 288    ugaccuaugaauugacagacaa gucauuguuaagagaugguguacagg
## 289    uggaauguaaggaaguguguga     aagcgccucuccgcggccccaa
##                                                    consensus.precursor.sequence
## 284      uuaaggcacgcggugaaugccaagagcggagccuacggcugcacuugaaggacacuagagcgucuuagaa
## 285                          aaacaucugguugguugagagaauuuuuuacuuagucucuggauguuuug
## 286 ucugugaacuagaaaccucuggaaagugguguacauuuuacuuccacccacaagucaaguuucagcucacuuaac
## 287               uugugcucuauucccgaagaagagauugcucuaucagagugagguaguagauuguauaguu
## 288                        gucauuguuaagagaugguguacaggaaaaugaccuaugaauugacagacaa
## 289                         uggaauguaaggaagugugugaucucguuaagcgccucuccgcggccccaa
##         precursor.coordinate
## 284   4:69968831..69968901:-
## 285 X:102636070..102636120:-
## 286 8:114110616..114110691:-
## 287   3:43466004..43466065:+
## 288    10:9704659..9704711:-
## 289   7:45976341..45976392:+
```

## Analysis

Now I can open the `novel_mirna_predicted.csv` file and filter by various thresholds

### 1. miRDeep2 score of 10 or more


```r
colnames(novelmirna)
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
dim(novelmirna)	
```

```
## [1] 1201   17
```

```r
score10novelmirna<-novelmirna[novelmirna$miRDeep2.score >= 10, ]

dim(score10novelmirna)
```

```
## [1] 369  17
```

```r
nrow(score10novelmirna)
```

```
## [1] 369
```

So, there are 369 predicted miRNA precursors with a miRDeep2 score > or = 10

Estimated false positives is 34 +/- 6 (obtained from `1_extracted_mirdeep2_core_output_stats.csv`)


```r
34 - 6
```

```
## [1] 28
```

```r
34 + 6
```

```
## [1] 40
```

```r
28/369
```

```
## [1] 0.07588076
```

```r
40/369
```

```
## [1] 0.1084011
```

### 2. Significant Randfold p-value

Now, subset this again into those that had a significant Randfold p value, indicating ability of secondary structure formation


```r
head(score10novelmirna$significant.randfold.p.value)
```

```
## [1] yes yes yes yes yes yes
## Levels: no yes
```

```r
sum(score10novelmirna$significant.randfold.p.value=="yes")
```

```
## [1] 283
```

```r
randfoldsigpval<-score10novelmirna[score10novelmirna$significant.randfold.p.value == "yes", ]
dim(randfoldsigpval)
```

```
## [1] 283  17
```

Also checking that the Rfam alert is null for these sequences, meaning there are no other species of small RNA identified in these putative novel miRNA sequences


```r
rfam<-randfoldsigpval$rfam.alert
sum(rfam=="-")
```

```
## [1] 283
```

```r
sum(rfam !="-")
```

```
## [1] 0
```

### 3. Do the putative novel sequences have a homologous human miRNA seed sequence?


```r
sum(randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-")
```

```
## [1] 114
```

This indicates that 114 of the sequences have a homologous human seed sequence


```r
homologseed<-randfoldsigpval[randfoldsigpval$example.miRBase.miRNA.with.the.same.seed != "-",]
```

Subset of the 100 sequences


```r
dim(homologseed)
```

```
## [1] 114  17
```

```r
write.table(sts, file="../1_core_output_stats.csv", sep = "\t", col.names = TRUE)
write.table(novelmirna, file="../2_predicted_novel_mirna.csv", sep = "\t", col.names=TRUE)
write.table(md, file = "../3_mature_mirna_detected.csv", sep = "\t", col.names=TRUE)
write.table(homologseed, "../4_predicted_novel_sigranfold_homologseed.txt", sep = "\t", col.names=TRUE)
```

