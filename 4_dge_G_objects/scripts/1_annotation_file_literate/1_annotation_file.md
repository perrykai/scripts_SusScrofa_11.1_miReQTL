**Script:** `0_create_annotation_file_for_dge_object.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/scripts`

**Date:**  12/4/17

**Input File Directory:**  
1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/reference_sequences`
2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/3_miRNA_expression_mx/`                

**Input File(s):** 
1. `ssc.gff3`
2. `1_mean_mature_mirna_exp.Rdata` 

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`

**Output File(s):** 

1. `1_precursor_mirna_annotation.Rdata`
2. `2_mature_mirna_annotation.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives
Create a data frame of gene annotation and assembly data with dimensions genes x categories, where genes here is the same as genes in count matrix
## Install libraries


```r
library(rtracklayer)
```

```
## Loading required package: methods
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: BiocGenerics
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
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

## Load data


```r
rm(list=ls())
```

Load expression matrix:


```r
load("../../3_miRNA_expression_mx/1_mean_mature_mirna_exp.Rdata")
```

Load annotation file:


```r
GFFfile <- "../../reference_sequences/ssc.gff3"
gff <- import.gff(GFFfile, version = "3")

slotNames(gff)
```

```
## [1] "seqnames"        "ranges"          "strand"          "elementMetadata"
## [5] "seqinfo"         "metadata"
```

```r
gff
```

```
## GRanges object with 761 ranges and 8 metadata columns:
##         seqnames                 ranges strand   |   source
##            <Rle>              <IRanges>  <Rle>   | <factor>
##     [1]     chr1     [1882209, 1882304]      -   |     <NA>
##     [2]     chr1     [1882219, 1882238]      -   |     <NA>
##     [3]     chr1     [2600232, 2600312]      -   |     <NA>
##     [4]     chr1     [2600282, 2600302]      -   |     <NA>
##     [5]     chr1     [2816327, 2816402]      -   |     <NA>
##     ...      ...                    ...    ... ...      ...
##   [757]     chrX [140861161, 140861182]      -   |     <NA>
##   [758]     chrX [141091501, 141091581]      -   |     <NA>
##   [759]     chrX [141091550, 141091569]      -   |     <NA>
##   [760]     chrX [141093619, 141093698]      -   |     <NA>
##   [761]     chrX [141093667, 141093686]      -   |     <NA>
##                             type     score     phase           ID
##                         <factor> <numeric> <integer>  <character>
##     [1] miRNA_primary_transcript      <NA>      <NA>    MI0031582
##     [2]                    miRNA      <NA>      <NA> MIMAT0037033
##     [3] miRNA_primary_transcript      <NA>      <NA>    MI0031588
##     [4]                    miRNA      <NA>      <NA> MIMAT0037039
##     [5] miRNA_primary_transcript      <NA>      <NA>    MI0031591
##     ...                      ...       ...       ...          ...
##   [757]                    miRNA      <NA>      <NA> MIMAT0025374
##   [758] miRNA_primary_transcript      <NA>      <NA>    MI0002410
##   [759]                    miRNA      <NA>      <NA> MIMAT0002116
##   [760] miRNA_primary_transcript      <NA>      <NA>    MI0002411
##   [761]                    miRNA      <NA>      <NA> MIMAT0002117
##                   Alias            Name Derives_from
##         <CharacterList>     <character>  <character>
##     [1]       MI0031582    ssc-mir-9815         <NA>
##     [2]    MIMAT0037033 ssc-miR-9815-3p    MI0031582
##     [3]       MI0031588    ssc-mir-9821         <NA>
##     [4]    MIMAT0037039 ssc-miR-9821-5p    MI0031588
##     [5]       MI0031591    ssc-mir-9824         <NA>
##     ...             ...             ...          ...
##   [757]    MIMAT0025374     ssc-miR-452    MI0022149
##   [758]       MI0002410   ssc-mir-105-1         <NA>
##   [759]    MIMAT0002116   ssc-miR-105-1    MI0002410
##   [760]       MI0002411   ssc-mir-105-2         <NA>
##   [761]    MIMAT0002117   ssc-miR-105-2    MI0002411
##   -------
##   seqinfo: 19 sequences from an unspecified genome; no seqlengths
```

## Analysis

Create a data frame of annotation data including the following columns:

Name, chr0, start, end, width, strand, type, Alias, Derives_from

Extract specific fields of information, build annotation data frame:


```r
mirannot <- data.frame(unlist(elementMetadata(gff)$Name),
				  as.data.frame(seqnames(gff)),
				  as.data.frame(ranges(gff)),
                  as.data.frame(strand(gff)),
                  unlist(elementMetadata(gff)$type),
                  unlist(elementMetadata(gff)$Alias),
                  unlist(elementMetadata(gff)$Derives_from)
                  )


colnames(mirannot) <- c("Name", 
					  "chr0",
					  "start", "end", "width", 
					  "strand", 
					  "type", 
					  "Alias", 
					  "Derives_from"
					  )
```

Check the class of each column in the data frame:


```r
colclass<-NULL

for (i in colnames(mirannot)) {
   colclass<-c(colclass,class(mirannot[,i]))
}

colclass
```

```
## [1] "factor"  "factor"  "integer" "integer" "integer" "factor"  "factor" 
## [8] "factor"  "factor"
```

```r
head(mirannot)
```

```
##              Name chr0   start     end width strand
## 1    ssc-mir-9815 chr1 1882209 1882304    96      -
## 2 ssc-miR-9815-3p chr1 1882219 1882238    20      -
## 3    ssc-mir-9821 chr1 2600232 2600312    81      -
## 4 ssc-miR-9821-5p chr1 2600282 2600302    21      -
## 5    ssc-mir-9824 chr1 2816327 2816402    76      -
## 6 ssc-miR-9824-5p chr1 2816372 2816392    21      -
##                       type        Alias Derives_from
## 1 miRNA_primary_transcript    MI0031582         <NA>
## 2                    miRNA MIMAT0037033    MI0031582
## 3 miRNA_primary_transcript    MI0031588         <NA>
## 4                    miRNA MIMAT0037039    MI0031588
## 5 miRNA_primary_transcript    MI0031591         <NA>
## 6                    miRNA MIMAT0037042    MI0031591
```

### Isolate precursor miRNA annotation/assembly information


```r
sum(mirannot$type=="miRNA_primary_transcript")
```

```
## [1] 341
```

```r
precursor.mirannot <- mirannot[(mirannot$type=="miRNA_primary_transcript"),]

dim(precursor.mirannot)
```

```
## [1] 341   9
```

```r
head(precursor.mirannot)
```

```
##            Name chr0    start      end width strand
## 1  ssc-mir-9815 chr1  1882209  1882304    96      -
## 3  ssc-mir-9821 chr1  2600232  2600312    81      -
## 5  ssc-mir-9824 chr1  2816327  2816402    76      -
## 7  ssc-mir-9831 chr1  4108414  4108490    77      +
## 9  ssc-mir-9857 chr1 12666392 12666461    70      -
## 11 ssc-mir-9787 chr1 12884779 12884869    91      -
##                        type     Alias Derives_from
## 1  miRNA_primary_transcript MI0031582         <NA>
## 3  miRNA_primary_transcript MI0031588         <NA>
## 5  miRNA_primary_transcript MI0031591         <NA>
## 7  miRNA_primary_transcript MI0031600         <NA>
## 9  miRNA_primary_transcript MI0031635         <NA>
## 11 miRNA_primary_transcript MI0031548         <NA>
```

```r
if (sum(precursor.mirannot$type == "miRNA") != 0) stop ("mature miRNA in primary miRNA annotation")
```

### Isolate mature miRNA annotation/assembly information


```r
sum(mirannot$type==("miRNA"))
```

```
## [1] 420
```

```r
mature.mirannot <- mirannot[(mirannot$type=="miRNA"),]
```

Order the rows by mature miRNA Name


```r
mature.mirannot <- mature.mirannot[order(mature.mirannot$Name),]

dim(mature.mirannot)
```

```
## [1] 420   9
```

```r
head(mature.mirannot)
```

```
##              Name  chr0     start       end width strand  type
## 434    ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## 664    ssc-let-7a  chr9  54378110  54378131    22      - miRNA
## 207    ssc-let-7c chr13 191559351 191559372    22      + miRNA
## 439 ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## 438 ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## 541    ssc-let-7e  chr6  51858385  51858406    22      + miRNA
##            Alias Derives_from
## 434 MIMAT0013865    MI0017984
## 664 MIMAT0013865    MI0013085
## 207 MIMAT0002151    MI0002445
## 439 MIMAT0025357    MI0022120
## 438 MIMAT0025356    MI0022120
## 541 MIMAT0013866    MI0013086
```

```r
tail(mature.mirannot)
```

```
##                Name  chr0     start       end width strand  type
## 698 ssc-miR-9859-3p  chrX  45362398  45362421    24      + miRNA
## 609 ssc-miR-9860-5p  chr7 106158401 106158421    21      - miRNA
## 558 ssc-miR-9861-5p  chr6  64986156  64986176    21      - miRNA
## 269 ssc-miR-9862-3p chr15 155563825 155563844    20      + miRNA
## 205     ssc-miR-99a chr13 191558621 191558642    22      + miRNA
## 539     ssc-miR-99b  chr6  51858223  51858244    22      + miRNA
##            Alias Derives_from
## 698 MIMAT0037078    MI0031637
## 609 MIMAT0037079    MI0031638
## 558 MIMAT0037080    MI0031639
## 269 MIMAT0037084    MI0031643
## 205 MIMAT0013896    MI0013114
## 539 MIMAT0006018    MI0007077
```

```r
if (sum(mature.mirannot$type == "miRNA_primary_transcript") != 0) stop ("primary miRNA in mature miRNA annotation")

sum(table(mature.mirannot$Name)>1)
```

```
## [1] 51
```

### Extract the duplicated mature miRNAs, order by miRNA name:


```r
dups<-data.frame(name=mature.mirannot$Name,dup1=duplicated(mature.mirannot$Name),dup2=duplicated(mature.mirannot$Name,fromLast=TRUE))
```

The "duplicated" function returns "TRUE" for the second and higher occurrence of a given element.

Executed "duplicated" on the data frame with normal order of rows, then reversed order of rows, then combined the logical index vectors.


```r
head(dups)
```

```
##            name  dup1  dup2
## 1    ssc-let-7a FALSE  TRUE
## 2    ssc-let-7a  TRUE FALSE
## 3    ssc-let-7c FALSE FALSE
## 4 ssc-let-7d-3p FALSE FALSE
## 5 ssc-let-7d-5p FALSE FALSE
## 6    ssc-let-7e FALSE FALSE
```

```r
rs<-rowSums(dups[,2:3])
```

As a check, calculated the rowSums of the data.frame and observed where the overlaps were.

Then, to calculate the total number of assembled duplicated miRNAs, subtract the overlaps from the total sum:


```r
rs
```

```
##   [1] 1 1 0 0 0 0 1 1 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1
##  [36] 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1
##  [71] 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [106] 0 0 0 0 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0
## [141] 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0
## [176] 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 0
## [211] 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [246] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0
## [281] 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 1 1
## [316] 0 0 0 0 0 1 1 0 1 1 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 2 2
## [351] 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 2 1 0 0 0 0 0 0 0 0 0
## [386] 0 0 1 1 0 0 0 0 1 2 2 2 2 1 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
```

```r
sum(rs)
```

```
## [1] 116
```

```r
sum(rs>1)
```

```
## [1] 7
```

```r
sum(rs)-(sum(rs>1))
```

```
## [1] 109
```

```r
dup.maturemir <- mature.mirannot[(duplicated(mature.mirannot$Name) | duplicated(mature.mirannot$Name, fromLast=TRUE)),]
dim(dup.maturemir)
```

```
## [1] 109   9
```

```r
head(dup.maturemir)
```

```
##            Name chr0     start       end width strand  type        Alias
## 434  ssc-let-7a chr3  44864443  44864464    22      + miRNA MIMAT0013865
## 664  ssc-let-7a chr9  54378110  54378131    22      - miRNA MIMAT0013865
## 436  ssc-let-7f chr3  44864810  44864831    22      + miRNA MIMAT0002152
## 709  ssc-let-7f chrX  51801782  51801803    22      - miRNA MIMAT0002152
## 46  ssc-miR-101 chr1 242988435 242988455    21      - miRNA MIMAT0010185
## 569 ssc-miR-101 chr6 135736164 135736184    21      + miRNA MIMAT0010185
##     Derives_from
## 434    MI0017984
## 664    MI0013085
## 436    MI0022121
## 709    MI0002446
## 46     MI0010678
## 569    MI0010679
```

```r
tail(dup.maturemir)
```

```
##                Name  chr0     start       end width strand  type
## 322 ssc-miR-9845-5p chr17  69657812  69657831    20      + miRNA
## 324 ssc-miR-9845-5p chr17  69661606  69661625    20      + miRNA
## 265 ssc-miR-9847-3p chr15 155366725 155366745    21      + miRNA
## 267 ssc-miR-9847-3p chr15 155444903 155444923    21      + miRNA
## 449 ssc-miR-9855-5p  chr3 103546153 103546173    21      + miRNA
## 451 ssc-miR-9855-5p  chr3 103707560 103707580    21      - miRNA
##            Alias Derives_from
## 322 MIMAT0037063    MI0031618
## 324 MIMAT0037063    MI0031619
## 265 MIMAT0037066    MI0031622
## 267 MIMAT0037066    MI0031623
## 449 MIMAT0037074    MI0031632
## 451 MIMAT0037074    MI0031633
```

### Extract the precursor IDs for duplicated mature miRNAs


```r
dup_derive <- as.data.frame(tapply(as.character(dup.maturemir$Derives_from), as.character(dup.maturemir$Name), paste))
dup_derive$Name <- rownames(dup_derive)
colnames(dup_derive) <- c("mult_precursors", "Name")

head(dup_derive)
```

```
##                   mult_precursors         Name
## ssc-let-7a   MI0017984, MI0013085   ssc-let-7a
## ssc-let-7f   MI0022121, MI0002446   ssc-let-7f
## ssc-miR-101  MI0010678, MI0010679  ssc-miR-101
## ssc-miR-103  MI0002448, MI0013105  ssc-miR-103
## ssc-miR-124a MI0002450, MI0010680 ssc-miR-124a
## ssc-miR-128  MI0013094, MI0002451  ssc-miR-128
```

```r
dim(dup_derive)
```

```
## [1] 51  2
```

```r
if ((nrow(dup_derive) == length(unique(dup.maturemir$Name))) != "TRUE") stop ("extraction of duplicated mature miRNAs did not work correctly")
```

### Build a data.frame containing the mature miRNAs, without duplicate information (all mature miRNAs present only once)



```r
sum(duplicated(mature.mirannot$Name)==FALSE)
```

```
## [1] 362
```

```r
single.mature.mirannot<-mature.mirannot[duplicated(mature.mirannot$Name)==FALSE,]

dim(single.mature.mirannot)
```

```
## [1] 362   9
```

```r
head(single.mature.mirannot)
```

```
##              Name  chr0     start       end width strand  type
## 434    ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## 207    ssc-let-7c chr13 191559351 191559372    22      + miRNA
## 439 ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## 438 ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## 541    ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## 436    ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##            Alias Derives_from
## 434 MIMAT0013865    MI0017984
## 207 MIMAT0002151    MI0002445
## 439 MIMAT0025357    MI0022120
## 438 MIMAT0025356    MI0022120
## 541 MIMAT0013866    MI0013086
## 436 MIMAT0002152    MI0022121
```

```r
if (sum(duplicated(single.mature.mirannot$Name))!=0) stop ("nonunique miRNA present in single.mature.mirannot")
```

### Merge the single.mature.mirannot and the dup_derive data frames to obtain the multiple precursor sequences for a duplicated mature miRNA:


```r
merged.single.mature.mirannot<-merge(single.mature.mirannot, dup_derive, by = "Name", all.x = TRUE)
rownames(merged.single.mature.mirannot)<-merged.single.mature.mirannot$Name

dim(merged.single.mature.mirannot)
```

```
## [1] 362  10
```

```r
head(merged.single.mature.mirannot)
```

```
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##                      Alias Derives_from      mult_precursors
## ssc-let-7a    MIMAT0013865    MI0017984 MI0017984, MI0013085
## ssc-let-7c    MIMAT0002151    MI0002445                   NA
## ssc-let-7d-3p MIMAT0025357    MI0022120                   NA
## ssc-let-7d-5p MIMAT0025356    MI0022120                   NA
## ssc-let-7e    MIMAT0013866    MI0013086                   NA
## ssc-let-7f    MIMAT0002152    MI0022121 MI0022121, MI0002446
```

```r
if (sum(merged.single.mature.mirannot$Name != single.mature.mirannot$Name) != 0) stop ("merged single miRNA names not the same as single miRNA names")
if (sum(rownames(merged.single.mature.mirannot) != merged.single.mature.mirannot$Name) != 0) stop ("merged.single.mature.mirannot rownames not the same as its miRNA column")
```

Determine which miRNAs lack assembly information and extract them:


```r
summary(match(rownames(no.zero.dfmeanrcround), merged.single.mature.mirannot$Name))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##    1.00   73.25  146.50  148.70  220.80  362.00      45
```

```r
length(rownames(no.zero.dfmeanrcround))
```

```
## [1] 335
```

```r
length(merged.single.mature.mirannot$Name)
```

```
## [1] 362
```

```r
noannot<-rownames(no.zero.dfmeanrcround)[is.na(match(rownames(no.zero.dfmeanrcround), merged.single.mature.mirannot$Name))]

length(noannot)
```

```
## [1] 45
```

The following 45 miRNAs lack assembly information:


```r
noannot
```

```
##  [1] "ssc-miR-1277"    "ssc-miR-1307"    "ssc-miR-140-3p" 
##  [4] "ssc-miR-140-5p"  "ssc-miR-151-3p"  "ssc-miR-151-5p" 
##  [7] "ssc-miR-153"     "ssc-miR-15a"     "ssc-miR-218"    
## [10] "ssc-miR-22-3p"   "ssc-miR-22-5p"   "ssc-miR-26a"    
## [13] "ssc-miR-299"     "ssc-miR-301"     "ssc-miR-323"    
## [16] "ssc-miR-342"     "ssc-miR-34a"     "ssc-miR-3613"   
## [19] "ssc-miR-369"     "ssc-miR-370"     "ssc-miR-376a-3p"
## [22] "ssc-miR-376a-5p" "ssc-miR-376b"    "ssc-miR-376c"   
## [25] "ssc-miR-381"     "ssc-miR-382"     "ssc-miR-411"    
## [28] "ssc-miR-424-3p"  "ssc-miR-424-5p"  "ssc-miR-4332"   
## [31] "ssc-miR-4338"    "ssc-miR-450a"    "ssc-miR-450b-3p"
## [34] "ssc-miR-450b-5p" "ssc-miR-450c-3p" "ssc-miR-450c-5p"
## [37] "ssc-miR-487b"    "ssc-miR-490"     "ssc-miR-503"    
## [40] "ssc-miR-542-3p"  "ssc-miR-542-5p"  "ssc-miR-652"    
## [43] "ssc-miR-758"     "ssc-miR-885-3p"  "ssc-miR-885-5p"
```

```r
noannot<-data.frame(Name=noannot, chr0=NA, start=NA, end=NA, width=NA, strand=NA, type=NA, Alias=NA, Derives_from=NA, mult_precursors=NA)
dim(noannot)
```

```
## [1] 45 10
```

```r
head(noannot)
```

```
##             Name chr0 start end width strand type Alias Derives_from
## 1   ssc-miR-1277   NA    NA  NA    NA     NA   NA    NA           NA
## 2   ssc-miR-1307   NA    NA  NA    NA     NA   NA    NA           NA
## 3 ssc-miR-140-3p   NA    NA  NA    NA     NA   NA    NA           NA
## 4 ssc-miR-140-5p   NA    NA  NA    NA     NA   NA    NA           NA
## 5 ssc-miR-151-3p   NA    NA  NA    NA     NA   NA    NA           NA
## 6 ssc-miR-151-5p   NA    NA  NA    NA     NA   NA    NA           NA
##   mult_precursors
## 1              NA
## 2              NA
## 3              NA
## 4              NA
## 5              NA
## 6              NA
```

Add the miRNAs lacking assembly information into the annotation data frame:

Extract the miRNA names from the count matrix:


```r
count.mirnas<-as.data.frame(rownames(no.zero.dfmeanrcround))
colnames(count.mirnas)<- "Name"
dim(count.mirnas)
```

```
## [1] 335   1
```

```r
head(count.mirnas)
```

```
##            Name
## 1    ssc-let-7a
## 2    ssc-let-7c
## 3 ssc-let-7d-3p
## 4 ssc-let-7d-5p
## 5    ssc-let-7e
## 6    ssc-let-7f
```

Merge the annotation file with the miRNA names from the count matrix:


```r
total.mature.annot<-merge(count.mirnas, merged.single.mature.mirannot, by = "Name", all.x = TRUE)
rownames(total.mature.annot)<-total.mature.annot$Name

dim(total.mature.annot)
```

```
## [1] 335  10
```

```r
head(total.mature.annot)
```

```
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##                      Alias Derives_from      mult_precursors
## ssc-let-7a    MIMAT0013865    MI0017984 MI0017984, MI0013085
## ssc-let-7c    MIMAT0002151    MI0002445                   NA
## ssc-let-7d-3p MIMAT0025357    MI0022120                   NA
## ssc-let-7d-5p MIMAT0025356    MI0022120                   NA
## ssc-let-7e    MIMAT0013866    MI0013086                   NA
## ssc-let-7f    MIMAT0002152    MI0022121 MI0022121, MI0002446
```

Check that the unassembled miRNAs bound correctly:


```r
match(noannot$Name, total.mature.annot$Name)
```

```
##  [1]  30  40  54  55  73  74  76  79 132 141 143 154 162 166 180 193 196
## [18] 198 205 206 211 212 213 214 217 218 220 225 226 233 236 237 238 239
## [35] 240 241 247 249 259 264 265 273 306 311 312
```

```r
total.mature.annot[match(noannot$Name, total.mature.annot$Name),]
```

```
##                            Name chr0 start end width strand type Alias
## ssc-miR-1277       ssc-miR-1277 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-1307       ssc-miR-1307 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-140-3p   ssc-miR-140-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-140-5p   ssc-miR-140-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-151-3p   ssc-miR-151-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-151-5p   ssc-miR-151-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-153         ssc-miR-153 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-15a         ssc-miR-15a <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-218         ssc-miR-218 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-22-3p     ssc-miR-22-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-22-5p     ssc-miR-22-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-26a         ssc-miR-26a <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-299         ssc-miR-299 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-301         ssc-miR-301 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-323         ssc-miR-323 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-342         ssc-miR-342 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-34a         ssc-miR-34a <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-3613       ssc-miR-3613 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-369         ssc-miR-369 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-370         ssc-miR-370 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-376a-3p ssc-miR-376a-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-376a-5p ssc-miR-376a-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-376b       ssc-miR-376b <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-376c       ssc-miR-376c <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-381         ssc-miR-381 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-382         ssc-miR-382 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-411         ssc-miR-411 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-424-3p   ssc-miR-424-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-424-5p   ssc-miR-424-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-4332       ssc-miR-4332 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-4338       ssc-miR-4338 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-450a       ssc-miR-450a <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-450b-3p ssc-miR-450b-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-450b-5p ssc-miR-450b-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-450c-3p ssc-miR-450c-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-450c-5p ssc-miR-450c-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-487b       ssc-miR-487b <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-490         ssc-miR-490 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-503         ssc-miR-503 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-542-3p   ssc-miR-542-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-542-5p   ssc-miR-542-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-652         ssc-miR-652 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-758         ssc-miR-758 <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-885-3p   ssc-miR-885-3p <NA>    NA  NA    NA   <NA> <NA>  <NA>
## ssc-miR-885-5p   ssc-miR-885-5p <NA>    NA  NA    NA   <NA> <NA>  <NA>
##                 Derives_from mult_precursors
## ssc-miR-1277            <NA>              NA
## ssc-miR-1307            <NA>              NA
## ssc-miR-140-3p          <NA>              NA
## ssc-miR-140-5p          <NA>              NA
## ssc-miR-151-3p          <NA>              NA
## ssc-miR-151-5p          <NA>              NA
## ssc-miR-153             <NA>              NA
## ssc-miR-15a             <NA>              NA
## ssc-miR-218             <NA>              NA
## ssc-miR-22-3p           <NA>              NA
## ssc-miR-22-5p           <NA>              NA
## ssc-miR-26a             <NA>              NA
## ssc-miR-299             <NA>              NA
## ssc-miR-301             <NA>              NA
## ssc-miR-323             <NA>              NA
## ssc-miR-342             <NA>              NA
## ssc-miR-34a             <NA>              NA
## ssc-miR-3613            <NA>              NA
## ssc-miR-369             <NA>              NA
## ssc-miR-370             <NA>              NA
## ssc-miR-376a-3p         <NA>              NA
## ssc-miR-376a-5p         <NA>              NA
## ssc-miR-376b            <NA>              NA
## ssc-miR-376c            <NA>              NA
## ssc-miR-381             <NA>              NA
## ssc-miR-382             <NA>              NA
## ssc-miR-411             <NA>              NA
## ssc-miR-424-3p          <NA>              NA
## ssc-miR-424-5p          <NA>              NA
## ssc-miR-4332            <NA>              NA
## ssc-miR-4338            <NA>              NA
## ssc-miR-450a            <NA>              NA
## ssc-miR-450b-3p         <NA>              NA
## ssc-miR-450b-5p         <NA>              NA
## ssc-miR-450c-3p         <NA>              NA
## ssc-miR-450c-5p         <NA>              NA
## ssc-miR-487b            <NA>              NA
## ssc-miR-490             <NA>              NA
## ssc-miR-503             <NA>              NA
## ssc-miR-542-3p          <NA>              NA
## ssc-miR-542-5p          <NA>              NA
## ssc-miR-652             <NA>              NA
## ssc-miR-758             <NA>              NA
## ssc-miR-885-3p          <NA>              NA
## ssc-miR-885-5p          <NA>              NA
```

```r
if (sum(rownames(no.zero.dfmeanrcround)!=rownames(total.mature.annot))!=0) stop ("rownames of count mx and annotation df not equal")
```

Combine the two columns of precursor information into one:


```r
multp<-sapply(total.mature.annot[,"mult_precursors"],paste,collapse=",")
singlep<-as.character(total.mature.annot[,"Derives_from"])
multp[multp=="NA"]<-singlep[multp=="NA"]
head(multp)
```

```
##            ssc-let-7a            ssc-let-7a            ssc-let-7a 
## "MI0017984,MI0013085"           "MI0002445"           "MI0022120" 
##            ssc-let-7a            ssc-let-7a            ssc-let-7f 
##           "MI0022120"           "MI0013086" "MI0022121,MI0002446"
```

```r
total.mature.annot$Precursors<-as.character(multp)

total.mature.annot2<-total.mature.annot[,-c(9,10)]

if (sum(rownames(total.mature.annot2)!=count.mirnas$Name)!=0) stop ("rownames of count mx and annotation df not equal")
if (sum(rownames(total.mature.annot2)!=rownames(no.zero.dfmeanrcround))!=0) stop ("rownames of count mx and annotation df not equal")

head(total.mature.annot2)
```

```
##                        Name  chr0     start       end width strand  type
## ssc-let-7a       ssc-let-7a  chr3  44864443  44864464    22      + miRNA
## ssc-let-7c       ssc-let-7c chr13 191559351 191559372    22      + miRNA
## ssc-let-7d-3p ssc-let-7d-3p  chr3  44867331  44867352    22      + miRNA
## ssc-let-7d-5p ssc-let-7d-5p  chr3  44867277  44867298    22      + miRNA
## ssc-let-7e       ssc-let-7e  chr6  51858385  51858406    22      + miRNA
## ssc-let-7f       ssc-let-7f  chr3  44864810  44864831    22      + miRNA
##                      Alias          Precursors
## ssc-let-7a    MIMAT0013865 MI0017984,MI0013085
## ssc-let-7c    MIMAT0002151           MI0002445
## ssc-let-7d-3p MIMAT0025357           MI0022120
## ssc-let-7d-5p MIMAT0025356           MI0022120
## ssc-let-7e    MIMAT0013866           MI0013086
## ssc-let-7f    MIMAT0002152 MI0022121,MI0002446
```

## Save data

Save both annotation of mature miRNAs for the eQTL analysis, and annotation of precursor sequences for reference. 


```r
save(precursor.mirannot, file = "../1_precursor_mirna_annotation.Rdata")
save(total.mature.annot2, file = "../2_mature_mirna_annotation.Rdata")
```

