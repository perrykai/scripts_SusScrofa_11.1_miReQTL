**Script:** `2_mirna_eqtl_summary.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`

**Date:**  12/6/17

**Input File Directory:**
1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`

2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`

**Input File(s):**

1. `1_gblup_results_summary.Rdata`

2. `2_gwa_results.Rdata`

3. `3_msuprp_mirna_gpdata.Rdata`

4. `4_normalized_dge_voom.Rdata`

5. `5_Z_G_miRNA.Rdata`

6. `6_mirna_precursor_annot_ssc11.csv`

**Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`

**Output File(s):** `2_miRNA_eqtl_summary.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Save data](#save-data)

## Objectives

The objective of this script is to summarize the number of eQTL peaks per miRNA output from the first eQTL analysis of the 174 F2 MSUPRP pig miRNA expression profiles.

## Install libraries



```r
rm(list=ls())
```

Load eqtl function Rdata containing stb function, which will summarize the eQTL results


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
ls()
```

```
##  [1] "absmap"     "add_legend" "AddPosGene" "distance"   "inrange"   
##  [6] "manhpt"     "peakrng"    "plot.GMA"   "sigpval"    "stb"       
## [11] "stb.nm"     "tbpos"      "zstandard"
```

```r
library(limma)
library(edgeR)
library(gwaR)
library(plyr)
library(corrplot)
```

```
## corrplot 0.84 loaded
```

```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts/")
```

## Load data

Load the gblup and eQTL output:


```r
load("../1_gblup_results_summary.Rdata")

load("../2_gwa_results.Rdata")

# #' Load the dge object to obtain the mature miRNA annotation:
# load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
```

Load the MSUPRP_miRNA gpdata object for the mapping information:


```r
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
```

Load the updated annotation file I created for my miR-eQTL miRNAs:


```r
annotation<-read.csv("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.csv", header=TRUE, na.strings="NA", colClasses=c("character","character","character","numeric","numeric","numeric","numeric","character","character","character"))
annotation<-annotation[,1:10]
save(annotation, file="../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
```

Load the Z matrix:


```r
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()
```

```
##  [1] "absmap"               "add_legend"           "AddPosGene"          
##  [4] "annotation"           "distance"             "G"                   
##  [7] "gblup.h2.se"          "inrange"              "manhpt"              
## [10] "MSUPRP_miRNA"         "peakrng"              "plot.GMA"            
## [13] "rst.gwa"              "sigpval"              "stb"                 
## [16] "stb.nm"               "summary_MSUPRP_miRNA" "summary.rst.gblup"   
## [19] "tbpos"                "Z"                    "zstandard"
```

## Analysis

### Summarize the heritability of the miRNAs, output from GBLUP:

The average heritability of all the miRNA expression profiles:


```r
mean(summary.rst.gblup$h2)
```

```
## [1] 0.119636
```

```r
summary(summary.rst.gblup$h2)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0000000 0.0000657 0.0852200 0.1196000 0.1852000 0.6309000
```

The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:

How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)


```r
sum(summary.rst.gblup$qvalue<0.05)
```

```
## [1] 46
```

Extract the significantly heritable miRNAs from the summary.rst.gblup dataset and calculate mean h2:


```r
sigh2<-summary.rst.gblup[summary.rst.gblup$qvalue<0.05,]
dim(sigh2)
```

```
## [1] 46 10
```

```r
mean(sigh2$h2)
```

```
## [1] 0.3350006
```

```r
summary(sigh2$h2)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1748  0.2553  0.3263  0.3350  0.3690  0.6309
```

Define the minimum p-value that is not significant (based on q-value < 0.05) as the threshold for plotting significant points


```r
summary(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.008129 0.082840 0.236200 0.267200 0.500000 0.500000
```

```r
sigthresh<-min(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh
```

```
## [1] 0.008128603
```

Plot h2 vs -log10(p-value) like before to determine trend in significance and h2:


```r
plot(summary.rst.gblup$h2, -log10(summary.rst.gblup$lrtpvalue),
    xlab = expression("Heritability"~(h^{2})),
    ylab = "-log10(p-value)",
    main = "Significance vs Heritability")
points(summary.rst.gblup$h2[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       -log10(summary.rst.gblup$lrtpvalue)[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       pch=19,col="red")
abline(a = -log10(sigthresh), b = 0, lty = 5)
```

![plot of chunk vol_sig_h2](figure/vol_sig_h2-1.tiff)

---

### The Summary Table of GWAS Results:

Assign the correct names to the different objects:


```r
map <- MSUPRP_miRNA$map
colnames(map)
```

```
## [1] "chr" "pos" "snp"
```

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
colnames(annotation)
```

```
##  [1] "miRNA"        "Ensembl.Gene" "Transcript"   "chr0"        
##  [5] "start"        "end"          "width"        "strand"      
##  [9] "type"         "Alias"
```

Annotation must contain columns "chr", "start", "end", and "strand"


```r
colnames(annotation)<-c("Name","gene", "transc", "chr","start","end","width","strand","type","Alias")
annotation$mid.mir<-round(rowMeans(annotation[,c("start", "end")], na.rm=TRUE))
head(annotation)
```

```
##              Name               gene               transc chr    start
## 1   ssc-let-7d-5p ENSSSCG00000019364 ENSSSCT00000020959.2   3 43468501
## 2      ssc-let-7g ENSSSCG00000019306 ENSSSCT00000020901.2  13 34406135
## 3     ssc-miR-128 ENSSSCG00000019503 ENSSSCT00000021098.2  15 16273557
## 4     ssc-miR-128 ENSSSCG00000019835 ENSSSCT00000021430.2  13 21038442
## 5 ssc-miR-1306-3p ENSSSCG00000019500 ENSSSCT00000021095.3  14 51473918
## 6  ssc-miR-140-5p ENSSSCG00000028710 ENSSSCT00000030513.2   6 17077517
##        end width strand                     type     Alias  mid.mir
## 1 43468596    96      + miRNA_primary_transcript MI0022120 43468548
## 2 34406222    88      - miRNA_primary_transcript MI0013087 34406178
## 3 16273638    82      - miRNA_primary_transcript MI0002451 16273598
## 4 21038525    84      + miRNA_primary_transcript MI0013094 21038484
## 5 51473997    80      + miRNA_primary_transcript MI0013148 51473958
## 6 17077610    94      - miRNA_primary_transcript MI0002437 17077564
```

Add the mid.mir to the map object for later use in manhattan plots (arrow where midpoint of miRNA precursor lies)

Need to remove the 2nd map position of miR-128 for this analysis:



```r
annotation<-annotation[-4,]
```

Extract the chromosome and the position of the miRNA precursor and build a data.frame to add to the map object later:


```r
mirpos<-data.frame(chr=annotation$chr,
	pos=annotation$mid.mir, row.names=annotation$Name)
str(mirpos)
```

```
## 'data.frame':	17 obs. of  2 variables:
##  $ chr: num  3 13 15 14 6 19 7 4 7 6 ...
##  $ pos: num  43468548 34406178 16273598 51473958 17077564 ...
```

```r
head(mirpos)
```

```
##                 chr      pos
## ssc-let-7d-5p     3 43468548
## ssc-let-7g       13 34406178
## ssc-miR-128      15 16273598
## ssc-miR-1306-3p  14 51473958
## ssc-miR-140-5p    6 17077564
## ssc-miR-1468     19 50337129
```

```r
dim(mirpos)
```

```
## [1] 17  2
```

```r
if(sum(mirpos$chr != annotation$chr, na.rm=TRUE) !=0) stop ("chr of miR did not add correctly")
if (sum(mirpos$pos != annotation$mid.mir, na.rm=TRUE) !=0) stop("mid-position of miRNA did not add correctly")

eqtlsum<-function(gwarst, map, annot, threshold=0.05, pergene=TRUE){
sigeqtl<-gwarst[gwarst$gwa.qval<threshold,]
mir<-sigeqtl$miRNA
head(mir)
snp<-as.character(sigeqtl$SNPid)
head(snp)
rownames(annot)<-annot$Name

eqtlrst<-data.frame(
	miRNA=mir, 
	chr.miR=annot[mir, "chr"],
	start.miR=annot[mir,"start"],
	end.miR=annot[mir,"end"],
	mid.miR=annot[mir,"mid.mir"],
	strand=as.character(annot[mir,"strand"]),
	miRBase.ID=as.character(annot[mir,"Alias"]),
	SNP=snp,
	chr.snp=map[match(snp, rownames(map)),"chr"],
	pos.snp=map[match(snp, rownames(map)),"pos"],
	snp.sign=as.character(sigeqtl$SNP.sign),
    pvalue=sigeqtl$gwa.pval,
	qvalue=sigeqtl$gwa.qval, 
    stringsAsFactors=FALSE
	)
if(pergene){
        id<-unique(eqtlrst[,"miRNA"])
        x<-list()

        for (i in id){
            a<-eqtlrst[eqtlrst[,"miRNA"]==i,]

            if (length(unique(a[,"chr.snp"]))==1){
                a<-a[order(a[,"pvalue"])[1],]

                } 
            else {

                b<-a[order(a[,"pvalue"]),]
                a<-list()

                    for (j in unique(b$chr.snp)){
                        a[[j]]<-b[b[,"chr.snp"]==j,][1,]
                }

            a<-do.call(rbind,a)
                }

        x[[i]]<-a

        }

    eqtlrst<-do.call(rbind,x)
    }
    rownames(eqtlrst)<-NULL
    return(eqtlrst)
}
```

Create the summary table, using pergene = TRUE to get gene-wise eQTL peaks:


```r
sum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=TRUE)
dim(sum.eqtl)
```

```
## [1] 22 13
```

```r
head(sum.eqtl)
```

```
##             miRNA chr.miR start.miR  end.miR  mid.miR strand miRBase.ID
## 1   ssc-let-7d-5p       3  43468501 43468596 43468548      +  MI0022120
## 2      ssc-let-7g      13  34406135 34406222 34406178      -  MI0013087
## 3     ssc-miR-128      15  16273557 16273638 16273598      -  MI0002451
## 4 ssc-miR-1306-3p      14  51473918 51473997 51473958      +  MI0013148
## 5  ssc-miR-140-5p       6  17077517 17077610 17077564      -  MI0002437
## 6  ssc-miR-140-5p       6  17077517 17077610 17077564      -  MI0002437
##           SNP chr.snp   pos.snp snp.sign       pvalue       qvalue
## 1 MARC0093624      15 122218534        - 3.093005e-07 0.0112251325
## 2 MARC0093624      15 122218534        - 1.643296e-07 0.0059638506
## 3 ALGA0023517       4  15332045        + 1.204780e-06 0.0437238842
## 4 H3GA0034702      12  52402985        + 1.689468e-07 0.0061314188
## 5 ALGA0117081       6  16973006        - 1.254506e-08 0.0004302663
## 6 ASGA0017748       4   7188213        - 1.146296e-06 0.0138671279
```

How many unique miRNAs have eQTL?


```r
length(unique(sum.eqtl$miRNA))
```

```
## [1] 17
```

```r
unique(as.character(sum.eqtl$miRNA))
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-140-5p"  "ssc-miR-1468"   
##  [7] "ssc-miR-184"     "ssc-miR-190b"    "ssc-miR-345-3p" 
## [10] "ssc-miR-429"     "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [13] "ssc-miR-874"     "ssc-miR-95"      "ssc-miR-9785-5p"
## [16] "ssc-miR-9810-3p" "ssc-miR-9843-3p"
```

miR-eQTL peaks at FDR<0.05:


```r
sum.eqtl
```

```
##              miRNA chr.miR start.miR   end.miR   mid.miR strand miRBase.ID
## 1    ssc-let-7d-5p       3  43468501  43468596  43468548      +  MI0022120
## 2       ssc-let-7g      13  34406135  34406222  34406178      -  MI0013087
## 3      ssc-miR-128      15  16273557  16273638  16273598      -  MI0002451
## 4  ssc-miR-1306-3p      14  51473918  51473997  51473958      +  MI0013148
## 5   ssc-miR-140-5p       6  17077517  17077610  17077564      -  MI0002437
## 6   ssc-miR-140-5p       6  17077517  17077610  17077564      -  MI0002437
## 7     ssc-miR-1468      19  50337088  50337170  50337129      -  MI0022160
## 8      ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 9      ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 10     ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 11    ssc-miR-190b       4  95540606  95540688  95540647      +  MI0017988
## 12  ssc-miR-345-3p       7 121193716 121193799 121193758      +  MI0013117
## 13     ssc-miR-429       6  63491921  63492001  63491961      +  MI0017991
## 14     ssc-miR-429       6  63491921  63492001  63491961      +  MI0017991
## 15 ssc-miR-6782-3p      NA        NA        NA       NaN   <NA>  MI0031620
## 16 ssc-miR-7135-3p       3  28371952  28372010  28371981      -  MI0023568
## 17     ssc-miR-874       2 139660118 139660201 139660160      -  MI0022157
## 18     ssc-miR-874       2 139660118 139660201 139660160      -  MI0022157
## 19      ssc-miR-95       8   3030934   3031014   3030974      +  MI0002436
## 20 ssc-miR-9785-5p      NA        NA        NA       NaN   <NA>  MI0031545
## 21 ssc-miR-9810-3p       4  83070363  83070457  83070410      +  MI0031577
## 22 ssc-miR-9843-3p       8 114110660 114110740 114110700      -  MI0031612
##            SNP chr.snp   pos.snp snp.sign       pvalue       qvalue
## 1  MARC0093624      15 122218534        - 3.093005e-07 1.122513e-02
## 2  MARC0093624      15 122218534        - 1.643296e-07 5.963851e-03
## 3  ALGA0023517       4  15332045        + 1.204780e-06 4.372388e-02
## 4  H3GA0034702      12  52402985        + 1.689468e-07 6.131419e-03
## 5  ALGA0117081       6  16973006        - 1.254506e-08 4.302663e-04
## 6  ASGA0017748       4   7188213        - 1.146296e-06 1.386713e-02
## 7  MARC0093624      15 122218534        - 5.550876e-07 2.014524e-02
## 8  ASGA0034057       7  50959933        - 3.850651e-11 1.996397e-07
## 9  M1GA0026172       6 169927307        + 1.294751e-05 1.423912e-02
## 10 DBWU0000430       3   9463123        - 4.560414e-05 4.137664e-02
## 11 ALGA0026452       4  87026192        - 4.686046e-07 1.700660e-02
## 12 H3GA0052416      15 121806256        - 1.705422e-06 3.581458e-02
## 13 ALGA0118516       6  63743210        + 2.207942e-09 5.359585e-06
## 14 ALGA0046283       8   7171284        + 7.776323e-05 3.101300e-02
## 15 DIAS0000707      10  24871779        - 8.537125e-09 1.549147e-04
## 16 ALGA0124095       3  28388183        - 2.905756e-06 1.172551e-02
## 17 ALGA0016550       2 139741258        + 7.919865e-14 2.874277e-09
## 18 ALGA0122273       3  61017470        - 5.676777e-05 2.064148e-02
## 19 MARC0093624      15 122218534        - 2.936732e-06 4.636241e-02
## 20 ALGA0121561       3   7321370        + 2.792947e-06 3.693384e-02
## 21 ALGA0030853       5  16390314        + 1.801621e-06 3.269221e-02
## 22 MARC0093624      15 122218534        + 8.887961e-08 3.225619e-03
```

```r
sum.eqtl<-sum.eqtl[order(sum.eqtl$miRNA, sum.eqtl$chr.snp),]
```

#### Summary of GWAS results at FDR < 0.05

Number of eQTL peaks per chromosome:


```r
table(sum.eqtl$chr.snp)
```

```
## 
## 10 12 15  2  3  4  5  6  7  8 
##  1  1  6  1  4  3  1  3  1  1
```

Names of associated miRNAs:


```r
table(sum.eqtl$miRNA)
```

```
## 
##   ssc-let-7d-5p      ssc-let-7g     ssc-miR-128 ssc-miR-1306-3p 
##               1               1               1               1 
##  ssc-miR-140-5p    ssc-miR-1468     ssc-miR-184    ssc-miR-190b 
##               2               1               3               1 
##  ssc-miR-345-3p     ssc-miR-429 ssc-miR-6782-3p ssc-miR-7135-3p 
##               1               2               1               1 
##     ssc-miR-874      ssc-miR-95 ssc-miR-9785-5p ssc-miR-9810-3p 
##               2               1               1               1 
## ssc-miR-9843-3p 
##               1
```

Chromosomes of associated miRNAs:


```r
table(sum.eqtl$chr.miR)
```

```
## 
##  2  3  4  6  7  8 13 14 15 19 
##  2  2  2  4  4  2  1  1  1  1
```

Names of associated peak markers:


```r
table(as.character(sum.eqtl$SNP))
```

```
## 
## ALGA0016550 ALGA0023517 ALGA0026452 ALGA0030853 ALGA0046283 ALGA0117081 
##           1           1           1           1           1           1 
## ALGA0118516 ALGA0121561 ALGA0122273 ALGA0124095 ASGA0017748 ASGA0034057 
##           1           1           1           1           1           1 
## DBWU0000430 DIAS0000707 H3GA0034702 H3GA0052416 M1GA0026172 MARC0093624 
##           1           1           1           1           1           5
```

---

### Determining the ranges of associated SNPs per eQTL peak on SSC15 (for ISAG abstract):

First, create the summary table at FDR 5% again, this time with pergene=F to identify all markers associated with each eQTL peak:


```r
fullsum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=FALSE)
dim(fullsum.eqtl)
```

```
## [1] 315  13
```

Summarize the number of SNPs associated with each miRNA eQTL (some have multiple peaks)


```r
numsnps<-by(fullsum.eqtl, as.character(fullsum.eqtl$miRNA), nrow)
numsnps<-ldply(numsnps, fun=NULL, id=names(numsnps))
colnames(numsnps) <- c("miRNA", "numsnps")
numsnps
```

```
##              miRNA numsnps
## 1    ssc-let-7d-5p       2
## 2       ssc-let-7g       2
## 3      ssc-miR-128       1
## 4  ssc-miR-1306-3p       1
## 5   ssc-miR-140-5p       3
## 6     ssc-miR-1468       1
## 7      ssc-miR-184      49
## 8     ssc-miR-190b       4
## 9   ssc-miR-345-3p       2
## 10     ssc-miR-429      91
## 11 ssc-miR-6782-3p       4
## 12 ssc-miR-7135-3p      14
## 13     ssc-miR-874     116
## 14      ssc-miR-95       3
## 15 ssc-miR-9785-5p      17
## 16 ssc-miR-9810-3p       2
## 17 ssc-miR-9843-3p       3
```

```r
sum(numsnps$numsnps)
```

```
## [1] 315
```

--- 

Examine the correlation of the SNPs within a peak with the same q-value


```r
id<-unique(as.character(fullsum.eqtl$miRNA))
pkmirqtl<-list()
corsnp<-list()
for(i in id){
    pkmirqtl[[i]]<-fullsum.eqtl[fullsum.eqtl$miRNA==i,]
    pkmirqtl[[i]]<-pkmirqtl[[i]][which(pkmirqtl[[i]]$qvalue == min(pkmirqtl[[i]]$qvalue)),]

if(length(as.character(unique(pkmirqtl[[i]]$chr.snp)))>1){
    pkmirqtl[[i]]<-pkmirqtl[[i]][pkmirqtl[[i]]$chr.snp==min(pkmirqtl[[i]]$chr.snp),]
}
if(sum(pkmirqtl[[i]]$qvalue == min(pkmirqtl[[i]]$qvalue))>1){
    corsnp[[i]]<-cor(Z[,as.character(pkmirqtl[[i]]$SNP)])
}
}

corsnp
```

```
## $`ssc-miR-140-5p`
##             ASGA0104343 ALGA0117081
## ASGA0104343   1.0000000   0.9932592
## ALGA0117081   0.9932592   1.0000000
## 
## $`ssc-miR-184`
##             ALGA0041952 M1GA0010398 ALGA0041972 ALGA0041993 H3GA0021767
## ALGA0041952           1           1           1           1           1
## M1GA0010398           1           1           1           1           1
## ALGA0041972           1           1           1           1           1
## ALGA0041993           1           1           1           1           1
## H3GA0021767           1           1           1           1           1
## ASGA0034057          -1          -1          -1          -1          -1
## DIAS0000025           1           1           1           1           1
##             ASGA0034057 DIAS0000025
## ALGA0041952          -1           1
## M1GA0010398          -1           1
## ALGA0041972          -1           1
## ALGA0041993          -1           1
## H3GA0021767          -1           1
## ASGA0034057           1          -1
## DIAS0000025          -1           1
## 
## $`ssc-miR-345-3p`
##             MARC0027291 H3GA0052416
## MARC0027291   1.0000000  -0.9879271
## H3GA0052416  -0.9879271   1.0000000
## 
## $`ssc-miR-429`
##             ASGA0094554 M1GA0008562 ASGA0090044 MARC0001121 MARC0082369
## ASGA0094554   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## M1GA0008562   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## ASGA0090044   0.9393866   0.9393866   1.0000000   1.0000000  -0.9393866
## MARC0001121   0.9393866   0.9393866   1.0000000   1.0000000  -0.9393866
## MARC0082369  -1.0000000  -1.0000000  -0.9393866  -0.9393866   1.0000000
## ALGA0118145   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## MARC0018157  -1.0000000  -1.0000000  -0.9393866  -0.9393866   1.0000000
## M1GA0024787  -1.0000000  -1.0000000  -0.9393866  -0.9393866   1.0000000
## ALGA0106537  -1.0000000  -1.0000000  -0.9393866  -0.9393866   1.0000000
## M1GA0009140   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## ASGA0030401   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## MARC0020138   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## MARC0011716   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## H3GA0056560  -1.0000000  -1.0000000  -0.9393866  -0.9393866   1.0000000
## M1GA0025298   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## MARC0048569   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## M1GA0026195   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## ASGA0095121   1.0000000   1.0000000   0.9393866   0.9393866  -1.0000000
## ALGA0118516  -0.9825149  -0.9825149  -0.9561155  -0.9561155   0.9825149
## ALGA0106326   0.9825149   0.9825149   0.9561155   0.9561155  -0.9825149
## ASGA0082593   0.9393866   0.9393866   1.0000000   1.0000000  -0.9393866
## MARC0030882   0.8981327   0.8981327   0.9540839   0.9540839  -0.8981327
## H3GA0054032   0.9401248   0.9401248   0.9082843   0.9082843  -0.9401248
##             ALGA0118145 MARC0018157 M1GA0024787 ALGA0106537 M1GA0009140
## ASGA0094554   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## M1GA0008562   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ASGA0090044   0.9393866  -0.9393866  -0.9393866  -0.9393866   0.9393866
## MARC0001121   0.9393866  -0.9393866  -0.9393866  -0.9393866   0.9393866
## MARC0082369  -1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0118145   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## MARC0018157  -1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## M1GA0024787  -1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0106537  -1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## M1GA0009140   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ASGA0030401   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## MARC0020138   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## MARC0011716   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## H3GA0056560  -1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## M1GA0025298   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## MARC0048569   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## M1GA0026195   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ASGA0095121   1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ALGA0118516  -0.9825149   0.9825149   0.9825149   0.9825149  -0.9825149
## ALGA0106326   0.9825149  -0.9825149  -0.9825149  -0.9825149   0.9825149
## ASGA0082593   0.9393866  -0.9393866  -0.9393866  -0.9393866   0.9393866
## MARC0030882   0.8981327  -0.8981327  -0.8981327  -0.8981327   0.8981327
## H3GA0054032   0.9401248  -0.9401248  -0.9401248  -0.9401248   0.9401248
##             ASGA0030401 MARC0020138 MARC0011716 H3GA0056560 M1GA0025298
## ASGA0094554   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## M1GA0008562   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## ASGA0090044   0.9393866   0.9393866   0.9393866  -0.9393866   0.9393866
## MARC0001121   0.9393866   0.9393866   0.9393866  -0.9393866   0.9393866
## MARC0082369  -1.0000000  -1.0000000  -1.0000000   1.0000000  -1.0000000
## ALGA0118145   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## MARC0018157  -1.0000000  -1.0000000  -1.0000000   1.0000000  -1.0000000
## M1GA0024787  -1.0000000  -1.0000000  -1.0000000   1.0000000  -1.0000000
## ALGA0106537  -1.0000000  -1.0000000  -1.0000000   1.0000000  -1.0000000
## M1GA0009140   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## ASGA0030401   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## MARC0020138   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## MARC0011716   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## H3GA0056560  -1.0000000  -1.0000000  -1.0000000   1.0000000  -1.0000000
## M1GA0025298   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## MARC0048569   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## M1GA0026195   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## ASGA0095121   1.0000000   1.0000000   1.0000000  -1.0000000   1.0000000
## ALGA0118516  -0.9825149  -0.9825149  -0.9825149   0.9825149  -0.9825149
## ALGA0106326   0.9825149   0.9825149   0.9825149  -0.9825149   0.9825149
## ASGA0082593   0.9393866   0.9393866   0.9393866  -0.9393866   0.9393866
## MARC0030882   0.8981327   0.8981327   0.8981327  -0.8981327   0.8981327
## H3GA0054032   0.9401248   0.9401248   0.9401248  -0.9401248   0.9401248
##             MARC0048569 M1GA0026195 ASGA0095121 ALGA0118516 ALGA0106326
## ASGA0094554   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## M1GA0008562   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## ASGA0090044   0.9393866   0.9393866   0.9393866  -0.9561155   0.9561155
## MARC0001121   0.9393866   0.9393866   0.9393866  -0.9561155   0.9561155
## MARC0082369  -1.0000000  -1.0000000  -1.0000000   0.9825149  -0.9825149
## ALGA0118145   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## MARC0018157  -1.0000000  -1.0000000  -1.0000000   0.9825149  -0.9825149
## M1GA0024787  -1.0000000  -1.0000000  -1.0000000   0.9825149  -0.9825149
## ALGA0106537  -1.0000000  -1.0000000  -1.0000000   0.9825149  -0.9825149
## M1GA0009140   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## ASGA0030401   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## MARC0020138   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## MARC0011716   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## H3GA0056560  -1.0000000  -1.0000000  -1.0000000   0.9825149  -0.9825149
## M1GA0025298   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## MARC0048569   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## M1GA0026195   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## ASGA0095121   1.0000000   1.0000000   1.0000000  -0.9825149   0.9825149
## ALGA0118516  -0.9825149  -0.9825149  -0.9825149   1.0000000  -1.0000000
## ALGA0106326   0.9825149   0.9825149   0.9825149  -1.0000000   1.0000000
## ASGA0082593   0.9393866   0.9393866   0.9393866  -0.9561155   0.9561155
## MARC0030882   0.8981327   0.8981327   0.8981327  -0.9138605   0.9138605
## H3GA0054032   0.9401248   0.9401248   0.9401248  -0.9565964   0.9565964
##             ASGA0082593 MARC0030882 H3GA0054032
## ASGA0094554   0.9393866   0.8981327   0.9401248
## M1GA0008562   0.9393866   0.8981327   0.9401248
## ASGA0090044   1.0000000   0.9540839   0.9082843
## MARC0001121   1.0000000   0.9540839   0.9082843
## MARC0082369  -0.9393866  -0.8981327  -0.9401248
## ALGA0118145   0.9393866   0.8981327   0.9401248
## MARC0018157  -0.9393866  -0.8981327  -0.9401248
## M1GA0024787  -0.9393866  -0.8981327  -0.9401248
## ALGA0106537  -0.9393866  -0.8981327  -0.9401248
## M1GA0009140   0.9393866   0.8981327   0.9401248
## ASGA0030401   0.9393866   0.8981327   0.9401248
## MARC0020138   0.9393866   0.8981327   0.9401248
## MARC0011716   0.9393866   0.8981327   0.9401248
## H3GA0056560  -0.9393866  -0.8981327  -0.9401248
## M1GA0025298   0.9393866   0.8981327   0.9401248
## MARC0048569   0.9393866   0.8981327   0.9401248
## M1GA0026195   0.9393866   0.8981327   0.9401248
## ASGA0095121   0.9393866   0.8981327   0.9401248
## ALGA0118516  -0.9561155  -0.9138605  -0.9565964
## ALGA0106326   0.9561155   0.9138605   0.9565964
## ASGA0082593   1.0000000   0.9540839   0.9082843
## MARC0030882   0.9540839   1.0000000   0.9553051
## H3GA0054032   0.9082843   0.9553051   1.0000000
## 
## $`ssc-miR-6782-3p`
##             ASGA0094215 DIAS0000707
## ASGA0094215           1          -1
## DIAS0000707          -1           1
## 
## $`ssc-miR-7135-3p`
##             MARC0056802 ALGA0018202 ASGA0014022 MARC0053267 ALGA0018219
## MARC0056802   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0018202   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ASGA0014022   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## MARC0053267   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0018219  -1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ALGA0112651   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0124243   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ASGA0097769   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ASGA0090944  -1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ALGA0018214  -1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ALGA0106253   1.0000000   1.0000000   1.0000000   1.0000000  -1.0000000
## ALGA0106252  -1.0000000  -1.0000000  -1.0000000  -1.0000000   1.0000000
## ALGA0124095   0.9856814   0.9856814   0.9856814   0.9856814  -0.9856814
## MARC0014155   0.9856814   0.9856814   0.9856814   0.9856814  -0.9856814
##             ALGA0112651 ALGA0124243 ASGA0097769 ASGA0090944 ALGA0018214
## MARC0056802   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ALGA0018202   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0014022   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## MARC0053267   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ALGA0018219  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0112651   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ALGA0124243   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0097769   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ASGA0090944  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0018214  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0106253   1.0000000   1.0000000   1.0000000  -1.0000000  -1.0000000
## ALGA0106252  -1.0000000  -1.0000000  -1.0000000   1.0000000   1.0000000
## ALGA0124095   0.9856814   0.9856814   0.9856814  -0.9856814  -0.9856814
## MARC0014155   0.9856814   0.9856814   0.9856814  -0.9856814  -0.9856814
##             ALGA0106253 ALGA0106252 ALGA0124095 MARC0014155
## MARC0056802   1.0000000  -1.0000000   0.9856814   0.9856814
## ALGA0018202   1.0000000  -1.0000000   0.9856814   0.9856814
## ASGA0014022   1.0000000  -1.0000000   0.9856814   0.9856814
## MARC0053267   1.0000000  -1.0000000   0.9856814   0.9856814
## ALGA0018219  -1.0000000   1.0000000  -0.9856814  -0.9856814
## ALGA0112651   1.0000000  -1.0000000   0.9856814   0.9856814
## ALGA0124243   1.0000000  -1.0000000   0.9856814   0.9856814
## ASGA0097769   1.0000000  -1.0000000   0.9856814   0.9856814
## ASGA0090944  -1.0000000   1.0000000  -0.9856814  -0.9856814
## ALGA0018214  -1.0000000   1.0000000  -0.9856814  -0.9856814
## ALGA0106253   1.0000000  -1.0000000   0.9856814   0.9856814
## ALGA0106252  -1.0000000   1.0000000  -0.9856814  -0.9856814
## ALGA0124095   0.9856814  -0.9856814   1.0000000   1.0000000
## MARC0014155   0.9856814  -0.9856814   1.0000000   1.0000000
## 
## $`ssc-miR-95`
##             MARC0027291 H3GA0052416 MARC0093624
## MARC0027291   1.0000000  -0.9879271   0.9018662
## H3GA0052416  -0.9879271   1.0000000  -0.8908862
## MARC0093624   0.9018662  -0.8908862   1.0000000
## 
## $`ssc-miR-9785-5p`
##             DRGA0003812 INRA0010427 ASGA0013843 ASGA0105049 ALGA0123606
## DRGA0003812   1.0000000  -1.0000000  -0.9151098   0.8515510   0.9216741
## INRA0010427  -1.0000000   1.0000000   0.9151098  -0.8515510  -0.9216741
## ASGA0013843  -0.9151098   0.9151098   1.0000000  -0.9424565  -0.8429895
## ASGA0105049   0.8515510  -0.8515510  -0.9424565   1.0000000   0.8965883
## ALGA0123606   0.9216741  -0.9216741  -0.8429895   0.8965883   1.0000000
## ALGA0123533  -0.8892990   0.8892990   0.8098592  -0.8629764  -0.9650816
## ASGA0099825  -0.4907968   0.4907968   0.5697135  -0.5244363  -0.4581610
## MARC0075820  -0.4311154   0.4311154   0.5047306  -0.4620722  -0.3683800
## ALGA0018363  -0.4210509   0.4210509   0.4935644  -0.4515275  -0.3930430
## ALGA0118443  -0.4976692   0.4976692   0.5978137  -0.5714427  -0.5001194
## MARC0063080   0.4976692  -0.4976692  -0.5978137   0.5714427   0.5001194
## MARC0058300   0.4976692  -0.4976692  -0.5978137   0.5714427   0.5001194
## ALGA0115191   0.4976692  -0.4976692  -0.5978137   0.5714427   0.5001194
## ASGA0014219   0.4382444  -0.4382444  -0.4868207   0.4646780   0.4276454
## ALGA0121561  -0.4054547   0.4054547   0.4742919  -0.4340794  -0.3788719
##             ALGA0123533 ASGA0099825 MARC0075820 ALGA0018363 ALGA0118443
## DRGA0003812  -0.8892990  -0.4907968  -0.4311154  -0.4210509  -0.4976692
## INRA0010427   0.8892990   0.4907968   0.4311154   0.4210509   0.4976692
## ASGA0013843   0.8098592   0.5697135   0.5047306   0.4935644   0.5978137
## ASGA0105049  -0.8629764  -0.5244363  -0.4620722  -0.4515275  -0.5714427
## ALGA0123606  -0.9650816  -0.4581610  -0.3683800  -0.3930430  -0.5001194
## ALGA0123533   1.0000000   0.4681922   0.3294758   0.3528962   0.4896591
## ASGA0099825   0.4681922   1.0000000   0.8374353   0.8321136   0.5798022
## MARC0075820   0.3294758   0.8374353   1.0000000   0.9622286   0.5328559
## ALGA0018363   0.3528962   0.8321136   0.9622286   1.0000000   0.5812150
## ALGA0118443   0.4896591   0.5798022   0.5328559   0.5812150   1.0000000
## MARC0063080  -0.4896591  -0.5798022  -0.5328559  -0.5812150  -1.0000000
## MARC0058300  -0.4896591  -0.5798022  -0.5328559  -0.5812150  -1.0000000
## ALGA0115191  -0.4896591  -0.5798022  -0.5328559  -0.5812150  -1.0000000
## ASGA0014219  -0.4322306  -0.8149094  -0.8908169  -0.9277458  -0.5922309
## ALGA0121561   0.3379881   0.8101413   0.9335425   0.9708320   0.5412495
##             MARC0063080 MARC0058300 ALGA0115191 ASGA0014219 ALGA0121561
## DRGA0003812   0.4976692   0.4976692   0.4976692   0.4382444  -0.4054547
## INRA0010427  -0.4976692  -0.4976692  -0.4976692  -0.4382444   0.4054547
## ASGA0013843  -0.5978137  -0.5978137  -0.5978137  -0.4868207   0.4742919
## ASGA0105049   0.5714427   0.5714427   0.5714427   0.4646780  -0.4340794
## ALGA0123606   0.5001194   0.5001194   0.5001194   0.4276454  -0.3788719
## ALGA0123533  -0.4896591  -0.4896591  -0.4896591  -0.4322306   0.3379881
## ASGA0099825  -0.5798022  -0.5798022  -0.5798022  -0.8149094   0.8101413
## MARC0075820  -0.5328559  -0.5328559  -0.5328559  -0.8908169   0.9335425
## ALGA0018363  -0.5812150  -0.5812150  -0.5812150  -0.9277458   0.9708320
## ALGA0118443  -1.0000000  -1.0000000  -1.0000000  -0.5922309   0.5412495
## MARC0063080   1.0000000   1.0000000   1.0000000   0.5922309  -0.5412495
## MARC0058300   1.0000000   1.0000000   1.0000000   0.5922309  -0.5412495
## ALGA0115191   1.0000000   1.0000000   1.0000000   0.5922309  -0.5412495
## ASGA0014219   0.5922309   0.5922309   0.5922309   1.0000000  -0.8974402
## ALGA0121561  -0.5412495  -0.5412495  -0.5412495  -0.8974402   1.0000000
## 
## $`ssc-miR-9810-3p`
##             MARC0021620 ALGA0030853
## MARC0021620           1           1
## ALGA0030853           1           1
```

---

### Extract peak range data from all miRNA eQTL peaks

I can obtain information on the range of each peak based on miRNA (adapted from Deborah's function "peakrng"):


```r
peakrngmir<-function(nmt, sumtb) {

# Positions of Snps in eQTL peak
    nmt <- nmt
    map.CI <- sumtb[sumtb$miRNA == nmt,c("SNP","chr.snp","pos.snp","qvalue")]

# Number of associated markers by chromosome
    chr <- table(as.character(map.CI$chr.snp))
    cat("Number of markers per chromosomal peak for",nmt,":")
    print(chr)
    idx <- as.numeric(names(chr))

    min.pos <- unlist(lapply(idx, function(x) min(map.CI[map.CI$chr.snp == x,"pos.snp"])))
    max.pos <- unlist(lapply(idx, function(x) max(map.CI[map.CI$chr.snp == x,"pos.snp"])))

    start.miR <- rep(sumtb[sumtb$miRNA == nmt,"start.miR"][1], length(idx))
    end.miR <- rep(sumtb[sumtb$miRNA == nmt,"end.miR"][1], length(idx))

# Identify the position of marker extreams for each peak
    peaks <- data.frame(miRNA=rep(nmt,length(idx)),
                chr.miR=rep(sumtb[sumtb$miRNA == nmt,"chr.miR"][1], length(idx)),
                start.miR=start.miR, end.miR=end.miR, range.miR=end.miR-start.miR,
                miRBase.ID=rep(sumtb[sumtb$miRNA == nmt,"miRBase.ID"][1], length(idx)),
                chr.snp=idx, range.peak=max.pos-min.pos,
                min.pos=min.pos, max.pos=max.pos, num.snp=as.vector(chr))

    return(peaks)
}
```

nmt = name transcript (in this case, miRNA); so, make a list of the miRNAs and loop through to get the peak range information for each miRNA with significant eQTL peaks
sumtb = output from the summary table function, with pergene = FALSE


```r
sigmirnames <- unique(as.character(sum.eqtl$miRNA))
sigmirnames
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-140-5p"  "ssc-miR-1468"   
##  [7] "ssc-miR-184"     "ssc-miR-190b"    "ssc-miR-345-3p" 
## [10] "ssc-miR-429"     "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [13] "ssc-miR-874"     "ssc-miR-95"      "ssc-miR-9785-5p"
## [16] "ssc-miR-9810-3p" "ssc-miR-9843-3p"
```

```r
mirpeaks<-data.frame(do.call(rbind, lapply(sigmirnames, peakrngmir, fullsum.eqtl)), sum.eqtl[,c("SNP","pos.snp","qvalue")])
```

```
## Number of markers per chromosomal peak for ssc-let-7d-5p :
## 15 
##  2 
## Number of markers per chromosomal peak for ssc-let-7g :
## 15 
##  2 
## Number of markers per chromosomal peak for ssc-miR-128 :
## 4 
## 1 
## Number of markers per chromosomal peak for ssc-miR-1306-3p :
## 12 
##  1 
## Number of markers per chromosomal peak for ssc-miR-140-5p :
## 4 6 
## 1 2 
## Number of markers per chromosomal peak for ssc-miR-1468 :
## 15 
##  1 
## Number of markers per chromosomal peak for ssc-miR-184 :
##  3  6  7 
##  2  1 46 
## Number of markers per chromosomal peak for ssc-miR-190b :
## 4 
## 4 
## Number of markers per chromosomal peak for ssc-miR-345-3p :
## 15 
##  2 
## Number of markers per chromosomal peak for ssc-miR-429 :
##  6  8 
## 90  1 
## Number of markers per chromosomal peak for ssc-miR-6782-3p :
## 10 
##  4 
## Number of markers per chromosomal peak for ssc-miR-7135-3p :
##  3 
## 14 
## Number of markers per chromosomal peak for ssc-miR-874 :
##   2   3 
## 115   1 
## Number of markers per chromosomal peak for ssc-miR-95 :
## 15 
##  3 
## Number of markers per chromosomal peak for ssc-miR-9785-5p :
##  3 
## 17 
## Number of markers per chromosomal peak for ssc-miR-9810-3p :
## 5 
## 2 
## Number of markers per chromosomal peak for ssc-miR-9843-3p :
## 15 
##  3
```

```r
mirpeaks
```

```
##              miRNA chr.miR start.miR   end.miR range.miR miRBase.ID
## 1    ssc-let-7d-5p       3  43468501  43468596        95  MI0022120
## 2       ssc-let-7g      13  34406135  34406222        87  MI0013087
## 3      ssc-miR-128      15  16273557  16273638        81  MI0002451
## 4  ssc-miR-1306-3p      14  51473918  51473997        79  MI0013148
## 5   ssc-miR-140-5p       6  17077517  17077610        93  MI0002437
## 6   ssc-miR-140-5p       6  17077517  17077610        93  MI0002437
## 7     ssc-miR-1468      19  50337088  50337170        82  MI0022160
## 8      ssc-miR-184       7  48345017  48345099        82  MI0002421
## 9      ssc-miR-184       7  48345017  48345099        82  MI0002421
## 10     ssc-miR-184       7  48345017  48345099        82  MI0002421
## 11    ssc-miR-190b       4  95540606  95540688        82  MI0017988
## 12  ssc-miR-345-3p       7 121193716 121193799        83  MI0013117
## 13     ssc-miR-429       6  63491921  63492001        80  MI0017991
## 14     ssc-miR-429       6  63491921  63492001        80  MI0017991
## 15 ssc-miR-6782-3p      NA        NA        NA        NA  MI0031620
## 16 ssc-miR-7135-3p       3  28371952  28372010        58  MI0023568
## 17     ssc-miR-874       2 139660118 139660201        83  MI0022157
## 18     ssc-miR-874       2 139660118 139660201        83  MI0022157
## 19      ssc-miR-95       8   3030934   3031014        80  MI0002436
## 20 ssc-miR-9785-5p      NA        NA        NA        NA  MI0031545
## 21 ssc-miR-9810-3p       4  83070363  83070457        94  MI0031577
## 22 ssc-miR-9843-3p       8 114110660 114110740        80  MI0031612
##    chr.snp range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 1       15     345721 121872813 122218534       2 MARC0093624 122218534
## 2       15     345721 121872813 122218534       2 MARC0093624 122218534
## 3        4          0  15332045  15332045       1 ALGA0023517  15332045
## 4       12          0  52402985  52402985       1 H3GA0034702  52402985
## 5        4          0   7188213   7188213       1 ASGA0017748   7188213
## 6        6     389705  16583301  16973006       2 ALGA0117081  16973006
## 7       15          0 122218534 122218534       1 MARC0093624 122218534
## 8        3  117042625   9463123 126505748       2 DBWU0000430   9463123
## 9        6          0 169927307 169927307       1 M1GA0026172 169927307
## 10       7   43229663  43560865  86790528      46 ASGA0034057  50959933
## 11       4   14432476  87026192 101458668       4 ALGA0026452  87026192
## 12      15      66557 121806256 121872813       2 H3GA0052416 121806256
## 13       6   23594459  42970844  66565303      90 ALGA0118516  63743210
## 14       8          0   7171284   7171284       1 ALGA0046283   7171284
## 15      10    9750125  18022870  27772995       4 DIAS0000707  24871779
## 16       3     704850  27716152  28421002      14 ALGA0124095  28388183
## 17       2   21713538 127788285 149501823     115 ALGA0016550 139741258
## 18       3          0  61017470  61017470       1 ALGA0122273  61017470
## 19      15     412278 121806256 122218534       3 MARC0093624 122218534
## 20       3   27074955   7321370  34396325      17 ALGA0121561   7321370
## 21       5      76120  16314194  16390314       2 ALGA0030853  16390314
## 22      15     412278 121806256 122218534       3 MARC0093624 122218534
##          qvalue
## 1  1.122513e-02
## 2  5.963851e-03
## 3  4.372388e-02
## 4  6.131419e-03
## 5  1.386713e-02
## 6  4.302663e-04
## 7  2.014524e-02
## 8  4.137664e-02
## 9  1.423912e-02
## 10 1.996397e-07
## 11 1.700660e-02
## 12 3.581458e-02
## 13 5.359585e-06
## 14 3.101300e-02
## 15 1.549147e-04
## 16 1.172551e-02
## 17 2.874277e-09
## 18 2.064148e-02
## 19 4.636241e-02
## 20 3.693384e-02
## 21 3.269221e-02
## 22 3.225619e-03
```

Separate the two miR-184 SNPs on SSC3:


```r
fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$chr.snp == "3",]
```

```
##          miRNA chr.miR start.miR  end.miR  mid.miR strand miRBase.ID
## 11 ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
## 12 ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
##            SNP chr.snp   pos.snp snp.sign       pvalue     qvalue
## 11 DBWU0000430       3   9463123        - 4.560414e-05 0.04137664
## 12 ASGA0016793       3 126505748        + 4.813730e-05 0.04159521
```

```r
map["DBWU0000430",]
```

```
##             chr     pos   snp
## DBWU0000430   3 9463123 [A/G]
```

```r
rownames(map[map$pos=="9463123",])
```

```
## [1] "DBWU0000430"
```

```r
mir184.1<-mirpeaks[8,]
mir184.1$range.peak<-0
mir184.1$max.pos<-mir184.1$min.pos
mir184.1$num.snp<-1
mir184.1
```

```
##         miRNA chr.miR start.miR  end.miR range.miR miRBase.ID chr.snp
## 8 ssc-miR-184       7  48345017 48345099        82  MI0002421       3
##   range.peak min.pos max.pos num.snp         SNP pos.snp     qvalue
## 8          0 9463123 9463123       1 DBWU0000430 9463123 0.04137664
```

```r
map["ASGA0016793",]
```

```
##             chr       pos   snp
## ASGA0016793   3 126505748 [A/G]
```

```r
rownames(map[map$pos=="126505748",])
```

```
## [1] "ASGA0016793"
```

```r
mir184.2<-mirpeaks[8,]
mir184.2$range.peak<-0
mir184.2$min.pos<-mir184.2$max.pos
mir184.2$num.snp<-1
mir184.2$SNP<-rownames(map[map$pos=="126505748",])
mir184.2$pos.snp<-map["ASGA0016793","pos"]
mir184.2$qvalue<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP == "ASGA0016793","qvalue"]
mir184.2
```

```
##         miRNA chr.miR start.miR  end.miR range.miR miRBase.ID chr.snp
## 8 ssc-miR-184       7  48345017 48345099        82  MI0002421       3
##   range.peak   min.pos   max.pos num.snp         SNP   pos.snp     qvalue
## 8          0 126505748 126505748       1 ASGA0016793 126505748 0.04159521
```

```r
mirpeaks[8,]<-mir184.1

mirpeaks[(nrow(mirpeaks)+1),]<-mir184.2

mirpeaks<-mirpeaks[order(mirpeaks$miRNA, mirpeaks$chr.snp),]
rownames(mirpeaks)<- NULL

mirpeaks
```

```
##              miRNA chr.miR start.miR   end.miR range.miR miRBase.ID
## 1    ssc-let-7d-5p       3  43468501  43468596        95  MI0022120
## 2       ssc-let-7g      13  34406135  34406222        87  MI0013087
## 3      ssc-miR-128      15  16273557  16273638        81  MI0002451
## 4  ssc-miR-1306-3p      14  51473918  51473997        79  MI0013148
## 5   ssc-miR-140-5p       6  17077517  17077610        93  MI0002437
## 6   ssc-miR-140-5p       6  17077517  17077610        93  MI0002437
## 7     ssc-miR-1468      19  50337088  50337170        82  MI0022160
## 8      ssc-miR-184       7  48345017  48345099        82  MI0002421
## 9      ssc-miR-184       7  48345017  48345099        82  MI0002421
## 10     ssc-miR-184       7  48345017  48345099        82  MI0002421
## 11     ssc-miR-184       7  48345017  48345099        82  MI0002421
## 12    ssc-miR-190b       4  95540606  95540688        82  MI0017988
## 13  ssc-miR-345-3p       7 121193716 121193799        83  MI0013117
## 14     ssc-miR-429       6  63491921  63492001        80  MI0017991
## 15     ssc-miR-429       6  63491921  63492001        80  MI0017991
## 16 ssc-miR-6782-3p      NA        NA        NA        NA  MI0031620
## 17 ssc-miR-7135-3p       3  28371952  28372010        58  MI0023568
## 18     ssc-miR-874       2 139660118 139660201        83  MI0022157
## 19     ssc-miR-874       2 139660118 139660201        83  MI0022157
## 20      ssc-miR-95       8   3030934   3031014        80  MI0002436
## 21 ssc-miR-9785-5p      NA        NA        NA        NA  MI0031545
## 22 ssc-miR-9810-3p       4  83070363  83070457        94  MI0031577
## 23 ssc-miR-9843-3p       8 114110660 114110740        80  MI0031612
##    chr.snp range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 1       15     345721 121872813 122218534       2 MARC0093624 122218534
## 2       15     345721 121872813 122218534       2 MARC0093624 122218534
## 3        4          0  15332045  15332045       1 ALGA0023517  15332045
## 4       12          0  52402985  52402985       1 H3GA0034702  52402985
## 5        4          0   7188213   7188213       1 ASGA0017748   7188213
## 6        6     389705  16583301  16973006       2 ALGA0117081  16973006
## 7       15          0 122218534 122218534       1 MARC0093624 122218534
## 8        3          0   9463123   9463123       1 DBWU0000430   9463123
## 9        3          0 126505748 126505748       1 ASGA0016793 126505748
## 10       6          0 169927307 169927307       1 M1GA0026172 169927307
## 11       7   43229663  43560865  86790528      46 ASGA0034057  50959933
## 12       4   14432476  87026192 101458668       4 ALGA0026452  87026192
## 13      15      66557 121806256 121872813       2 H3GA0052416 121806256
## 14       6   23594459  42970844  66565303      90 ALGA0118516  63743210
## 15       8          0   7171284   7171284       1 ALGA0046283   7171284
## 16      10    9750125  18022870  27772995       4 DIAS0000707  24871779
## 17       3     704850  27716152  28421002      14 ALGA0124095  28388183
## 18       2   21713538 127788285 149501823     115 ALGA0016550 139741258
## 19       3          0  61017470  61017470       1 ALGA0122273  61017470
## 20      15     412278 121806256 122218534       3 MARC0093624 122218534
## 21       3   27074955   7321370  34396325      17 ALGA0121561   7321370
## 22       5      76120  16314194  16390314       2 ALGA0030853  16390314
## 23      15     412278 121806256 122218534       3 MARC0093624 122218534
##          qvalue
## 1  1.122513e-02
## 2  5.963851e-03
## 3  4.372388e-02
## 4  6.131419e-03
## 5  1.386713e-02
## 6  4.302663e-04
## 7  2.014524e-02
## 8  4.137664e-02
## 9  4.159521e-02
## 10 1.423912e-02
## 11 1.996397e-07
## 12 1.700660e-02
## 13 3.581458e-02
## 14 5.359585e-06
## 15 3.101300e-02
## 16 1.549147e-04
## 17 1.172551e-02
## 18 2.874277e-09
## 19 2.064148e-02
## 20 4.636241e-02
## 21 3.693384e-02
## 22 3.269221e-02
## 23 3.225619e-03
```

Add this change to sum.eqtl object:


```r
mir184.1
```

```
##         miRNA chr.miR start.miR  end.miR range.miR miRBase.ID chr.snp
## 8 ssc-miR-184       7  48345017 48345099        82  MI0002421       3
##   range.peak min.pos max.pos num.snp         SNP pos.snp     qvalue
## 8          0 9463123 9463123       1 DBWU0000430 9463123 0.04137664
```

```r
mir184.2
```

```
##         miRNA chr.miR start.miR  end.miR range.miR miRBase.ID chr.snp
## 8 ssc-miR-184       7  48345017 48345099        82  MI0002421       3
##   range.peak   min.pos   max.pos num.snp         SNP   pos.snp     qvalue
## 8          0 126505748 126505748       1 ASGA0016793 126505748 0.04159521
```

```r
tmp<-sum.eqtl[8,1:7]
tmp<-cbind(tmp, 
    SNP=as.character(mir184.2$SNP),
    chr.snp=mir184.2$chr.snp, 
    pos.snp=mir184.2$pos.snp, 
    snp.sign=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "snp.sign"], 
    pvalue=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "pvalue"],
    qvalue=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "qvalue"])

tmp
```

```
##          miRNA chr.miR start.miR  end.miR  mid.miR strand miRBase.ID
## 10 ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
##            SNP chr.snp   pos.snp snp.sign      pvalue     qvalue
## 10 ASGA0016793       3 126505748        + 4.81373e-05 0.04159521
```

```r
sum.eqtl<-rbind(sum.eqtl, tmp)
sum.eqtl<-sum.eqtl[order(sum.eqtl$miRNA, sum.eqtl$chr.snp),]
rownames(sum.eqtl)<-NULL
sum.eqtl
```

```
##              miRNA chr.miR start.miR   end.miR   mid.miR strand miRBase.ID
## 1    ssc-let-7d-5p       3  43468501  43468596  43468548      +  MI0022120
## 2       ssc-let-7g      13  34406135  34406222  34406178      -  MI0013087
## 3      ssc-miR-128      15  16273557  16273638  16273598      -  MI0002451
## 4  ssc-miR-1306-3p      14  51473918  51473997  51473958      +  MI0013148
## 5   ssc-miR-140-5p       6  17077517  17077610  17077564      -  MI0002437
## 6   ssc-miR-140-5p       6  17077517  17077610  17077564      -  MI0002437
## 7     ssc-miR-1468      19  50337088  50337170  50337129      -  MI0022160
## 8      ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 9      ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 10     ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 11     ssc-miR-184       7  48345017  48345099  48345058      +  MI0002421
## 12    ssc-miR-190b       4  95540606  95540688  95540647      +  MI0017988
## 13  ssc-miR-345-3p       7 121193716 121193799 121193758      +  MI0013117
## 14     ssc-miR-429       6  63491921  63492001  63491961      +  MI0017991
## 15     ssc-miR-429       6  63491921  63492001  63491961      +  MI0017991
## 16 ssc-miR-6782-3p      NA        NA        NA       NaN   <NA>  MI0031620
## 17 ssc-miR-7135-3p       3  28371952  28372010  28371981      -  MI0023568
## 18     ssc-miR-874       2 139660118 139660201 139660160      -  MI0022157
## 19     ssc-miR-874       2 139660118 139660201 139660160      -  MI0022157
## 20      ssc-miR-95       8   3030934   3031014   3030974      +  MI0002436
## 21 ssc-miR-9785-5p      NA        NA        NA       NaN   <NA>  MI0031545
## 22 ssc-miR-9810-3p       4  83070363  83070457  83070410      +  MI0031577
## 23 ssc-miR-9843-3p       8 114110660 114110740 114110700      -  MI0031612
##            SNP chr.snp   pos.snp snp.sign       pvalue       qvalue
## 1  MARC0093624      15 122218534        - 3.093005e-07 1.122513e-02
## 2  MARC0093624      15 122218534        - 1.643296e-07 5.963851e-03
## 3  ALGA0023517       4  15332045        + 1.204780e-06 4.372388e-02
## 4  H3GA0034702      12  52402985        + 1.689468e-07 6.131419e-03
## 5  ASGA0017748       4   7188213        - 1.146296e-06 1.386713e-02
## 6  ALGA0117081       6  16973006        - 1.254506e-08 4.302663e-04
## 7  MARC0093624      15 122218534        - 5.550876e-07 2.014524e-02
## 8  DBWU0000430       3   9463123        - 4.560414e-05 4.137664e-02
## 9  ASGA0016793       3 126505748        + 4.813730e-05 4.159521e-02
## 10 M1GA0026172       6 169927307        + 1.294751e-05 1.423912e-02
## 11 ASGA0034057       7  50959933        - 3.850651e-11 1.996397e-07
## 12 ALGA0026452       4  87026192        - 4.686046e-07 1.700660e-02
## 13 H3GA0052416      15 121806256        - 1.705422e-06 3.581458e-02
## 14 ALGA0118516       6  63743210        + 2.207942e-09 5.359585e-06
## 15 ALGA0046283       8   7171284        + 7.776323e-05 3.101300e-02
## 16 DIAS0000707      10  24871779        - 8.537125e-09 1.549147e-04
## 17 ALGA0124095       3  28388183        - 2.905756e-06 1.172551e-02
## 18 ALGA0016550       2 139741258        + 7.919865e-14 2.874277e-09
## 19 ALGA0122273       3  61017470        - 5.676777e-05 2.064148e-02
## 20 MARC0093624      15 122218534        - 2.936732e-06 4.636241e-02
## 21 ALGA0121561       3   7321370        + 2.792947e-06 3.693384e-02
## 22 ALGA0030853       5  16390314        + 1.801621e-06 3.269221e-02
## 23 MARC0093624      15 122218534        + 8.887961e-08 3.225619e-03
```

```r
sum.eqtl[sum.eqtl$miRNA=="ssc-miR-184",]
```

```
##          miRNA chr.miR start.miR  end.miR  mid.miR strand miRBase.ID
## 8  ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
## 9  ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
## 10 ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
## 11 ssc-miR-184       7  48345017 48345099 48345058      +  MI0002421
##            SNP chr.snp   pos.snp snp.sign       pvalue       qvalue
## 8  DBWU0000430       3   9463123        - 4.560414e-05 4.137664e-02
## 9  ASGA0016793       3 126505748        + 4.813730e-05 4.159521e-02
## 10 M1GA0026172       6 169927307        + 1.294751e-05 1.423912e-02
## 11 ASGA0034057       7  50959933        - 3.850651e-11 1.996397e-07
```

```r
mirpeaks[mirpeaks$miRNA=="ssc-miR-184",]
```

```
##          miRNA chr.miR start.miR  end.miR range.miR miRBase.ID chr.snp
## 8  ssc-miR-184       7  48345017 48345099        82  MI0002421       3
## 9  ssc-miR-184       7  48345017 48345099        82  MI0002421       3
## 10 ssc-miR-184       7  48345017 48345099        82  MI0002421       6
## 11 ssc-miR-184       7  48345017 48345099        82  MI0002421       7
##    range.peak   min.pos   max.pos num.snp         SNP   pos.snp
## 8           0   9463123   9463123       1 DBWU0000430   9463123
## 9           0 126505748 126505748       1 ASGA0016793 126505748
## 10          0 169927307 169927307       1 M1GA0026172 169927307
## 11   43229663  43560865  86790528      46 ASGA0034057  50959933
##          qvalue
## 8  4.137664e-02
## 9  4.159521e-02
## 10 1.423912e-02
## 11 1.996397e-07
```

```r
rm(mir184.1)
rm(mir184.2)
rm(tmp)
```

---

### Creating Manhattan plots of the six miRNA with the highest numbers of associated SNP markers (for ISAG poster)

First, convert the map positions to absolute map positions using Deborah's function "absmap"

Notice how the map object's positions are relative to chromosome:


```r
head(map)
```

```
##             chr    pos   snp
## MARC0044150   1 205163 [A/G]
## ASGA0000014   1 261794 [A/C]
## H3GA0000026   1 309120 [A/G]
## ASGA0000021   1 289363 [A/C]
## ALGA0000009   1 408640 [T/C]
## ALGA0000014   1 381075 [T/C]
```

```r
dim(map)
```

```
## [1] 36292     3
```

Add the map positions of the miRNA with eQTL


```r
map.full<-rbind(map[,1:2], mirpos[sigmirnames,])
dim(map.full)
```

```
## [1] 36309     2
```

```r
head(map.full)
```

```
##             chr    pos
## MARC0044150   1 205163
## ASGA0000014   1 261794
## H3GA0000026   1 309120
## ASGA0000021   1 289363
## ALGA0000009   1 408640
## ALGA0000014   1 381075
```

```r
if(nrow(map.full) - nrow(map) != length(sigmirnames)) stop ("miRNA map positions not added correctly")
```

Use the absmap function to convert the chromosomal map positions to absolute map positions:


```r
absposmap<-absmap(map.full)
head(absposmap)
```

```
## MARC0044150 ASGA0000005 ASGA0000014 ASGA0000021 H3GA0000026 H3GA0000032 
##      205163      239523      261794      289363      309120      373626
```

```r
tail(absposmap)
```

```
##  ALGA0098927  ALGA0098928  ASGA0080447  ASGA0080449  M1GA0023446 
##   2258421445   2258498427   2258598096   2258759259   2258777666 
## ssc-miR-1468 
##   2309114795
```

Notice it didn't include the miRNA with no map position (ssc-miR-140-5p)


```r
absposmap[which(names(absposmap) %in% sigmirnames)]
```

```
##     ssc-miR-874 ssc-miR-7135-3p   ssc-let-7d-5p ssc-miR-9810-3p 
##       413948022       454014796       469111363       641413652 
##    ssc-miR-190b  ssc-miR-140-5p     ssc-miR-429     ssc-miR-184 
##       653883889       810720748       857135145      1012764421 
##  ssc-miR-345-3p      ssc-miR-95 ssc-miR-9843-3p      ssc-let-7g 
##      1085613121      1089201828      1200281554      1607599735 
## ssc-miR-1306-3p     ssc-miR-128    ssc-miR-1468 
##      1832773780      1939265163      2309114795
```

Divide by 1e6 to get absolute position in Mb (nicer x-axis for plots)


```r
head(absposmap/1e6)
```

```
## MARC0044150 ASGA0000005 ASGA0000014 ASGA0000021 H3GA0000026 H3GA0000032 
##    0.205163    0.239523    0.261794    0.289363    0.309120    0.373626
```

Use sigpval function (from DV) to calculate the significant p-value cutoff for plotting manhattan plots for each miRNA

Then, provide the vector of miRNA names and the absposmap object to the manhpt function (also from Deborah's func_eqtl.Rdata) and loop through the vector of miRNA names to create the Manhattan plots:
## Visualize

Correlation plots for miReQTL SNPs with same qvalue within peak:



```r
for(i in 1:length(names(corsnp))){
    corrplot(corsnp[[i]], method="number",mar=c(2,2,2,2), cl.cex=0.7, tl.cex=0.8, tl.col="black", title=names(pkmirqtl[i]))
}
```

![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-1.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-2.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-3.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-4.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-5.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-6.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-7.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-8.tiff)![plot of chunk corrplot_mirSNP](figure/corrplot_mirSNP-9.tiff)

Look again at the plots with many SNPs with similar qvalues:


```r
for(i in names(pkmirqtl)[c(10, 12, 15)]){
    corrplot.mixed(corsnp[[i]], tl.pos="lt", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title=names(pkmirqtl[i]))
}
```

![plot of chunk corrplot_mirSNP2](figure/corrplot_mirSNP2-1.tiff)![plot of chunk corrplot_mirSNP2](figure/corrplot_mirSNP2-2.tiff)![plot of chunk corrplot_mirSNP2](figure/corrplot_mirSNP2-3.tiff)

---


```r
sigmirnames
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-140-5p"  "ssc-miR-1468"   
##  [7] "ssc-miR-184"     "ssc-miR-190b"    "ssc-miR-345-3p" 
## [10] "ssc-miR-429"     "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [13] "ssc-miR-874"     "ssc-miR-95"      "ssc-miR-9785-5p"
## [16] "ssc-miR-9810-3p" "ssc-miR-9843-3p"
```

```r
for(i in sigmirnames){
	pvalcutoff<-sigpval(i, pvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.pval"], qvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.qval"], fdr=0.05)
	pvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.pval"]
	names(pvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=pvals,map=map.full,annotation=NULL,pvalues=TRUE,cutoff=pvalcutoff, arrow=TRUE)

}
```

![plot of chunk man_plot_pval](figure/man_plot_pval-1.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-2.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-3.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-4.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-5.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-6.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-7.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-8.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-9.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-10.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-11.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-12.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-13.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-14.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-15.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-16.tiff)![plot of chunk man_plot_pval](figure/man_plot_pval-17.tiff)

The Manhattan Plots of q-values should be put on equal y-axes for comparison on poster.

Three of the peaks are strong enough signals to be on a y-axis of 0-10 (-log10(qval))

The remaining 14 peaks are less strong, on a y-axis of 0-4



```r
sigmirnames4<-sigmirnames[c(1:4,6,8,9,11,12,14:17)]
sigmirnames4
```

```
##  [1] "ssc-let-7d-5p"   "ssc-let-7g"      "ssc-miR-128"    
##  [4] "ssc-miR-1306-3p" "ssc-miR-1468"    "ssc-miR-190b"   
##  [7] "ssc-miR-345-3p"  "ssc-miR-6782-3p" "ssc-miR-7135-3p"
## [10] "ssc-miR-95"      "ssc-miR-9785-5p" "ssc-miR-9810-3p"
## [13] "ssc-miR-9843-3p"
```

```r
sigmirnames8<-sigmirnames[c(5,7,10)]
sigmirnames8
```

```
## [1] "ssc-miR-140-5p" "ssc-miR-184"    "ssc-miR-429"
```

```r
sigmirnames10<-sigmirnames[13]
sigmirnames10
```

```
## [1] "ssc-miR-874"
```

```r
for(i in sigmirnames4){
	qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
	names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,4))

}
```

![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-1.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-2.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-3.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-4.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-5.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-6.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-7.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-8.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-9.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-10.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-11.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-12.tiff)![plot of chunk man_plot_qval_y4](figure/man_plot_qval_y4-13.tiff)

```r
for(i in sigmirnames8){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,8))

}
```

![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-1.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-2.tiff)![plot of chunk man_plot_qval_y8](figure/man_plot_qval_y8-3.tiff)

```r
for(i in sigmirnames10){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,10))

}
```

![plot of chunk man_plot_qval_y10](figure/man_plot_qval_y10-1.tiff)

## Save data


```r
save(sum.eqtl, fullsum.eqtl, absposmap, map.full, mirpeaks, file = "../3_eqtl_summary_tables_maps.Rdata")
write.table(sum.eqtl, file="../3_eqtl_summary.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(fullsum.eqtl, file="../3_eqtl_fullsummary.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

