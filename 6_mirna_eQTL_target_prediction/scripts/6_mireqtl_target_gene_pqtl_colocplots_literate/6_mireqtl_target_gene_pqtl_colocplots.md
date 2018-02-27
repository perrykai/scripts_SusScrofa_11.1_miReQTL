**Script:** `6_mireqtl_target_gene_pqtl_colocplots.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`

**Date:**  `12/18/17`

**Input File Directory:**  

1. `/mnt/research/pigeqtl/analyses/eQTL/paper/output`

1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/`

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/5_gblup_gwa_eqtl`

1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/6_mirna_eQTL_target_prediction`

**Input File(s):** 

1. `funct_eqtl.Rdata`

1. `pQTL_60K.Rdata`

1. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`

1. `10_mrna_mirna_corr_rst.Rdata` `13_target_mrna_coloc_pqtl.Rdata`

**Output File Directory:** 

`/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts/6_mireqtl_target_gene_pqtl_colocplots_literate/figure`

**Output File(s):** 

`ssc-miR-874-1.tiff`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

This scripts creates a graphical representation of the co-localization of ssc-miR-874 target genes with pQTL
by first plotting the manhattan plot for ssc-miR-874 eQTL and then overlappling the pQTL with co-localized
ssc-miR-874 target genes. The targets genes apping within each pQTL peak is highlighted in red.

## Install libraries

## Load data

Clear environment


```r
rm(list=ls())
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
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] knitr_1.17
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5    tools_3.2.0     stringi_1.1.1   methods_3.2.0  
## [5] stringr_1.2.0   evaluate_0.10.1
```

```r
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")
```

Load required R objects


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")

load("../../5_gblup_gwa_eqtl/2_gwa_results.Rdata")
load("../../5_gblup_gwa_eqtl/3_eqtl_summary_tables_maps.Rdata")
load("../10_mrna_mirna_corr_rst.Rdata")
load("../13_target_mrna_coloc_pqtl.Rdata")

ls()
```

```
##  [1] "absmap"              "absposmap"           "add_legend"         
##  [4] "AddPosGene"          "annot"               "coloc"              
##  [7] "distance"            "fullsum.eqtl"        "GBLUP"              
## [10] "GWAS"                "inrange"             "manhpt"             
## [13] "map.full"            "mirpeaks"            "negcoloc"           
## [16] "peakrng"             "plot.GMA"            "QTL"                
## [19] "QTLpeaks"            "rst.corR"            "rst.gwa"            
## [22] "sig.mrnaR"           "sig.mrnaR.names"     "sig.neg.mrnaR"      
## [25] "sig.neg.mrnaR.names" "sigpval"             "stb"                
## [28] "stb.nm"              "sum.eqtl"            "summary.sigR"       
## [31] "tbpos"               "zstandard"
```

## Analysis

Identify map positions of negatively correlated genes for miR-874:

First, calculate the midpoint of the gene:


```r
annot$pos<-annot$start+(round((annot$end - annot$start)/2))
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
##                   qvalue       pos
## XLOC_001644 3.894488e-02 120877013
## XLOC_002556 9.883916e-03  24550278
## XLOC_004063 3.613759e-02  36465364
## XLOC_014709 3.613759e-02  72260951
## XLOC_007579 5.830717e-05 128964105
## XLOC_002971 4.546732e-02  46881157
```

```r
annot<-annot[annot$miRNA=="ssc-miR-874",]
annot<-annot[annot$cor<0,]
annot<-annot[annot$chr!="AEMK02000682.1",]
annot$chr<-gsub("X", "19", annot$chr)
head(annot)
```

```
##             chr     start       end  width strand          ID    genes
## XLOC_000043   1   6771661   6891206 119546      + XLOC_000043   AGPAT4
## XLOC_000073   1  14794909  14851995  57087      + XLOC_000073    RMND1
## XLOC_000121   1  29835667  29904022  68356      + XLOC_000121  SLC2A12
## XLOC_000475   1 122138991 122394729 255739      + XLOC_000475  FAM227B
## XLOC_000500   1 127814360 127830785  16426      + XLOC_000500 CATSPER2
## XLOC_000649   1 162936246 163065872 129627      + XLOC_000649   ATP8B1
##                     locus       miRNA        cor      pvalue      qvalue
## XLOC_000043 RLOC_00000068 ssc-miR-874 -0.1693319 0.001202162 0.007901230
## XLOC_000073 RLOC_00000130 ssc-miR-874 -0.1469880 0.004937962 0.009874884
## XLOC_000121 RLOC_00000243 ssc-miR-874 -0.1519533 0.003660677 0.009441877
## XLOC_000475 RLOC_00000917 ssc-miR-874 -0.1450894 0.005524492 0.010400658
## XLOC_000500 RLOC_00000979 ssc-miR-874 -0.1053669 0.043896596 0.030371229
## XLOC_000649 RLOC_00001272 ssc-miR-874 -0.1342826 0.010226623 0.013962173
##                   pos
## XLOC_000043   6831433
## XLOC_000073  14823452
## XLOC_000121  29869845
## XLOC_000475 122266860
## XLOC_000500 127822572
## XLOC_000649 163001059
```

```r
dim(annot)
```

```
## [1] 332  13
```

miRNA being analyzed (manual input)


```r
X <- "ssc-miR-874"
```

Extract the results for the miRNA X


```r
miX <- negcoloc[negcoloc$miRNA == X,]
miX$chr<-as.numeric(as.character(miX$chr))
str(miX)
```

```
## 'data.frame':	52 obs. of  7 variables:
##  $ chr    : num  1 1 2 2 2 2 2 2 2 2 ...
##  $ ID     : chr  "XLOC_002346" "XLOC_002373" "XLOC_011766" "XLOC_011767" ...
##  $ genes  : Factor w/ 22233 levels "5_8S_RRNA","5S_RRNA",..: 21266 18022 1800 16858 16019 1972 19199 16019 1800 16858 ...
##  $ miRNA  : Factor w/ 5 levels "ssc-let-7d-5p",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ cor    : num  -0.103 -0.14 -0.12 -0.137 -0.155 ...
##  $ pheno  : Factor w/ 14 levels "car_bf10","WBS",..: 1 1 2 2 2 2 2 2 3 3 ...
##  $ cor.mag: chr  "-" "-" "-" "-" ...
```

```r
# Extracts the gene name, chromosome and position for all genes significantly corrlated to miRNA X
tmp <- lapply(unique(as.character(miX$genes)), function(x) coloc[coloc$ID %in%
	unique(miX[miX$genes == x,"ID"]) & as.character(coloc$miRNA) == X,
	c("ID", "chr", "start", "end", "qvalue")])
# Obtains the mid position of each target gene transcript
tmp <- lapply(tmp, function(x) unique(data.frame(geneID=x$ID, chr=x$chr,
	pos=unlist(lapply(1:nrow(x), function(y) mean(unlist(x[y,c("start", "end")])))),
	cor.qval=x$qvalue)))
targets.miX <- do.call(rbind, tmp)
targets.miX$chr<-as.numeric(as.character(targets.miX$chr))
rownames(targets.miX) <- targets.miX$geneID
targets.miX[order(targets.miX$chr),]
```

```
##                  geneID chr       pos    cor.qval
## XLOC_002346 XLOC_002346   1 271450290 0.031691123
## XLOC_002373 XLOC_002373   1 273796420 0.012162627
## XLOC_011766 XLOC_011766   2   5101364 0.019611261
## XLOC_011767 XLOC_011767   2   5123510 0.012795924
## XLOC_011794 XLOC_011794   2   6512744 0.009055118
## XLOC_013030 XLOC_013030   2   6532744 0.019474731
## XLOC_013026 XLOC_013026   2   6437333 0.018768508
## XLOC_013029 XLOC_013029   2   6481902 0.015775143
## XLOC_019568 XLOC_019568   6  69581599 0.008344526
## XLOC_021684 XLOC_021684   7  29654606 0.009441877
## XLOC_021691 XLOC_021691   7  29994014 0.009874884
## XLOC_021695 XLOC_021695   7  30345765 0.014626639
## XLOC_022575 XLOC_022575   7  28710899 0.007715510
## XLOC_022612 XLOC_022612   7  32095806 0.017144718
## XLOC_022617 XLOC_022617   7  32360737 0.048976188
## XLOC_022870 XLOC_022870   7  74858888 0.008344526
## XLOC_022875 XLOC_022875   7  74989130 0.009358773
## XLOC_023083 XLOC_023083   7  89146434 0.048976188
## XLOC_023143 XLOC_023143   7  98070058 0.025709409
## XLOC_023161 XLOC_023161   7 100399206 0.019256573
## XLOC_003117 XLOC_003117  11   3540408 0.007715510
## XLOC_008680 XLOC_008680  15  57458678 0.019897284
## XLOC_009274 XLOC_009274  15  57487704 0.008599413
## XLOC_008685 XLOC_008685  15  59639081 0.026845682
## XLOC_008850 XLOC_008850  15 104595833 0.007815367
## XLOC_008918 XLOC_008918  15 118859800 0.021676537
## XLOC_009276 XLOC_009276  15  59298208 0.040503591
## XLOC_009284 XLOC_009284  15  60411548 0.010970398
## XLOC_009349 XLOC_009349  15  77520240 0.031966317
## XLOC_009537 XLOC_009537  15 121339509 0.018249445
## XLOC_009607 XLOC_009607  15 137718665 0.016460629
```

Check that the number of targets matches the number of significantly correlated
target genes for miRNA X: should be TRUE


```r
length(table(as.character(miX$ID))) == nrow(targets.miX)
```

```
## [1] TRUE
```

Add the position of the significantly correlated target genes for miRNA X to the map.full object


```r
map.full <- rbind(map.full, targets.miX[,c("chr", "pos")])
map.full <- rbind(map.full, annot[,c("chr", "pos")])
map.full$chr<-as.numeric(map.full$chr)
```

Use the absmap function to convert the chromosomal map positions to absolute map positions:


```r
absposmap <- absmap(map.full)/1e6
```

Obtain phenotypic QTL markers only for the pQTL with a co-localized target gene


```r
pheno <- unique(as.character(miX$pheno))
sig.pheno <- lapply(pheno, function(x) GWAS$qvalue[GWAS$qvalue[,x]<0.05, x])
names(sig.pheno) <- pheno
```

Check markers associated to a pQTL that are not in the absposmap


```r
lapply(sig.pheno, function(x) names(x)[!names(x) %in% names(absposmap)])
```

```
## $car_bf10
## character(0)
## 
## $WBS
## character(0)
## 
## $tenderness
## [1] "ASGA0025952" "H3GA0016568" "H3GA0016570" "M1GA0007871" "ALGA0032493"
## 
## $overtend
## [1] "ALGA0032416" "ASGA0025952" "H3GA0016568" "H3GA0016570" "M1GA0007871"
## [6] "ALGA0032493"
## 
## $lma_16wk
## character(0)
## 
## $dress_ptg
## character(0)
## 
## $num_ribs
## [1] "INRA0027603"
## 
## $protein
## [1] "MARC0089139"
## 
## $ph_24h
## character(0)
## 
## $driploss
## character(0)
## 
## $cook_yield
## character(0)
## 
## $juiciness
## character(0)
```

```r
unlist(lapply(sig.pheno, length))
```

```
##   car_bf10        WBS tenderness   overtend   lma_16wk  dress_ptg 
##        102         13         20         25          4         72 
##   num_ribs    protein     ph_24h   driploss cook_yield  juiciness 
##         91         97         35         81         68          5
```

Remove markers associated to a pQTL that is not co-localized with the target genes of miRNA X


```r
sig.pheno <- lapply(pheno, function(x) sig.pheno[[x]][map.full[names(sig.pheno[[x]]), "chr"] %in%
	as.character(miX[as.character(miX$pheno) == x,"chr"])])
names(sig.pheno) <- pheno
unlist(lapply(sig.pheno, length))
```

```
##   car_bf10        WBS tenderness   overtend   lma_16wk  dress_ptg 
##          7         13         14         18          4         72 
##   num_ribs    protein     ph_24h   driploss cook_yield  juiciness 
##         90         96         35         81         65          5
```

Create list with targets co-localized with pQTL (used to plot pQTL names in manhatan plot)


```r
bychr <- lapply(unique(as.character(miX$chr)), function(x) miX[as.character(miX$chr) == x,])
bychr <- lapply(bychr, function(x) data.frame(chr=rep(unique(x$chr), length(unique(x$pheno))),
	pheno=unique(as.character(x$pheno)),
	do.call(rbind, lapply(unique(as.character(x$pheno)), function(y) data.frame(
		start=min(absposmap[names(sig.pheno[[y]])[map.full[names(sig.pheno[[y]]),"chr"]
			%in% unique(x$chr)]]),
		end=max(absposmap[names(sig.pheno[[y]])[map.full[names(sig.pheno[[y]]),"chr"]
			%in% unique(x$chr)]]),
		snp.qval=min(sig.pheno[[y]][map.full[names(sig.pheno[[y]]),"chr"] %in% unique(x$chr)]))))))
```

Manually specify were the text for each pQTL will be (left, right or top of peak)


```r
# Position next to peak
pos <- c("top", "top", "top", "left", "top", "right")
# Position of text next to peak (pos object in R plots)
Spos <- c("top", "top", "top", "left", "top", "right")
```

Generate the manhatan plot for the ssc-miR-874 eQTL


```r
qvals <- rst.gwa[rst.gwa$miRNA==X, "gwa.qval"]
names(qvals) <- rst.gwa[rst.gwa$miRNA==X, "SNPid"]
```

Position of miRNA X eQTL


```r
eqtlX <- map.full[names(qvals)[qvals < 0.05],]
pos_eqtlX <- unlist(lapply(unique(eqtlX$chr), function(x) paste("SSC", x, ":",
	ifelse(nrow(eqtlX[eqtlX$chr == x,]) == 1,
		paste(round(eqtlX[eqtlX$chr == x, "pos"]/1e6, 3), " Mb", sep=""),
		paste(round(min(eqtlX[eqtlX$chr == x, "pos"])/1e6, 3),
	"-", end=round(max(eqtlX[eqtlX$chr == x, "pos"])/1e6, 3), " Mb", sep="")), sep="")))
```

## Visualize

Plot eQTL and target genes co-localizing with pQTL


```r
par(oma=c(2,2,2,2))
manhpt(nm=X,abspos=absposmap,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,20))
```

```
## $transcript
## [1] "ssc-miR-874"
## 
## $max
## [1] 8.5415
## 
## $num.snp
## [1] 116
```

```r
points(absposmap[as.character(rownames(annot))], -log10(annot$qvalue), pch=20, col="grey60", cex=1)
lapply(pheno, function(x) points(absposmap[names(sig.pheno[[x]])],
	-log10(sig.pheno[[x]]), col="chocolate1", pch=1))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
## 
## [[10]]
## NULL
## 
## [[11]]
## NULL
## 
## [[12]]
## NULL
```

```r
# Highlight targets co-localized with pQTL
points(absposmap[as.character(targets.miX$geneID)], -log10(targets.miX$cor.qval), pch=20, col="black", cex=1.5)

# Highlight eQTL peak
text(x=0, y=19, labels="eQTL Chromosome: Position", font=2, pos=4, cex=2)
ifelse(length(pos_eqtlX) > 1, lapply(1:length(pos_eqtlX),
	function(i) text(x=0, y=19 - (i * 1.5), labels=pos_eqtlX[i], pos=4, cex=2)),
	text(x=0, y=19 - 1.5, labels=pos_eqtlX, pos=4, cex=2))
```

```
## [[1]]
## NULL
```

```r
text(absposmap[X], max(-log10(qvals)) + 1, labels="eQTL", font=2, cex=2)
# Add chromosome to pQTL peaks
lapply(1:length(bychr), function(z)
	text(x=ifelse(pos[z] == "left", min(bychr[[z]]$start) - 5, ifelse(pos[z] == "right",
			max(bychr[[z]]$end) + 10, min(bychr[[z]]$start))),
		y=ifelse(pos[z]=="top", max(-log10(bychr[[z]]$snp.qval)) + (nrow(bychr[[z]]) * 2.5),
			max(-log10(bychr[[z]]$snp.qval))),
	labels=paste("SSC", unique(bychr[[z]]$chr), ":", sep=""),
	pos=ifelse(Spos[z] == "right", 4, 2), font=2, cex=2))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
```

```r
# Add phenotype names to pQTL peaks
lapply(1:length(bychr), function(z) lapply(1:nrow(bychr[[z]]), function(i)
	text(x=ifelse(pos[z] == "left", min(bychr[[z]]$start) - 5, ifelse(pos[z] == "right",
			max(bychr[[z]]$end) + 10, min(bychr[[z]]$start))),
		y=ifelse(pos[z]=="top", (max(-log10(bychr[[z]]$snp.qval)) + (nrow(bychr[[z]]) * 2.5)) - (1 * i),
			max(-log10(bychr[[z]]$snp.qval)) - (1 * i)),
	labels=bychr[[z]][i,"pheno"], pos=ifelse(Spos[z] == "right", 4, 2), cex=2)))
```

```
## [[1]]
## [[1]][[1]]
## NULL
## 
## 
## [[2]]
## [[2]][[1]]
## NULL
## 
## [[2]][[2]]
## NULL
## 
## [[2]][[3]]
## NULL
## 
## 
## [[3]]
## [[3]][[1]]
## NULL
## 
## 
## [[4]]
## [[4]][[1]]
## NULL
## 
## [[4]][[2]]
## NULL
## 
## 
## [[5]]
## [[5]][[1]]
## NULL
## 
## 
## [[6]]
## [[6]][[1]]
## NULL
## 
## [[6]][[2]]
## NULL
## 
## [[6]][[3]]
## NULL
## 
## [[6]][[4]]
## NULL
## 
## [[6]][[5]]
## NULL
## 
## [[6]][[6]]
## NULL
## 
## [[6]][[7]]
## NULL
## 
## [[6]][[8]]
## NULL
```

```r
add_legend("bottomright",legend=c("miR-874 SNP Assoc.", "pQTL SNP Assoc.", "Target Genes", "Co-localized Targets"),
	pch=c(1,1), bty="o",ncol=2, cex=1.2)
	add_legend("bottomright",legend=c("miR-874 SNP Assoc.", "pQTL SNP Assoc.", "Target Genes", "Co-localized Targets"),
	pch=19, 
	col=c("deepskyblue", "chocolate1", "grey60", "black"), bty="n",ncol=2, cex=1.2)
```

<img src="figure/ssc-miR-874-1.tiff" title="plot of chunk ssc-miR-874" alt="plot of chunk ssc-miR-874" style="display: block; margin: auto;" />

