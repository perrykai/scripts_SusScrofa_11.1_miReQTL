**Script:** `9_test_coloc_mireqtl_pqtl.R`

**Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`

**Date:**  12-19-17

**Input File Directory:**  

1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/` `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`

2. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`

3. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects`

4. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl`

**Input File(s):** 

1. `pQTL_60k.Rdata`, `funct_eqtl.Rdata`

2. `MSUPRP_gpData_Ss11.Rdata`

3. `6_mirna_precursor_annot_ssc11.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`

4. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`

5. `9_mireqtl_pqtl_coloc_peaks.Rdata`, `13_Z_full.Rdata`

**Output File Directory:** 

`/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/7_miRNA_eQTL_colocalization/`

**Output File(s):** 

`14_sig_coloc.Rdata`

**Table of contents:**

1. [Objectives](#objectives)
2. [Install libraries](#install-libraries)
3. [Load data](#load-data)
4. [Analysis](#analysis)
5. [Visualize](#visualize)
6. [Save data](#save-data)

## Objectives

Given the identified co-localized miR-eQTL with pQTL, test the the effect of the top pQTL marker
on the expression of its colocalized miR-eQTL. Results include a dataframe indicating if any significant SNPs remain after fixing peak SNP, 
and the amount & proportion of variance explained by the top SNP.

## Install libraries

Clear environment


```r
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")
```

Required Packages


```r
library(limma)
library(edgeR)
library(gwaR)
library(regress)
library(qvalue)
library(methods)
```

## Load data


```r
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
```

Load DV functions:


```r
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
load("../../5_gblup_gwa_eqtl/2_gwa_results.Rdata")
load("../../5_gblup_gwa_eqtl/3_eqtl_summary_tables_maps.Rdata")
load("../9_mireqtl_pqtl_coloc_peaks.Rdata")
load("../13_Z_full.Rdata")

ls()
```

```
##  [1] "absmap"               "absposmap"            "add_legend"          
##  [4] "AddPosGene"           "annotation"           "dge"                 
##  [7] "distance"             "fullsum.eqtl"         "G"                   
## [10] "GBLUP"                "GWAS"                 "inrange"             
## [13] "manhpt"               "map.full"             "mirpeaks"            
## [16] "MSUPRP"               "MSUPRP168"            "MSUPRP_miRNA"        
## [19] "peakrng"              "plot.GMA"             "QTL"                 
## [22] "QTLpeaks"             "regul"                "rst.gwa"             
## [25] "sigpval"              "stb"                  "stb.nm"              
## [28] "sum.eqtl"             "summary_MSUPRP_miRNA" "tbpos"               
## [31] "v"                    "wtcen"                "Z"                   
## [34] "Zfull"                "zstandard"
```

## Analysis

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

Full Marker Map


```r
mapA <- data.frame(chr=MSUPRP$map$chr, pos=MSUPRP$map$pos)
rownames(mapA) <- colnames(MSUPRP$geno)
mapA$chr<-gsub("X", "19", mapA$chr)
mapA$chr<-gsub("Y", "20", mapA$chr)
mapA$chr<-as.numeric(mapA$chr)
mapA$pos <- mapA$pos/1e6
dim(mapA)
```

```
## [1] 43130     2
```

Retain marker information for all pQTL


```r
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(mapA[names(sig[[x]]),],
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

Get pQTL peak positions


```r
pos.pqtl <- lapply(sig, function(x) do.call(rbind, lapply(unique(x$chr),
	function(y) data.frame(chr=y, min=min(x[x$chr == y,"pos"]), max=max(x[x$chr == y,"pos"])))))

# Total number of markers in each pQTL peak
snpP <- lapply(pos.pqtl, function(x) unlist(lapply(1:nrow(x),
			function(y) sum(mapA[mapA[,"chr"] == x[y,"chr"],"pos"] >= x[y,"min"] &
			mapA[mapA[,"chr"] == x[y,"chr"],"pos"] <= x[y,"max"]))))
for(i in 1:length(snpP)){
	names(snpP[[i]]) <- unique(pos.pqtl[[i]]$chr)
}

snpP
```

```
## $bf10_10wk
##  6 
## 36 
## 
## $lrf_10wk
##   6  12 
## 576   1 
## 
## $bf10_13wk
##   3   6 
##   1 109 
## 
## $lrf_13wk
##   6 
## 625 
## 
## $bf10_16wk
## 6 
## 1 
## 
## $lma_16wk
## 6 
## 6 
## 
## $lrf_16wk
##   5   6 
##  25 148 
## 
## $bf10_19wk
## 6 
## 1 
## 
## $lrf_19wk
##  6 
## 38 
## 
## $bf10_22wk
##  6 
## 35 
## 
## $lrf_22wk
##   6 
## 187 
## 
## $dress_ptg
##   7  11 
## 370 114 
## 
## $cook_yield
##  15   5   8 
## 470  74   1 
## 
## $WBS
##  2 15 
## 73 77 
## 
## $juiciness
## 15 
## 54 
## 
## $tenderness
##  2  3  5 15 
## 73  1 10 77 
## 
## $overtend
##   2   3   5  15 
##  33   1  21 119 
## 
## $driploss
##  15 
## 612 
## 
## $ph_24h
##  15 
## 484 
## 
## $car_length
## 7 
## 3 
## 
## $num_ribs
##   7 
## 995 
## 
## $last_lum
##   4   6   9  10 
##   1 150   2   1 
## 
## $car_bf10
##   1   6 
## 106 661 
## 
## $loin
##   6  11 
## 198   1 
## 
## $protein
##   15 
## 1214
```

Phenotypic QTL co-localized with expression QTL


```r
coloc <- regul[!regul$colocalized.pqtl == "",]
coloc
```

```
##    chr.snp         SNP   pos.snp     pval.snp     qval.snp   min.pos
## 2        4 ALGA0026452  87026192 4.686046e-07 1.700660e-02  87026192
## 1        7 ASGA0034057  50959933 3.850651e-11 1.996397e-07  43560865
## 9       15 H3GA0052416 121806256 1.705422e-06 3.581458e-02 121806256
## 7       15 MARC0093624 122218534 3.093005e-07 1.122513e-02 121872813
## 8       15 MARC0093624 122218534 1.643296e-07 5.963851e-03 121872813
## 10      15 MARC0093624 122218534 2.936732e-06 4.636241e-02 121806256
## 12      15 MARC0093624 122218534 8.887961e-08 3.225619e-03 121806256
## 16      15 MARC0093624 122218534 5.550876e-07 2.014524e-02 122218534
##      max.pos range.peak num.snp           miRNA chr.miR start.miR
## 2  101458668   14432476       4    ssc-miR-190b       4  95540606
## 1   86790528   43229663      46     ssc-miR-184       7  48345017
## 9  121872813      66557       2  ssc-miR-345-3p       7 121193716
## 7  122218534     345721       2   ssc-let-7d-5p       3  43468501
## 8  122218534     345721       2      ssc-let-7g      13  34406135
## 10 122218534     412278       3      ssc-miR-95       8   3030934
## 12 122218534     412278       3 ssc-miR-9843-3p       8 114110660
## 16 122218534          0       1    ssc-miR-1468      19  50337088
##      end.miR range.miR         h2    lrtpvalue       qvalue regulator
## 2   95540688        82 0.50100208 3.441217e-05 1.700660e-02       cis
## 1   48345099        82 0.63092373 4.682400e-09 1.996397e-07       cis
## 9  121193799        83 0.41455537 1.474400e-05 3.581458e-02      tran
## 7   43468596        95 0.11285941 4.512942e-02 1.122513e-02      tran
## 8   34406222        87 0.09900429 5.157767e-02 5.963851e-03      tran
## 10   3031014        80 0.41584487 2.144450e-07 4.636241e-02      tran
## 12 114110740        80 0.31124288 3.548349e-04 3.225619e-03      tran
## 16  50337170        82 0.19263335 3.002969e-03 2.014524e-02      tran
##                                                               colocalized.pqtl
## 2                                                                     last_lum
## 1                                                                     num_ribs
## 9  juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
## 7  juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
## 8  juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
## 10 juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
## 12 juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
## 16 juiciness, driploss, ph_24h, protein, cook_yield, WBS, tenderness, overtend
```

 eQTL colocalized with meat quality traits only


```r
mq <- coloc[grep("juiciness",coloc$colocalized.pqtl),]
```

SNP associated to eQTL co-localized with pQTL


```r
# xloc ids
idx <- coloc$miRNA
length(idx)
```

```
## [1] 8
```

```r
snp.col <- lapply(idx, function(x) fullsum.eqtl[fullsum.eqtl$miRNA == x & fullsum.eqtl$chr.snp == coloc[coloc$miRNA == x, "chr.snp"],])
names(snp.col) <- idx

# Add pvalues (GWAS results including PRKAG3 SNP)
for(i in names(snp.col)){
	x<-rst.gwa[rst.gwa$miRNA == i,]

	snp.col[[i]]$pvalue<-x[match(as.character(snp.col[[i]]$SNP), as.character(x$SNPid)), "gwa.pval"]
}
```

List of colocalized pQTL (phenotype name) per gene expression


```r
pheno <- lapply(idx,
	function(x) strsplit(as.character(coloc[grep(x, coloc$miRNA),"colocalized.pqtl"]),", ")[[1]])
names(pheno) <- idx
```

Check number of SNPs in common between eQTL co-localized with pQTL


```r
com.snp <- lapply(idx, function(x) unlist(lapply(pheno[[x]], function(y)
	sum(rownames(sig[[y]]) %in% as.character(snp.col[[x]]$SNP)))))
names(com.snp) <- idx
for (i in idx){
	names(com.snp[[i]]) <- pheno[[i]]
}
check <- unlist(lapply(com.snp, function(x) sum(x)))
check[check > 0]
```

```
##     ssc-miR-184  ssc-miR-345-3p   ssc-let-7d-5p      ssc-let-7g 
##               7              16              16              16 
##      ssc-miR-95 ssc-miR-9843-3p    ssc-miR-1468 
##              24              24               8
```

Data: response, covariates, snp


```r
# Covariates
X <- data.frame(sex=as.factor(MSUPRP_miRNA$covar$sex),
	 selcrit=as.factor(paste(MSUPRP_miRNA$covar$Status, MSUPRP_miRNA$covar$selcrit,
	 sep="-")))
rownames(X) <- MSUPRP_miRNA$covar$id

# Data Frame
dataM <- data.frame(t(v$E[idx,rownames(X)]), Zfull[rownames(X),], X)
dim(dataM)
```

```
## [1]   174 42705
```

### Function to test the effect of fixing a SNP on a known pQTL

* rsp = response variable

* qtl = data frame or list of data frames containing information on associated SNPs per pQTL (chromosome, position, snp effect, pvalue, qvalue). If snps is TRUE the qtl is a scalar or vector with the names of the SNP to test in the GBLUP model

* eqtl = matrix with information on associated SNPs for response eQTL (gene, SNP, chromosome snp, position snp, pvalue, qvalue) Note: do not include snps that are not on the same chromosome as the colocalized pQTL

* coloc = scalar or vector containing the name of the co-localized QTL (must match the names in the qtl list), if snp is TRUE the

* design = model design

* dataM = matrix containg the data for the response, covariates and fixed effects

* G = genomic relationship matrix

* vdata = list containing additional random effects (default NULL)

* wt = matrix of weights (use the matrix of precision weights to model the mean variance relationship of gene expression profiles)

* top = logical, True: fix only the top significant marker in the QTL peak (TRUE, default), False: test the effect of all the markers in the co-localized QTL region

* Z = standardized SNP matrix

* fdr = false discovery rate cutoff (default 0.01)

* snps = logical, TRUE if the SNPs to be tested are known (qtl should contain the SNP names), FALSE if SNPs are to be selected from the QTL matrix

* ... = additional parameters to pass to the gblup function in gwaR

### Annotation of code:

1. loop through coloc object; if snps=TRUE, snp comes from vector of names of SNPs to be fixed (qtl); if snps=FALSE, take first element of qtl list and extract SNP names.

2. Check that there are only SNP in eQTL from one chromosome, and that it's the same chromosome as the pQTL.

3. If only testing top SNP, take SNP from qtl with minimum pvalue, if not (top=FALSE), take all SNPs from qtl to test.

4. If testing multiple SNPs (all colocalized SNPs), use loop to separate analyses by SNP in region.

5. Build design/model including tested SNP, run gblup

6. Extract SNP being tested from Z matrix, run gwas, get pvalues

7. Remove any null/NA pvalues obtained from gwas

8. Get qvalues, and qvalue for SNP being tested; remove NAs.

9. Standardize the SNP effect, determine how much variance the SNP explains & prop. of variance explained by SNP.

10. Build pv.anova object, which summarizes the variance explained by the tested SNP, if that amount is significant, the gblup variance components, and how many sig. SNPs remain after fixing top SNP.

11. Return these results!


### coloc.test Function:


```r
coloc.test <- function(rsp, qtl, eqtl, coloc, design, dataM, G, vdata=NULL, wt, top=TRUE, Z, fdr=0.01, snps=FALSE, ...){
	pv.anova <- NULL

	for (i in coloc){

		if (snps == TRUE){

			snp <- qtl

		} else {

			ifelse(is.list(qtl), tmp <- qtl[[i]], tmp <- qtl)

			tmp <- tmp[tmp$chr == unique(eqtl$chr.snp),]

			ifelse(top==TRUE, tmp<-tmp[tmp$pvalue == min(tmp$pvalue),], tmp<-tmp)

			snp <- rownames(tmp)
		}

		if(length(snp) > 1){
			for(j in snp){
				cat(paste("testing effect of sub ",j," for ",i,"...",sep=""),"\n")

			#' Design
			ifelse(is.list(design), design2 <- c(paste("~", design[[1]], " + ", j, sep="")[2],design[[2]]),
				design2 <- paste("~", design, " + ", j, sep="")[2])
			design2 <- as.formula(design2)

			gb <- gblup(rsp=rsp, data=dataM,
				design=design2, G=G, vdata=vdata, wt=wt, ...)

			ifelse(j %in% rownames(Z), z <- Z[-(grep(j, rownames(Z))),], z <- Z)

			gwa <- gwas(gblup=gb, x=z)
			pval <- getpvalue(gwa, log.p=F, is.z=F)

			if(sum(is.na(pval)) > 0){
				pval <- pval[!is.na(pval)]
			}

			qval <- qvalue(pval)$qvalue
			qvsnp <- qval[as.character(eqtl$SNP)][!is.na(qval[as.character(eqtl$SNP)])]

			VS <- gb$coef[j,"Estimate"]^2 * var(dataM[,j])
			propVS <- (VS) / sum(VS, gb$sigma)

			pv.anova[[i]][[j]] <- data.frame(snp=j,anova(gb)[j,], varG=gb$sigma["G"],
				varE=gb$sigma[2], varS=VS, propVS=propVS, sigsnp=sum(qvsnp < fdr))
			}
			pv.anova[[i]] <- do.call(rbind, pv.anova[[i]])

		} else {

		cat(paste("testing effect of ",snp," for ",i,"...",sep=""),"\n")

		#' Design
		ifelse(is.list(design), design2 <- c(paste("~", design[[1]], " + ", snp, sep="")[2], design[[2]]),
			design2 <- paste("~", design, " + ", snp, sep="")[2])
		design2 <- as.formula(design2)

		gb <- gblup(rsp=rsp, data=dataM,
				design=design2, G=G, vdata=vdata, wt=wt, ...)

		ifelse(snp %in% rownames(Z), z <- Z[-(grep(snp, rownames(Z))),], z <- Z)

		gwa <- gwas(gblup=gb, x=z)
		pval <- getpvalue(gwa, log.p=F, is.z=F)

		if(sum(is.na(pval)) > 0){
			pval <- pval[!is.na(pval)]
		}

		qval <- qvalue(pval)$qvalue
		qvsnp <- qval[as.character(eqtl$SNP)][!is.na(qval[as.character(eqtl$SNP)])]

		VS <- gb$coef[snp,"Estimate"]^2 * var(dataM[,snp])
		propVS <- (VS) / sum(VS, gb$sigma)

		pv.anova[[i]] <- data.frame(snp=snp,anova(gb)[snp,], varG=gb$sigma["G"],
			varE=gb$sigma[2], varS=VS, propVS=propVS, sigsnp=sum(qvsnp < fdr))
		}
	}
	pv <- do.call(rbind, pv.anova)
	return(pv)
}

idx<-gsub("-",".", idx)
names(pheno)<-idx
names(com.snp)<-idx
names(snp.col)<-idx
colnames(wtcen)<-gsub("-",".", colnames(wtcen))

design <- ~sex + selcrit
```

Co-localization test: Fix top peak marker per pQTL and fit GWAS model per gene expression


```r
testCo <- lapply(idx, function(x)
	coloc.test(rsp=x, qtl=sig, eqtl=snp.col[[x]], coloc=pheno[[x]], top=TRUE,
		design=design, dataM=dataM, G=G, wt=wtcen, Z=t(Zfull), fdr=0.05, pos=c(T,T)))
```

```
## testing effect of ASGA0092651 for last_lum... 
## testing effect of ALGA0043983 for num_ribs... 
## testing effect of H3GA0052416 for juiciness... 
## testing effect of MARC0093624 for driploss... 
## testing effect of MARC0093624 for ph_24h... 
## testing effect of MARC0093624 for protein... 
## testing effect of MARC0093624 for cook_yield... 
## testing effect of MARC0093624 for WBS... 
## testing effect of H3GA0052416 for tenderness... 
## testing effect of H3GA0052416 for overtend... 
## testing effect of H3GA0052416 for juiciness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for driploss... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for ph_24h... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for protein... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for cook_yield... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for WBS... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for tenderness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for overtend... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for juiciness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for driploss... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for ph_24h... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for protein... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for cook_yield... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for WBS... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for tenderness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for overtend... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for juiciness... 
## testing effect of MARC0093624 for driploss... 
## testing effect of MARC0093624 for ph_24h... 
## testing effect of MARC0093624 for protein... 
## testing effect of MARC0093624 for cook_yield... 
## testing effect of MARC0093624 for WBS... 
## testing effect of H3GA0052416 for tenderness... 
## testing effect of H3GA0052416 for overtend... 
## testing effect of H3GA0052416 for juiciness... 
## testing effect of MARC0093624 for driploss... 
## testing effect of MARC0093624 for ph_24h... 
## testing effect of MARC0093624 for protein... 
## testing effect of MARC0093624 for cook_yield... 
## testing effect of MARC0093624 for WBS... 
## testing effect of H3GA0052416 for tenderness... 
## testing effect of H3GA0052416 for overtend... 
## testing effect of H3GA0052416 for juiciness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for driploss... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for ph_24h... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for protein... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for cook_yield... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of MARC0093624 for WBS... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for tenderness... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
## testing effect of H3GA0052416 for overtend... 
## Warning: solution lies close to zero for some positive variance components, their standard errors may not be valid
```

```r
names(testCo) <- idx
```

Summarize results


```r
rstCo <- do.call(rbind, lapply(idx, function(x) data.frame(gene=x,
	pheno=unlist(lapply(rownames(testCo[[x]]), function(y) strsplit(y,"[.]")[[1]][1])),
	testCo[[x]][,c("snp", "p.value", "propVS","sigsnp")])))
rownames(rstCo) <- NULL

sigCo <- rstCo[rstCo$sigsnp == 0,]
sigEff <- rstCo[rstCo$p.value < 0.05 & rstCo$sigsnp > 0,]

table(as.character(sigCo$gene))
```

```
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468    ssc.miR.190b 
##               8               8               8               1 
##  ssc.miR.345.3p      ssc.miR.95 ssc.miR.9843.3p 
##               8               8               8
```

```r
rownames(dge$genes)<-gsub("-",".", rownames(dge$genes))
dge$genes[as.character(unique(sigCo$gene)),]
```

```
##                            Name  chr0     start       end width strand
## ssc.miR.190b       ssc-miR-190b  chr4 104457820 104457841    22      +
## ssc.miR.345.3p   ssc-miR-345-3p  chr7 128658348 128658368    21      +
## ssc.let.7d.5p     ssc-let-7d-5p  chr3  44867277  44867298    22      +
## ssc.let.7g           ssc-let-7g chr13  37599007  37599028    22      +
## ssc.miR.95           ssc-miR-95  chr8   4277286   4277307    22      +
## ssc.miR.9843.3p ssc-miR-9843-3p  chr8 122371914 122371935    22      -
## ssc.miR.1468       ssc-miR-1468  chrX  56757096  56757117    22      -
##                  type        Alias Precursors
## ssc.miR.190b    miRNA MIMAT0020588  MI0017988
## ssc.miR.345.3p  miRNA MIMAT0013900  MI0013117
## ssc.let.7d.5p   miRNA MIMAT0025356  MI0022120
## ssc.let.7g      miRNA MIMAT0013867  MI0013087
## ssc.miR.95      miRNA MIMAT0002142  MI0002436
## ssc.miR.9843.3p miRNA MIMAT0037061  MI0031612
## ssc.miR.1468    miRNA MIMAT0025386  MI0022160
```

Proportion of variance explained by top pQTL SNP


```r
sigComiR <- lapply(as.character(unique(sigCo$gene)), function(x)
	table(as.character(sigCo$pheno)[grep(x, sigCo$gene)]))
names(sigComiR) <- as.character(unique(sigCo$gene))
summary(sigCo$propVS*100)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   13.11   15.80   18.08   17.46   18.66   21.96
```

```r
sigEffmiR <- lapply(as.character(unique(sigEff$gene)), function(x)
	table(as.character(sigEff$pheno)[grep(x, sigEff$gene)]))
names(sigEffmiR) <- as.character(unique(sigEff$gene))
summary(sigEff$propVS*100)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   8.796   8.796   8.796   8.796   8.796   8.796
```

Investigate correlation of peak SNP when duplicates exist in miR-pheno combinations:



```r
lapply(unique(as.character(sigCo$pheno)), function(x) sigCo[sigCo$pheno==x,])
```

```
## [[1]]
##           gene    pheno         snp      p.value    propVS sigsnp
## 1 ssc.miR.190b last_lum ASGA0092651 8.688955e-07 0.1566973      0
## 
## [[2]]
##               gene     pheno         snp      p.value    propVS sigsnp
## 3   ssc.miR.345.3p juiciness H3GA0052416 3.794223e-08 0.1821112      0
## 11   ssc.let.7d.5p juiciness H3GA0052416 7.257963e-07 0.1409226      0
## 19      ssc.let.7g juiciness H3GA0052416 1.505415e-07 0.1399211      0
## 27      ssc.miR.95 juiciness H3GA0052416 4.589335e-07 0.1580258      0
## 35 ssc.miR.9843.3p juiciness H3GA0052416 6.128928e-09 0.1866145      0
## 43    ssc.miR.1468 juiciness H3GA0052416 2.628575e-06 0.1310848      0
## 
## [[3]]
##               gene    pheno         snp      p.value    propVS sigsnp
## 4   ssc.miR.345.3p driploss MARC0093624 1.973029e-07 0.1739949      0
## 12   ssc.let.7d.5p driploss MARC0093624 2.721450e-08 0.1811097      0
## 20      ssc.let.7g driploss MARC0093624 7.936065e-10 0.1807520      0
## 28      ssc.miR.95 driploss MARC0093624 2.553800e-07 0.1743376      0
## 36 ssc.miR.9843.3p driploss MARC0093624 9.750118e-10 0.2196241      0
## 44    ssc.miR.1468 driploss MARC0093624 1.023683e-08 0.1866595      0
## 
## [[4]]
##               gene  pheno         snp      p.value    propVS sigsnp
## 5   ssc.miR.345.3p ph_24h MARC0093624 1.973029e-07 0.1739949      0
## 13   ssc.let.7d.5p ph_24h MARC0093624 2.721450e-08 0.1811097      0
## 21      ssc.let.7g ph_24h MARC0093624 7.936065e-10 0.1807520      0
## 29      ssc.miR.95 ph_24h MARC0093624 2.553800e-07 0.1743376      0
## 37 ssc.miR.9843.3p ph_24h MARC0093624 9.750118e-10 0.2196241      0
## 45    ssc.miR.1468 ph_24h MARC0093624 1.023683e-08 0.1866595      0
## 
## [[5]]
##               gene   pheno         snp      p.value    propVS sigsnp
## 6   ssc.miR.345.3p protein MARC0093624 1.973029e-07 0.1739949      0
## 14   ssc.let.7d.5p protein MARC0093624 2.721450e-08 0.1811097      0
## 22      ssc.let.7g protein MARC0093624 7.936065e-10 0.1807520      0
## 30      ssc.miR.95 protein MARC0093624 2.553800e-07 0.1743376      0
## 38 ssc.miR.9843.3p protein MARC0093624 9.750118e-10 0.2196241      0
## 46    ssc.miR.1468 protein MARC0093624 1.023683e-08 0.1866595      0
## 
## [[6]]
##               gene      pheno         snp      p.value    propVS sigsnp
## 7   ssc.miR.345.3p cook_yield MARC0093624 1.973029e-07 0.1739949      0
## 15   ssc.let.7d.5p cook_yield MARC0093624 2.721450e-08 0.1811097      0
## 23      ssc.let.7g cook_yield MARC0093624 7.936065e-10 0.1807520      0
## 31      ssc.miR.95 cook_yield MARC0093624 2.553800e-07 0.1743376      0
## 39 ssc.miR.9843.3p cook_yield MARC0093624 9.750118e-10 0.2196241      0
## 47    ssc.miR.1468 cook_yield MARC0093624 1.023683e-08 0.1866595      0
## 
## [[7]]
##               gene pheno         snp      p.value    propVS sigsnp
## 8   ssc.miR.345.3p   WBS MARC0093624 1.973029e-07 0.1739949      0
## 16   ssc.let.7d.5p   WBS MARC0093624 2.721450e-08 0.1811097      0
## 24      ssc.let.7g   WBS MARC0093624 7.936065e-10 0.1807520      0
## 32      ssc.miR.95   WBS MARC0093624 2.553800e-07 0.1743376      0
## 40 ssc.miR.9843.3p   WBS MARC0093624 9.750118e-10 0.2196241      0
## 48    ssc.miR.1468   WBS MARC0093624 1.023683e-08 0.1866595      0
## 
## [[8]]
##               gene      pheno         snp      p.value    propVS sigsnp
## 9   ssc.miR.345.3p tenderness H3GA0052416 3.794223e-08 0.1821112      0
## 17   ssc.let.7d.5p tenderness H3GA0052416 7.257963e-07 0.1409226      0
## 25      ssc.let.7g tenderness H3GA0052416 1.505415e-07 0.1399211      0
## 33      ssc.miR.95 tenderness H3GA0052416 4.589335e-07 0.1580258      0
## 41 ssc.miR.9843.3p tenderness H3GA0052416 6.128928e-09 0.1866145      0
## 49    ssc.miR.1468 tenderness H3GA0052416 2.628575e-06 0.1310848      0
## 
## [[9]]
##               gene    pheno         snp      p.value    propVS sigsnp
## 10  ssc.miR.345.3p overtend H3GA0052416 3.794223e-08 0.1821112      0
## 18   ssc.let.7d.5p overtend H3GA0052416 7.257963e-07 0.1409226      0
## 26      ssc.let.7g overtend H3GA0052416 1.505415e-07 0.1399211      0
## 34      ssc.miR.95 overtend H3GA0052416 4.589335e-07 0.1580258      0
## 42 ssc.miR.9843.3p overtend H3GA0052416 6.128928e-09 0.1866145      0
## 50    ssc.miR.1468 overtend H3GA0052416 2.628575e-06 0.1310848      0
```

```r
lapply(unique(as.character(sigCo$pheno)), function(x) table(as.character(sigCo[sigCo$pheno==x,"snp"])))
```

```
## [[1]]
## 
## ASGA0092651 
##           1 
## 
## [[2]]
## 
## H3GA0052416 
##           6 
## 
## [[3]]
## 
## MARC0093624 
##           6 
## 
## [[4]]
## 
## MARC0093624 
##           6 
## 
## [[5]]
## 
## MARC0093624 
##           6 
## 
## [[6]]
## 
## MARC0093624 
##           6 
## 
## [[7]]
## 
## MARC0093624 
##           6 
## 
## [[8]]
## 
## H3GA0052416 
##           6 
## 
## [[9]]
## 
## H3GA0052416 
##           6
```

```r
lapply(unique(as.character(sigCo$pheno)), function(x) table(as.character(sigCo[sigCo$pheno==x,"gene"])))
```

```
## [[1]]
## 
## ssc.miR.190b 
##            1 
## 
## [[2]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[3]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[4]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[5]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[6]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[7]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[8]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1 
## 
## [[9]]
## 
##   ssc.let.7d.5p      ssc.let.7g    ssc.miR.1468  ssc.miR.345.3p 
##               1               1               1               1 
##      ssc.miR.95 ssc.miR.9843.3p 
##               1               1
```

## Save data


```r
save(testCo, sigCo, sigEff, file="../6_sig_coloc.Rdata")
```

