#' **Script:** `9_test_coloc_mireqtl_pqtl.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12-19-17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/` `/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/`
#' 
#' 2. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`
#' 
#' 3. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects`
#' 
#' 4. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl`
#' 
#' **Input File(s):** 
#' 
#' 1. `pQTL_60k.Rdata`, `funct_eqtl.Rdata`
#' 
#' 2. `MSUPRP_gpData_Ss11.Rdata`
#' 
#' 3. `6_mirna_precursor_annot_ssc11.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`
#' 
#' 4. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 5. `9_mireqtl_pqtl_coloc_peaks.Rdata`, `13_Z_full.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' `/mnt/research/pigeqtl/analyses/microRNA/2_mirna_characterization_expression/7_miRNA_eQTL_colocalization/`
#' 
#' **Output File(s):** 
#' 
#' `14_sig_coloc.Rdata`
#' 
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Visualize](#visualize)
#' 6. [Save data](#save-data)
#' 
#' ## Objectives
#' 
#' Given the identified co-localized miR-eQTL with pQTL, test the the effect of the top pQTL marker
#' on the expression of its colocalized miR-eQTL. Results include a dataframe indicating if any significant SNPs remain after fixing peak SNP, 
#' and the amount & proportion of variance explained by the top SNP.
#' 

#' ## Install libraries
#' 
#' Clear environment
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

#' Required Packages
library(limma)
library(edgeR)
library(gwaR)
library(regress)
library(qvalue)
library(methods)

#' ## Load data
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")

#' Load DV functions:
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
#' ## Analysis
#' 
#' ### Phenotypic QTL MSUPRP
length(QTL)
names(QTL)
QTL[[1]]

#' pQTL peak information
dim(QTLpeaks)
QTLpeaks$pheno.chr<-as.factor(make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE))
rownames(QTLpeaks)<-make.names(paste(QTLpeaks$pheno, QTLpeaks$chr, sep="."), unique=TRUE)
head(QTLpeaks)
table(QTLpeaks$pheno)
table(QTLpeaks$chr)

#' q-values pQTL GWAS
qval <- GWAS$qvalue
dim(qval)

#' p-values pQTL GWAS
pval <- GWAS$pvalue
dim(pval)

#' Standardized SNP effects pQTL GWA
sdEff <- GWAS$Estimate
dim(sdEff)


#' Full Marker Map
mapA <- data.frame(chr=MSUPRP$map$chr, pos=MSUPRP$map$pos)
rownames(mapA) <- colnames(MSUPRP$geno)
mapA$chr<-gsub("X", "19", mapA$chr)
mapA$chr<-gsub("Y", "20", mapA$chr)
mapA$chr<-as.numeric(mapA$chr)
mapA$pos <- mapA$pos/1e6
dim(mapA)

#' Retain marker information for all pQTL
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(mapA[names(sig[[x]]),],
	std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]
length(sig)
names(sig)

#' Get pQTL peak positions
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

#' Phenotypic QTL co-localized with expression QTL
coloc <- regul[!regul$colocalized.pqtl == "",]
coloc
#'  eQTL colocalized with meat quality traits only
mq <- coloc[grep("juiciness",coloc$colocalized.pqtl),]

#' SNP associated to eQTL co-localized with pQTL
# xloc ids
idx <- coloc$miRNA
length(idx)
snp.col <- lapply(idx, function(x) fullsum.eqtl[fullsum.eqtl$miRNA == x & fullsum.eqtl$chr.snp == coloc[coloc$miRNA == x, "chr.snp"],])
names(snp.col) <- idx

# Add pvalues (GWAS results including PRKAG3 SNP)
for(i in names(snp.col)){
	x<-rst.gwa[rst.gwa$miRNA == i,]

	snp.col[[i]]$pvalue<-x[match(as.character(snp.col[[i]]$SNP), as.character(x$SNPid)), "gwa.pval"]
}


#' List of colocalized pQTL (phenotype name) per gene expression
pheno <- lapply(idx,
	function(x) strsplit(as.character(coloc[grep(x, coloc$miRNA),"colocalized.pqtl"]),", ")[[1]])
names(pheno) <- idx

#' Check number of SNPs in common between eQTL co-localized with pQTL
com.snp <- lapply(idx, function(x) unlist(lapply(pheno[[x]], function(y)
	sum(rownames(sig[[y]]) %in% as.character(snp.col[[x]]$SNP)))))
names(com.snp) <- idx
for (i in idx){
	names(com.snp[[i]]) <- pheno[[i]]
}
check <- unlist(lapply(com.snp, function(x) sum(x)))
check[check > 0]

#' Data: response, covariates, snp
# Covariates
X <- data.frame(sex=as.factor(MSUPRP_miRNA$covar$sex),
	 selcrit=as.factor(paste(MSUPRP_miRNA$covar$Status, MSUPRP_miRNA$covar$selcrit,
	 sep="-")))
rownames(X) <- MSUPRP_miRNA$covar$id

# Data Frame
dataM <- data.frame(t(v$E[idx,rownames(X)]), Zfull[rownames(X),], X)
dim(dataM)

#' ### Function to test the effect of fixing a SNP on a known pQTL
#' 
#' * rsp = response variable
#' 
#' * qtl = data frame or list of data frames containing information on associated SNPs per pQTL (chromosome, position, snp effect, pvalue, qvalue). If snps is TRUE the qtl is a scalar or vector with the names of the SNP to test in the GBLUP model
#' 
#' * eqtl = matrix with information on associated SNPs for response eQTL (gene, SNP, chromosome snp, position snp, pvalue, qvalue) Note: do not include snps that are not on the same chromosome as the colocalized pQTL
#' 
#' * coloc = scalar or vector containing the name of the co-localized QTL (must match the names in the qtl list), if snp is TRUE the
#' 
#' * design = model design
#' 
#' * dataM = matrix containg the data for the response, covariates and fixed effects
#' 
#' * G = genomic relationship matrix
#' 
#' * vdata = list containing additional random effects (default NULL)
#' 
#' * wt = matrix of weights (use the matrix of precision weights to model the mean variance relationship of gene expression profiles)
#' 
#' * top = logical, True: fix only the top significant marker in the QTL peak (TRUE, default), False: test the effect of all the markers in the co-localized QTL region
#' 
#' * Z = standardized SNP matrix
#' 
#' * fdr = false discovery rate cutoff (default 0.01)
#' 
#' * snps = logical, TRUE if the SNPs to be tested are known (qtl should contain the SNP names), FALSE if SNPs are to be selected from the QTL matrix
#' 
#' * ... = additional parameters to pass to the gblup function in gwaR
#' 
#' ### Annotation of code:
#' 
#' 1. loop through coloc object; if snps=TRUE, snp comes from vector of names of SNPs to be fixed (qtl); if snps=FALSE, take first element of qtl list and extract SNP names.
#' 
#' 2. Check that there are only SNP in eQTL from one chromosome, and that it's the same chromosome as the pQTL.
#' 
#' 3. If only testing top SNP, take SNP from qtl with minimum pvalue, if not (top=FALSE), take all SNPs from qtl to test.
#'
#' 4. If testing multiple SNPs (all colocalized SNPs), use loop to separate analyses by SNP in region.
#' 
#' 5. Build design/model including tested SNP, run gblup
#' 
#' 6. Extract SNP being tested from Z matrix, run gwas, get pvalues
#' 
#' 7. Remove any null/NA pvalues obtained from gwas
#' 
#' 8. Get qvalues, and qvalue for SNP being tested; remove NAs.
#' 
#' 9. Standardize the SNP effect, determine how much variance the SNP explains & prop. of variance explained by SNP.
#' 
#' 10. Build pv.anova object, which summarizes the variance explained by the tested SNP, if that amount is significant, the gblup variance components, and how many sig. SNPs remain after fixing top SNP.
#' 
#' 11. Return these results!
#' 
#' 
#' ### coloc.test Function:
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

#' Co-localization test: Fix top peak marker per pQTL and fit GWAS model per gene expression
testCo <- lapply(idx, function(x)
	coloc.test(rsp=x, qtl=sig, eqtl=snp.col[[x]], coloc=pheno[[x]], top=TRUE,
		design=design, dataM=dataM, G=G, wt=wtcen, Z=t(Zfull), fdr=0.05, pos=c(T,T)))
names(testCo) <- idx

#' Summarize results
rstCo <- do.call(rbind, lapply(idx, function(x) data.frame(gene=x,
	pheno=unlist(lapply(rownames(testCo[[x]]), function(y) strsplit(y,"[.]")[[1]][1])),
	testCo[[x]][,c("snp", "p.value", "propVS","sigsnp")])))
rownames(rstCo) <- NULL

sigCo <- rstCo[rstCo$sigsnp == 0,]
sigEff <- rstCo[rstCo$p.value < 0.05 & rstCo$sigsnp > 0,]

table(as.character(sigCo$gene))
rownames(dge$genes)<-gsub("-",".", rownames(dge$genes))
dge$genes[as.character(unique(sigCo$gene)),]

#' Proportion of variance explained by top pQTL SNP
sigComiR <- lapply(as.character(unique(sigCo$gene)), function(x)
	table(as.character(sigCo$pheno)[grep(x, sigCo$gene)]))
names(sigComiR) <- as.character(unique(sigCo$gene))
summary(sigCo$propVS*100)

sigEffmiR <- lapply(as.character(unique(sigEff$gene)), function(x)
	table(as.character(sigEff$pheno)[grep(x, sigEff$gene)]))
names(sigEffmiR) <- as.character(unique(sigEff$gene))
summary(sigEff$propVS*100)

#' Investigate correlation of peak SNP when duplicates exist in miR-pheno combinations:
#' 
lapply(unique(as.character(sigCo$pheno)), function(x) sigCo[sigCo$pheno==x,])

lapply(unique(as.character(sigCo$pheno)), function(x) table(as.character(sigCo[sigCo$pheno==x,"snp"])))

lapply(unique(as.character(sigCo$pheno)), function(x) table(as.character(sigCo[sigCo$pheno==x,"gene"])))

#' ## Save data
save(testCo, sigCo, sigEff, file="../6_sig_coloc.Rdata")