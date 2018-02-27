#' **Script:** `4_mRNA_miRNA_corr.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  12/13/17
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/`
#' 
#' 2. `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11`, `/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/`
#' 
#' 3. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' 4. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`
#' 
#' 5. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/`
#' 
#' **Input File(s):** 
#' 
#' 1. `voom.Rdata`
#' 
#' 2. `classified_peaks.Rdata` (need to ask DV, only for co-localizing )
#' 
#' 3. `G_Z_matrix_Ss11.Rdata`, `MSUPRP_gpData_Ss11.Rdata`
#' 
#' 4. `eQTL60K.Rdata`
#' 
#' 5. `5_filtered_targets_exp_rst.Rdata`
#' 
#' 6. `3_msuprp_mirna_gpdata.Rdata`
#' 
#' 7. `4_normalized_dge_voom.Rdata`
#' 
#' 8. `5_Z_G_miRNA.Rdata`
#' 
#' 9. `1_hsa_mirna_names.txt`, `1_mirna_names.txt`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction`
#' 
#' **Output File(s):** 
#' 
#' 1. `1.5_ssc_hsa_hom_mir.txt`
#' 
#' 1. `8_DAVID_cor_target_genes.txt`
#' 
#' 1. `9_DAVID_neg_cor_target_genes.txt`
#' 
#' 1. `10_mrna_mirna_corr_rst.Rdata`
#' 
#' 1. `11_mrna_mirna_corr_char.Rdata`
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
#' The objective of this script is to correlate the expression of the miR-eQTL miRNAs with their putative target mRNAs. 
#' The mRNA and miRNA expression data will separately be adjusted for the same fixed (sex, selection criteria) and random (population stratification) effects as in the eQTL analysis, accounting for the mean-variance relationship, as well. 
#' The adjusted expression data will then be correlated between the miRNA and its putative target mRNAs, utilizing the residual expression data (see Ponsuksili et al. 2013 [BMC Genomics] for details). 
#' Correlations found to be significant will then be characterized for their overlap with mRNA-eQTLs, and pQTLs in this dataset.
#' 
#' ## Install libraries

rm(list=ls())

library(regress)
library(gwaR)
library(limma)
library(edgeR)
library(parallel)
library(qvalue)
library(corrplot)

#' ## Load data

load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/voom/voom.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/G_matrix_Ss11/G_Z_matrix_Ss11.Rdata")
ls()

#' Rename R object to differentiate between mRNA and microRNA
Mdge <- dge
Mv <- v
Mwcen <- wcen
MG <- G

#' Remove R objects that will not be used in this analysis
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168")))
ls()

#' Load microRNA targets
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/5_filtered_targets_exp_rst.Rdata")
ls()

#' Summary of target genes per microRNA
do.call(rbind, lapply(targets.exp.sum, function(x) x[[1]]))

#' mRNA annotation
annot <- Mdge$genes
head(annot)

#' List of mRNA target genes per microRNA
#' 
#' Some genes have more than one transcript expressed per target gene, meaning the number of potential 'target' transcripts will increase.
targets.mrna <- lapply(targets.exp, function(x) 
	rownames(annot)[as.character(annot$genes) %in% x$external_gene_name])
names(targets.mrna) <- names(targets.exp)
head(targets.mrna[[1]])
str(targets.mrna)

#' Load microRNA expression data
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()

#' Rename R object to differentiate microRNA data from mRNA data
dge.mi <- dge
v.mi <- v
wcen.mi <- wtcen
G.mi <- G

#' Retain only the objects needed for the analysis
rm(list=setdiff(ls(), c("Mdge","Mv","Mwcen","MG","MSUPRP168",
	"dge.mi","v.mi","wcen.mi","G.mi","MSUPRP_miRNA",
	"targets.mrna","targets.exp")))
ls()

#' ## Analysis
#' 
#' miRNA human-pig homologs
hsa.mirna <- as.character(read.table("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/1_hsa_mirna_names.txt")[,1])
hsa.mirna <- c(hsa.mirna[1], hsa.mirna)
ssc.mirna <- apply(read.table("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/1_mirna_names.txt"), 1, function(x) paste("ssc-", x, sep=""))
ssc.hsa.hom.mir <- data.frame(ssc=ssc.mirna[1:15], hsa=hsa.mirna)
ssc.hsa.hom.mir

#' Save miRNA human-pig homologes table
write.table(ssc.hsa.hom.mir, "/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/1.5_ssc_hsa_hom_mir.txt", col.names=TRUE, quote=FALSE, row.names=FALSE, sep="\t")


#' Reorder targets.mrna to match order in ssc.hsa.hom.mir
targets.mrna <- targets.mrna[as.character(ssc.hsa.hom.mir[-1,"hsa"])]

#' Rename targets.mrna R object to match ssc miRNA id (currently has hsa id names)
x <- unlist(sapply(names(targets.mrna), 
	function(x) as.character(ssc.hsa.hom.mir[as.character(ssc.hsa.hom.mir$hsa) %in% x,"ssc"])))
x

#' Add ssc names to targets.mrna and duplicate let-7-5p/98-5p target genes to account for the two ssc miRNA 
targets.mrna <- c(targets.mrna[1], targets.mrna)
names(targets.mrna) <- x
str(targets.mrna)

#' Function to extract residual values from regress function after correcting for sex, growth group, population structure and mean variance relationship
fit.vals <- function(rsp, data, design, G, vdata = NULL, wt = NULL, ...){
	x <- gblup(rsp=rsp, data=dataM, design=design, G=G, vdata=vdata, wt=wt, ...)
	ehat <- x$ehat
	rownames(ehat) <- names(x$model$y)
	colnames(ehat) <- rsp
	return(ehat)
}

#' Estimate the **miRNA** expression (log-cpm) after correcting for 
#' sex, selection criteria, population structure and mean variance relationship
idx <- as.character(ssc.hsa.hom.mir$ssc)
length(idx)

str(MSUPRP_miRNA$covar[,c("sex","growth_group")])
X <- MSUPRP_miRNA$covar[,c("sex","growth_group")]
rownames(X) <- MSUPRP_miRNA$covar$id

dataM <- cbind(t(v.mi$E[idx,]), X[colnames(v.mi$E),])

design <- c(~ sex + growth_group)

fit.mi <- do.call(cbind, lapply(idx, function(x) fit.vals(rsp=x, data=dataM, design=design, G=G.mi, 
	vdata = NULL, wt = wcen.mi, pos=c(T,T))))
dim(fit.mi)

#' Estimate the **mRNA** expression (log-cpm) after correcting for 
#' sex, selection criteria, population structure and mean variance relationship
# Target mRNA
idx <- unique(unlist(targets.mrna))
length(idx)

str(MSUPRP168$covar[,c("sex","selcrit")])
X <- MSUPRP168$covar[,c("sex","selcrit")]
rownames(X) <- MSUPRP168$covar$id

dataM <- cbind(t(Mv$E[idx,]), X[colnames(Mv$E),])

design <- c(~ sex + selcrit)

Mfit <- do.call(cbind, mclapply(idx, function(x) fit.vals(rsp=x, data=dataM, design=design, G=MG, 
	vdata = NULL, wt = Mwcen, pos=c(T,T)), mc.cores=10))
dim(Mfit)


#' Correlate adjusted miRNA expression with adjusted mRNA expression
# Reduced fitted matrices for animals with both mRNA and miRNA fitted values
anim <- rownames(Mfit)[rownames(Mfit) %in% rownames(fit.mi)]
# miRNA
rfit.mi <- fit.mi[anim,]
dim(rfit.mi)
# mRNA
rMfit <- Mfit[anim,]
dim(rMfit)


#' Function to perform correlation analysis of mRNA and miRNA espression and export a summary of correlation analysis
cor.exp <- function(target, mRNA, miRNA, ...){
	x <- cor.test(mRNA, miRNA, ...)
	rst <- data.frame(cor=x$estimate, 
		z=x$statistic, 
		pvalue=x$p.value)
	rownames(rst) <- target
	return(rst)
}

#' ---
#' Run corrlation analysis
#' 
#' Correlation between miRNA and their target mRNA
rst.corR <- lapply(names(targets.mrna), function(x) do.call(rbind,mclapply(targets.mrna[[x]], 
		function(y) cor.exp(target=y, mRNA=rMfit[,y], 
		miRNA=rfit.mi[,x], alternative="two.sided", method="kendall"), mc.cores=10)))
names(rst.corR) <- names(targets.mrna)

#' False discovering rate multiple test correction
for (i in names(rst.corR)){
	rst.corR[[i]] <- data.frame(rst.corR[[i]], qvalue=qvalue(rst.corR[[i]]$pvalue)$qvalue, pi0=qvalue(rst.corR[[i]]$pvalue)$pi0)
}

#' Summary of significantly correlated targets at FDR < 0.05:
thres <- 0.05
summary.sigR <- do.call(rbind, lapply(rst.corR, function(x) data.frame(tot=nrow(x), 
	sig=sum(x$qvalue < thres), prop.sig=sum(x$qvalue < thres)/nrow(x) * 100, 
	neg=sum(x$cor < 0), prop.neg=sum(x$cor < 0)/nrow(x) *100,
	sig.neg=sum(x$cor < 0 & x$qvalue < thres),
	prop.sig.neg=sum(x$cor < 0 & x$qvalue < thres)/nrow(x) * 100,
	pi0=unique(x$pi0))))
summary.sigR

#' Observe the histogram of p-values for the correlation analysis:
tst<-lapply(rst.corR, function(x) hist(x$pvalue))
lapply(tst, function(x) x$counts)

#' Recognize that the correlations may not be incredibly strong: many miRNAs can target one gene, 
#' and the effect of the miRNA on a given target varies. 
#' 
#' List of mRNA XLOC IDs that are significantly correlated per miRNA
sig.mrnaR <- lapply(rst.corR, function(x) rownames(x)[x$qvalue < thres])

head(as.character(rownames(Mdge$genes)))
Mdge$genes$geneXID<-as.character(rownames(Mdge$genes))

rst.corR[[2]][sig.mrnaR[[1]],]
rst.corR[[1]][sig.mrnaR[[1]],]

#' Extract the correlations for those significant miRNA-mRNA pairs
sig.mrnacorR <- lapply(rst.corR, function(x) x[x$qvalue < thres, "cor"])

#' Extract the names of the significant target mRNAs for use with DAVID:
sig.mrnaR.names <- do.call(rbind, lapply(names(sig.mrnaR), function(x) data.frame(miRNA=rep(x, length(sig.mrnaR[[x]])), 
    Mdge$genes[sig.mrnaR[[x]], c("geneXID", "genes")])))
rownames(sig.mrnaR.names) <- NULL
dim(sig.mrnaR.names)
head(sig.mrnaR.names)

#' These sums should match the summary.sigR "sig.neg" column
unlist(lapply(sig.mrnacorR, function(x) sum(x < 0)))

#' Extract the significant negatively-correlated miRNA-mRNA pairs:
sig.neg.mrnaR <- lapply(rst.corR, function(x) rownames(x)[x$qvalue < thres & x$cor<0])
#' Extract the names of those genes for use in DAVID analysis:
sig.neg.mrnaR.names <- do.call(rbind, lapply(names(sig.neg.mrnaR), function(x) data.frame(miRNA=rep(x, length(sig.neg.mrnaR[[x]])), 
    Mdge$genes[sig.neg.mrnaR[[x]], c("geneXID", "genes")])))
rownames(sig.neg.mrnaR.names) <- NULL
dim(sig.neg.mrnaR.names)
head(sig.neg.mrnaR.names)

#' Names of unique miRNA target genes per miRNA:
unique(as.character(sig.neg.mrnaR.names$miRNA))
# lapply(unique(as.character(sig.neg.mrnaR.names$miRNA)), function(x) unique(as.character(sig.neg.mrnaR.names[sig.neg.mrnaR.names$miRNA==x, "genes"])))

#' ----
#' 
#' Code received from DV 6/29/17
#' 
#' Correlation among miRNA expression 
cor.mirna <- cor(t(v.mi$E[names(targets.mrna),]))
#+ corr_miRNA_exp, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(cor.mirna, method="number", mar=c(2,2,2,2), number.cex=0.7, tl.cex=0.8, tl.col="black", tl.srt=45, title="Correlation of microRNA Expression")

#' Xloc IDs of common targets per miRNA
com.targets <- lapply(names(targets.mrna), function(x) lapply(targets.mrna, 
    function(y) targets.mrna[[x]][targets.mrna[[x]] %in% y]))
names(com.targets)<-names(targets.mrna)

#' Number of common targets per miRNA
comM <- do.call(cbind, lapply(com.targets, 
    function(x) unlist(lapply(x, function(y) length(y)))))
colnames(comM) <- rownames(comM)

#' Proportion of common targets per miRNA
comMperc <- do.call(cbind, lapply(names(targets.mrna), function(x) unlist(lapply(targets.mrna, 
    function(y) sum(targets.mrna[[x]] %in% y)/length(targets.mrna[[x]])))))
colnames(comMperc) <- rownames(comMperc)
#+ miRNA_targets_common, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(comMperc, method="number", mar=c(2,2,2,2), number.cex=0.7, tl.cex=0.8, tl.col="black", tl.srt=45, title="Proportion of microRNA Targets in Common")

#' ---
#' Number of times a single mRNA target is found significantly correlated to a miRNA
idx <- rownames(summary.sigR)[summary.sigR$sig > 0]
tail(sort(table(unlist(lapply(idx, function(x) rownames(rst.corR[[x]])[rst.corR[[x]]$qvalue < thres])))),n=50)

#' Identify target mRNAs that have significant correlation with multiple miRNAs:
com.sig <- lapply(idx, function(x) lapply(com.targets[[x]][idx], 
    function(y) y[y %in% rownames(rst.corR[[x]])[rst.corR[[x]]$qvalue < thres]]))
names(com.sig) <- idx
com.sig <- lapply(com.sig, function(x) lapply(names(x), 
    function(y) x[[y]][x[[y]] %in% rownames(rst.corR[[y]])[rst.corR[[y]]$qvalue < thres]]))
for (i in idx){
    names(com.sig[[i]]) <- idx 
    com.sig[[i]] <- com.sig[[i]][-grep(i, names(com.sig[[i]]))]
}

#' Check that this function worked correctly:
rst.corR[["ssc-miR-874"]][grep("XLOC_009194", rownames(rst.corR[["ssc-miR-874"]])),]
rst.corR[["ssc-miR-6782-3p"]][grep("XLOC_009194", rownames(rst.corR[["ssc-miR-6782-3p"]])),]

#' List of mRNA significantly correlated to miRNA and also have an associated eQTL
# load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/classified_peaks.Rdata")
# sig.mrna.eqtlR<-lapply(sig.mrnaR, function(x) regul[regul$gene %in% x,])
# sig.mrna.eqtlR

#' List of correlation results only for mRNA corrlelated to miRNA and have an associated eQTL
# lapply(1:length(sig.mrnaR), function(x) rst.corR[[x]][sig.mrnaR[[x]][sig.mrnaR[[x]] %in% regul$gene],])

#' ## Visualize
#' 
#' Plot the correlation vs the significance of the miRNA with its putative targets:
#+ miR-874_cor, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
# plot(rst.corR$"ssc-miR-874"$cor, 
# 	-log10(rst.corR$"ssc-miR-874"$qvalue),
# 	pch=20, 
# 	xlab="miRNA-mRNA Correlation",
# 	ylab="Significance (-log10(qvalue))", 
# 	main="ssc-miR-874", 
# 	cex.main=1.7,
# 	cex.lab=1.3,
# 	cex.axis=1.3)
# points(rst.corR$"ssc-miR-874"$cor[-log10(rst.corR$"ssc-miR-874"$qvalue)>(-log10(thres))], 
# 	-log10(rst.corR$"ssc-miR-874"$qvalue)[-log10(rst.corR$"ssc-miR-874"$qvalue)>(-log10(thres))], 
# 	pch=20,
# 	col="red")
# abline(a=-log10(thres), b=0, lty=5)
# text(rst.corR$"ssc-miR-874"["XLOC_011914","cor"], 
# 	-log10(rst.corR$"ssc-miR-874"["XLOC_011914","qvalue"]), 
# 	labels=paste(sig.mrna.eqtlR$"ssc-miR-874"$gene.name,
# 		format(round(rst.corR$"ssc-miR-874"["XLOC_011914","cor"], 3))),
# 	pos=4,
# 	col="blue")
#' 

#' Code to run with DV's "runR" hack:
# 4_mRNA_miRNA_corr.R nodes=1:ppn=11,walltime=00:15:00,mem=20G

#' ## Save Data
write.table(sig.mrnaR.names, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/8_DAVID_cor_target_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(sig.neg.mrnaR.names, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/9_DAVID_neg_cor_target_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
save(rst.corR, summary.sigR, sig.mrnaR, sig.neg.mrnaR, sig.mrnaR.names, sig.neg.mrnaR.names, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/10_mrna_mirna_corr_rst.Rdata")
save(cor.mirna, com.targets, comM, comMperc, com.sig, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/11_mrna_mirna_corr_char.Rdata")
save(rMfit, rfit.mi, file="/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/12_mrna_mirna_resid_exp.Rdata")