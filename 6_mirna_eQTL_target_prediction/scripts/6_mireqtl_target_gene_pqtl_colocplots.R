#' **Script:** `6_mireqtl_target_gene_pqtl_colocplots.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts`
#' 
#' **Date:**  `12/18/17`
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/pigeqtl/analyses/eQTL/paper/output`
#' 
#' 1. `/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/`
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/5_gblup_gwa_eqtl`
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/6_mirna_eQTL_target_prediction`
#' 
#' **Input File(s):** 
#' 
#' 1. `funct_eqtl.Rdata`
#' 
#' 1. `pQTL_60K.Rdata`
#' 
#' 1. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 1. `10_mrna_mirna_corr_rst.Rdata` `13_target_mrna_coloc_pqtl.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts/6_mireqtl_target_gene_pqtl_colocplots_literate/figure`
#' 
#' **Output File(s):** 
#' 
#' `ssc-miR-874-1.tiff`
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
#' This scripts creates a graphical representation of the co-localization of ssc-miR-874 target genes with pQTL
#' by first plotting the manhattan plot for ssc-miR-874 eQTL and then overlappling the pQTL with co-localized
#' ssc-miR-874 target genes. The targets genes apping within each pQTL peak is highlighted in red.
#' 
#' ## Install libraries
#' 

#' ## Load data
#' 
#' Clear environment
rm(list=ls())

#' Session Information
sessionInfo()

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/6_mirna_eQTL_target_prediction/scripts")

#' Load required R objects
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")

load("../../5_gblup_gwa_eqtl/2_gwa_results.Rdata")
load("../../5_gblup_gwa_eqtl/3_eqtl_summary_tables_maps.Rdata")
load("../10_mrna_mirna_corr_rst.Rdata")
load("../13_target_mrna_coloc_pqtl.Rdata")

ls()
#' ## Analysis
#' 

#' Identify map positions of negatively correlated genes for miR-874:
#' 
#' First, calculate the midpoint of the gene:
annot$pos<-annot$start+(round((annot$end - annot$start)/2))
head(annot)

annot<-annot[annot$miRNA=="ssc-miR-874",]
annot<-annot[annot$cor<0,]
annot<-annot[annot$chr!="AEMK02000682.1",]
annot$chr<-gsub("X", "19", annot$chr)
head(annot)
dim(annot)

#' miRNA being analyzed (manual input)
X <- "ssc-miR-874"

#' Extract the results for the miRNA X
miX <- negcoloc[negcoloc$miRNA == X,]
miX$chr<-as.numeric(as.character(miX$chr))
str(miX)
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

#' Check that the number of targets matches the number of significantly correlated
#' target genes for miRNA X: should be TRUE
length(table(as.character(miX$ID))) == nrow(targets.miX)

#' Add the position of the significantly correlated target genes for miRNA X to the map.full object
map.full <- rbind(map.full, targets.miX[,c("chr", "pos")])
map.full <- rbind(map.full, annot[,c("chr", "pos")])
map.full$chr<-as.numeric(map.full$chr)

#' Use the absmap function to convert the chromosomal map positions to absolute map positions:
absposmap <- absmap(map.full)/1e6

#' Obtain phenotypic QTL markers only for the pQTL with a co-localized target gene
pheno <- unique(as.character(miX$pheno))
sig.pheno <- lapply(pheno, function(x) GWAS$qvalue[GWAS$qvalue[,x]<0.05, x])
names(sig.pheno) <- pheno

#' Check markers associated to a pQTL that are not in the absposmap
lapply(sig.pheno, function(x) names(x)[!names(x) %in% names(absposmap)])
unlist(lapply(sig.pheno, length))

#' Remove markers associated to a pQTL that is not co-localized with the target genes of miRNA X
sig.pheno <- lapply(pheno, function(x) sig.pheno[[x]][map.full[names(sig.pheno[[x]]), "chr"] %in%
	as.character(miX[as.character(miX$pheno) == x,"chr"])])
names(sig.pheno) <- pheno
unlist(lapply(sig.pheno, length))

#' Create list with targets co-localized with pQTL (used to plot pQTL names in manhatan plot)
bychr <- lapply(unique(as.character(miX$chr)), function(x) miX[as.character(miX$chr) == x,])
bychr <- lapply(bychr, function(x) data.frame(chr=rep(unique(x$chr), length(unique(x$pheno))),
	pheno=unique(as.character(x$pheno)),
	do.call(rbind, lapply(unique(as.character(x$pheno)), function(y) data.frame(
		start=min(absposmap[names(sig.pheno[[y]])[map.full[names(sig.pheno[[y]]),"chr"]
			%in% unique(x$chr)]]),
		end=max(absposmap[names(sig.pheno[[y]])[map.full[names(sig.pheno[[y]]),"chr"]
			%in% unique(x$chr)]]),
		snp.qval=min(sig.pheno[[y]][map.full[names(sig.pheno[[y]]),"chr"] %in% unique(x$chr)]))))))

#' Manually specify were the text for each pQTL will be (left, right or top of peak)
# Position next to peak
pos <- c("top", "top", "top", "left", "top", "right")
# Position of text next to peak (pos object in R plots)
Spos <- c("top", "top", "top", "left", "top", "right")

#' Generate the manhatan plot for the ssc-miR-874 eQTL
qvals <- rst.gwa[rst.gwa$miRNA==X, "gwa.qval"]
names(qvals) <- rst.gwa[rst.gwa$miRNA==X, "SNPid"]

#' Position of miRNA X eQTL
eqtlX <- map.full[names(qvals)[qvals < 0.05],]
pos_eqtlX <- unlist(lapply(unique(eqtlX$chr), function(x) paste("SSC", x, ":",
	ifelse(nrow(eqtlX[eqtlX$chr == x,]) == 1,
		paste(round(eqtlX[eqtlX$chr == x, "pos"]/1e6, 3), " Mb", sep=""),
		paste(round(min(eqtlX[eqtlX$chr == x, "pos"])/1e6, 3),
	"-", end=round(max(eqtlX[eqtlX$chr == x, "pos"])/1e6, 3), " Mb", sep="")), sep="")))

#' ## Visualize
#' 
#' Plot eQTL and target genes co-localizing with pQTL
#+ ssc-miR-874, fig.align='center', fig.width=20, fig.height=9, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
par(oma=c(2,2,2,2))
manhpt(nm=X,abspos=absposmap,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,20))
points(absposmap[as.character(rownames(annot))], -log10(annot$qvalue), pch=20, col="grey60", cex=1)
lapply(pheno, function(x) points(absposmap[names(sig.pheno[[x]])],
	-log10(sig.pheno[[x]]), col="chocolate1", pch=1))
# Highlight targets co-localized with pQTL
points(absposmap[as.character(targets.miX$geneID)], -log10(targets.miX$cor.qval), pch=20, col="black", cex=1.5)

# Highlight eQTL peak
text(x=0, y=19, labels="eQTL Chromosome: Position", font=2, pos=4, cex=2)
ifelse(length(pos_eqtlX) > 1, lapply(1:length(pos_eqtlX),
	function(i) text(x=0, y=19 - (i * 1.5), labels=pos_eqtlX[i], pos=4, cex=2)),
	text(x=0, y=19 - 1.5, labels=pos_eqtlX, pos=4, cex=2))
text(absposmap[X], max(-log10(qvals)) + 1, labels="eQTL", font=2, cex=2)
# Add chromosome to pQTL peaks
lapply(1:length(bychr), function(z)
	text(x=ifelse(pos[z] == "left", min(bychr[[z]]$start) - 5, ifelse(pos[z] == "right",
			max(bychr[[z]]$end) + 10, min(bychr[[z]]$start))),
		y=ifelse(pos[z]=="top", max(-log10(bychr[[z]]$snp.qval)) + (nrow(bychr[[z]]) * 2.5),
			max(-log10(bychr[[z]]$snp.qval))),
	labels=paste("SSC", unique(bychr[[z]]$chr), ":", sep=""),
	pos=ifelse(Spos[z] == "right", 4, 2), font=2, cex=2))
# Add phenotype names to pQTL peaks
lapply(1:length(bychr), function(z) lapply(1:nrow(bychr[[z]]), function(i)
	text(x=ifelse(pos[z] == "left", min(bychr[[z]]$start) - 5, ifelse(pos[z] == "right",
			max(bychr[[z]]$end) + 10, min(bychr[[z]]$start))),
		y=ifelse(pos[z]=="top", (max(-log10(bychr[[z]]$snp.qval)) + (nrow(bychr[[z]]) * 2.5)) - (1 * i),
			max(-log10(bychr[[z]]$snp.qval)) - (1 * i)),
	labels=bychr[[z]][i,"pheno"], pos=ifelse(Spos[z] == "right", 4, 2), cex=2)))

add_legend("bottomright",legend=c("miR-874 SNP Assoc.", "pQTL SNP Assoc.", "Target Genes", "Co-localized Targets"),
	pch=c(1,1), bty="o",ncol=2, cex=1.2)
	add_legend("bottomright",legend=c("miR-874 SNP Assoc.", "pQTL SNP Assoc.", "Target Genes", "Co-localized Targets"),
	pch=19, 
	col=c("deepskyblue", "chocolate1", "grey60", "black"), bty="n",ncol=2, cex=1.2)
