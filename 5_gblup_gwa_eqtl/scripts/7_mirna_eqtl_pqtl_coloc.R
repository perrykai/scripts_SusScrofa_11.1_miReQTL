#' **Script:** `7_mirna_eqtl_pqtl_coloc.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  12/18/17
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
#' 1. `pQTL_60K.Rdata`, `inrange_function.Rdata`
#' 
#' 2. `MSUPRP_gpData_Ss11.Rdata`
#' 
#' 3. `6_mirna_precursor_annot_ssc11.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`
#' 
#' 4. `1_gblup_results_summary.Rdata`, `2_gwa_results.Rdata`, `4_miRNA_eQTL_local_distal_regulators.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `9_mireqtl_pqtl_coloc_peaks.Rdata`
#' 
#' 2. `10_mireqtl_pqtl_coloc_peaks.txt`
#' 
#' 3. `11_pQTL_summary.txt`
#' 
#' 4. `12_summary_pheno_174.txt`
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
#' The objective of this script is to investigate if axny miRNA eQTL co-localize with pQTL.
#' 
#' ## Install libraries
setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

rm(list=ls())

library(methods)
library(limma)
library(edgeR)

#' ## Load data
load("/mnt/research/ernstc_lab/RNAseq_ASE/QTL/pQTL60k/pQTL_60k.Rdata")
load("/mnt/research/ernstc_lab/RNAseq_ASE/SNP60K_Ss11/SNP60_Ss11_Map/MSUPRP_gpData_Ss11.Rdata")
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/corrected-Z/inrange_function.Rdata")

#' Load DV functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
load("../../5_gblup_gwa_eqtl/1_gblup_results_summary.Rdata")
load("../../5_gblup_gwa_eqtl/2_gwa_results.Rdata")
load("../../5_gblup_gwa_eqtl/4_miRNA_eQTL_local_distal_regulators.Rdata")

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

#' Marker Map
map<-MSUPRP$map
str(map)
map$chr<-gsub("X", "19", map$chr)
map$chr<-gsub("Y", "20", map$chr)
str(map)
map$chr<-as.numeric(map$chr)

#' Retain marker information for all pQTL
sig <- apply(qval, 2, function(x) x[x < 0.05])
sig <- lapply(names(sig), function(x) data.frame(map[names(sig[[x]]),],
    std.eff=sdEff[names(sig[[x]]),x], pvalue=pval[names(sig[[x]]),x], qvalue=sig[[x]]))
names(sig) <- colnames(qval)
sig <- sig[unlist(lapply(sig, nrow)) > 0]

#' 25 of the 67 traits have pQTL:
length(sig)
names(sig)

#' If a pQTL contains more than one peak split each peak
idx <- lapply(sig, function(x) as.numeric(names(table(x$chr))))
idx <- idx[unlist(lapply(idx, function(x) length(x) > 1))]

mp <- lapply(1:length(idx), function(x) lapply(1:length(idx[[x]]), function(y) sig[[names(idx)[x]]][sig[[names(idx)[x]]]$chr == idx[[x]][y],]))
names(mp) <- names(idx)

for (i in 1:length(mp)){
names(mp[[i]]) <- paste(names(mp)[[i]], 1:length(mp[[i]]), sep=".")
}

qtl <- sig[!names(sig) %in% names(mp)]
for (i in 1:length(mp)){
qtl <- c(qtl, mp[[i]])
}

# Number of associated markers per QTL peak
length(qtl)
unlist(lapply(qtl, nrow))

#' Co-localization eQTL with pQTL
# Check all eQTL peaks within a pQTL
sig <- cbind(regul[,c(7:14)], regul[,c(1:6)])
head(sig)

#' qtl = pQTL peaks
#' map = sig refers to range of miR-eQTL peaks
#' 
#' map = rqtl refers to range of pQTL peaks
#' The inrange function is then used to identify if the pQTL peak lies within a miR-eQTL, 
#' or overlaps the miR-eQTL peak on the left or the right
#' 
#' 
#' Identify miRNA eQTL where the position of the peak miR-eQTL SNP lies between the beginning and end of the pQTL peak:
eqtl.pqtl <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="pos.snp", range=NULL))
unlist(lapply(eqtl.pqtl, nrow))[unlist(lapply(eqtl.pqtl, nrow)) > 0]

eqtl.pqtl[names(unlist(lapply(eqtl.pqtl, nrow))[unlist(lapply(eqtl.pqtl, nrow)) > 0])]

#' rqtl = range of pQTL peaks:
rqtl <- do.call(rbind,lapply(qtl, function(x) data.frame(chr=unique(x$chr),
        min=min(x$pos), max=max(x$pos))))

#' ### pQTL peaks lying completely within a miR-eQTL peak:
#' 
#' (start of pQTL peak is after start of miR-eQTL peak, and 
#' end of pQTL peak is before end of miR-eQTL peak)
PinE <- lapply(1:nrow(sig), function(x) inrange(chr=sig[x, "chr.snp"], start=sig[x, "min.pos"],
        end=sig[x, "max.pos"], map=rqtl, single=NULL, range=c(start="min", end="max")))
names(PinE) <- 1:nrow(sig)
unlist(lapply(PinE, nrow))[unlist(lapply(PinE, nrow)) > 0]

PinE <- do.call(rbind, lapply(names(PinE), function(x)
        data.frame(pqtl=rownames(PinE[[x]]), PinE[[x]], eqtl=rep(x, nrow(PinE[[x]])))))
rownames(PinE) <- NULL

PinE <- lapply(names(qtl), function(x) sig[as.character(PinE$eqtl[grep(x, PinE$pqtl)]),])
names(PinE) <- names(qtl)

unlist(lapply(PinE, nrow))[unlist(lapply(PinE, nrow)) > 0]
PinE[names(unlist(lapply(PinE, nrow))[unlist(lapply(PinE, nrow)) > 0])]

#' ### pQTL peaks overlapping left side of miR-eQTL:
#' 
#' (minimum position of miR-eQTL lies after start of pQTL peak, and before end of pQTL peak)
left <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="min.pos", range=NULL))
#' Number of pQTL overlapping left side of miR-eQTL:
sum(lapply(left, nrow) >0)

unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0]
names(unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0])

#' Notice it's the same 6 miRNAs for all of these pQTL:
left[names(unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0])]

#' ### pQTL peaks overlapping right side of miR-eQTL:
#' 
#' (maximum position of miR-eQTL lies after start of pQTL peak, and before end of pQTL peak)
right <- lapply(qtl, function(x) inrange(chr=unique(x$chr), start=min(x$pos),
        end=max(x$pos), map=sig, single="max.pos", range=NULL))
unlist(lapply(right, nrow))[unlist(lapply(right, nrow)) > 0]
names(unlist(lapply(right, nrow))[unlist(lapply(right, nrow)) > 0])

right[names(unlist(lapply(right, nrow))[unlist(lapply(right, nrow)) > 0])]

#' Notice mostly the same pQTL traits between left and right;
names(unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0])
names(unlist(lapply(right, nrow))[unlist(lapply(right, nrow)) > 0])

#' See if all the miR-eQTL peaks end before the pQTL peak ends:
#' 
#' If the minimum of the pQTL peak is less than the minimum of the miR-eQTL peak, 
#' and the maximum of the pQTL peak is greater than the maximum of the miR-eQTL peak, 
#' then the miR-eQTL peak lies completely within the pQTL peak.
lrmir<-sig[sig$miRNA %in% left$juiciness$miRNA,]
lrmir

lrpqtl<-rqtl[names(unlist(lapply(left, nrow))[unlist(lapply(left, nrow)) > 0]),]

for(i in 1:nrow(lrmir)){
    tst<-data.frame(lb=lrmir[i,"min.pos"]>lrpqtl$min,
        ub=lrmir[i,"max.pos"]<lrpqtl$max)
}
rownames(tst)<-rownames(lrpqtl)
tst

#' Notice for juiciness there is a "false"; investigate this:
for(j in 1:nrow(lrpqtl)){
for(i in 1:nrow(lrmir)){
    print(c(rownames(lrpqtl)[j], lrmir$miRNA[i], lrpqtl[j,"min"]<lrmir[i,"min.pos"], lrpqtl[j,"max"]>lrmir[i,"max.pos"]))
}}

#' For juiciness, the MARC0093624 SNP gets involved, in combination with single-SNP miR-eQTL peaks:
#' 
#' ### Merge all pQTL co-localizing with miR-eQTL
coloc <- lapply(names(eqtl.pqtl), function(x) rbind(eqtl.pqtl[[x]],
        PinE[[x]][!rownames(PinE[[x]]) %in% rownames(eqtl.pqtl[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]],
        left[[x]][!rownames(left[[x]]) %in% rownames(coloc[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

coloc <- lapply(names(coloc), function(x) rbind(coloc[[x]],
        right[[x]][!rownames(right[[x]]) %in% rownames(coloc[[x]]),]))
names(coloc) <- names(eqtl.pqtl)

# Matrix of pQTL and colocalized miR-eQTL
qtlM <- do.call(rbind, lapply(names(qtl), function(x)
        data.frame(chr=unique(qtl[[x]]$chr), lowb=min(qtl[[x]]$pos), upperb=max(qtl[[x]]$pos),
                eQTL=nrow(coloc[[x]]))))
rownames(qtlM) <- names(qtl)

# Number of pQTL with co-localized miR-eQTL
sum(qtlM$eQTL > 0)

# Summary of colocalized QTL:
coloc <- coloc[unlist(lapply(coloc, nrow)) > 0]
unlist(lapply(coloc, nrow))


#' Create table with pQTL information
sum.qtl <- do.call(rbind, lapply(qtl, function(x) data.frame(chr=unique(x$chr), start=min(x$pos), end=max(x$pos),
        SNP=rownames(x)[min(x$pvalue) == x$pvalue][1], x[min(x$pvalue) == x$pvalue,][1,
        c("pos", "std.eff", "pvalue", "qvalue")], nSNP=nrow(x))))
sum.qtl <- sum.qtl[order(sum.qtl$chr),]

# Add heritability of each pQTL and number of co-localized eQTL
idx <- unlist(lapply(rownames(sum.qtl), function(x) strsplit(x, "[.]")[[1]][1]))
sum.qtl <- data.frame(sum.qtl, h2=GBLUP[idx,"h2"], pval.h2=GBLUP[idx,"pvalue"],
        qval.h2=GBLUP[idx,"qvalue"], eQTL=unlist(lapply(rownames(sum.qtl),
        function(x) ifelse(x %in% names(coloc), nrow(coloc[[x]]), 0))))
dim(sum.qtl)

#' ---
#' Order colocalized eQTL per chromosome
coloc <- coloc[rownames(sum.qtl)]
coloc <- coloc[!is.na(names(coloc))]

# Number of eQTL co-localized with pQTL
length(unlist(lapply(coloc,rownames)))

# Number of unique eQTL co-localized with pQTL
length(unique(unlist(lapply(coloc,rownames))))

#' Add colocalized eQTL-pQTL infomation to eQTL table (regul R object)
lst <- do.call(rbind,lapply(names(coloc), function(x)
        data.frame(eqtl=rownames(coloc[[x]]), pqtl=strsplit(x,"[.]")[[1]][1])))

regul <- cbind(regul,
        qval.snp=unlist(lapply(1:nrow(regul),
               function(x) rst.gwa[(rst.gwa$miRNA==regul$miRNA[x] & rst.gwa$SNPid==regul$SNP[x]),"gwa.qval"])),
        pval.snp=unlist(lapply(1:nrow(regul),
               function(x) rst.gwa[(rst.gwa$miRNA==regul$miRNA[x] & rst.gwa$SNPid==regul$SNP[x]),"gwa.pval"])),
        summary.rst.gblup[regul$miRNA,c("h2","lrtpvalue","qvalue")],
        colocalized.pqtl=unlist(lapply(rownames(regul),
                function(x) paste(as.character(lst[lst$eqtl == x, "pqtl"]), collapse=', '))))

#' Reorder columns in regul data frame
regul <- regul[,c("chr.snp", "SNP", "pos.snp", "pval.snp", "qval.snp", "min.pos", "max.pos", "range.peak", "num.snp",
        "miRNA", "chr.miR", "start.miR", "end.miR", "range.miR", "h2", "lrtpvalue", "qvalue", "regulator", "colocalized.pqtl")]

head(regul)

#' Order regul matrix by peak chromosome and position
tmp <- lapply(1:18, function(x) regul[regul$chr.snp == x,])
regul <- do.call(rbind, lapply(tmp, function(x) x[order(x$pos.snp),]))

head(regul)

#' ---
#' Estimate the absolute position of markers and gene expressions
#+ warning=FALSE
posann <- (annotation$end+annotation$start)/2e6
names(posann) <- as.character(annotation$miRNA)
mapZ <- map[,1:2]
mapZ$pos <- mapZ$pos/1e6
mapt <- data.frame(chr=gsub("chr","",annotation$chr0),pos=posann)
rownames(mapt)<-make.names(annotation$miRNA,unique=T)
mapt <- mapt[!is.na(as.numeric(as.character(mapt$chr))), ]

map <- rbind(mapZ, mapt)
map$chr <- as.numeric(map$chr)
map <- map[order(map$chr), ]

#' Absolute postions of genes and markers
abspos <- absmap(map)
head(abspos)
tail(abspos)

#' Summary statistics for each phenotype and selected 174 animals
pheno <- MSUPRP$pheno[rownames(MSUPRP_miRNA$pheno[,,1]),,1]
sumstat <- data.frame(N=apply(pheno,2, function(x) length(x) - sum(is.na(x))),
        Mean=apply(pheno,2, mean, na.rm=T), SD=apply(pheno,2, sd, na.rm=T))

#' ## Visualize
#' 
#' #### Manhattan Plots: miRNA eQTL colocalized with pQTL
idx <- lapply(coloc, function(x) as.character(x$miRNA))

# miR-eQTL qvalues
qvalE <- lapply(idx, function(x) data.frame(rst.gwa[(rst.gwa$miRNA%in%x), c("miRNA","SNPid","gwa.qval")]))
traits<-names(qvalE)

#' Check that the SNPid maintained the same order throughout the qvalE data.frame:
for(i in 1:length(qvalE)){
    print(sum(rep(rownames(MSUPRP_miRNA$map), length(unique(qvalE[[i]]$miRNA))) != qvalE[[i]]$SNPid))
}

qvalE<- lapply(names(qvalE), function(x) as.data.frame(split(qvalE[[x]]$"gwa.qval", qvalE[[x]]$"miRNA"), col.names=idx[[x]]))
names(qvalE)<-traits

for(i in 1:length(qvalE)){
    rownames(qvalE[[i]])<-rownames(MSUPRP_miRNA$map)
}

# pQTL qvalues
nms <- unique(unlist(lapply(names(coloc), function(x) strsplit(x, "[.]")[[1]][1])))
qvalP <- qval[,nms]

#' Growth Phenotypes
idx <- c("last_lum")

#+ growth, fig.align='center', echo=FALSE
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 10))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

#' Carcass composition phenotypes
idx <- c("ph_24h","num_ribs")
#+ carcass, fig.align='center', echo=FALSE
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 10))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

#' Meat quality phenotypes
idx <- c("WBS","tenderness","overtend", "juiciness")
#+ quality, fig.align='center', echo=FALSE
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 6))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        points(abspos[rownames(qvalP)], -log10(qvalP[,i]), pch=20,
                col=c("deepskyblue","blue")[(mps%%2)+1])
        abline(h=-log10(0.05),col="red")
}

idx<-c("cook_yield","driploss","protein")
#+ quality_2, fig.align='center', echo=FALSE
for(i in idx){
        manhpt(nm=i, abspos=abspos, rst=qvalP, map=map[rownames(qvalP),], fdr=0.05,
                pvalues=FALSE, pQTL=TRUE, ylim=c(0, 20))
        x <- do.call(cbind,lapply(grep(i, names(qvalE)), function(x) qvalE[[x]]))
        lapply(colnames(x), function(j) points(abspos[rownames(x)[x[,j] < 0.05]],
                -log10(x[,j][x[,j] < 0.05]), col="orange", pch=1))
        mps <- as.numeric(as.factor(map[rownames(qvalP),]$chr))
        abline(h=-log10(0.05),col="red")
}

#' ## Save data
#' 
#' Save the colocalized miRNA eQTL pQTL peaks:
save(regul, file="../9_mireqtl_pqtl_coloc_peaks.Rdata")
#' Write table of colocalized miRNA eQTL pQTL peaks:
write.table(regul, quote=F, col.names=T, row.names=F, sep="\t", file="../10_mireqtl_pqtl_coloc_peaks.txt")
#' Save pQTL table to text file
write.table(sum.qtl, quote=F, row.names=T, col.names=T, sep="\t",file="../11_pQTL_summary.txt")
#' Save summary statistics for each phenotype and selected 174 animals
write.table(sumstat, quote=F, col.names=T, row.names=T, sep="\t", file="../12_summary_pheno_174.txt")