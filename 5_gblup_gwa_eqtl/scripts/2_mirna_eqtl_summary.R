#' **Script:** `2_mirna_eqtl_summary.R`
#'
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#'
#' **Date:**  12/6/17
#'
#' **Input File Directory:**
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#'
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`
#'
#' **Input File(s):**
#'
#' 1. `1_gblup_results_summary.Rdata`
#'
#' 2. `2_gwa_results.Rdata`
#' 
#' 3. `3_msuprp_mirna_gpdata.Rdata`
#'
#' 4. `4_normalized_dge_voom.Rdata`
#' 
#' 5. `5_Z_G_miRNA.Rdata`
#'
#' 6. `6_mirna_precursor_annot_ssc11.csv`
#' 
#' **Output File Directory:** `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#'
#' **Output File(s):** `2_miRNA_eqtl_summary.Rdata`
#'
#' **Table of contents:**
#'
#' 1. [Objectives](#objectives)
#' 2. [Install libraries](#install-libraries)
#' 3. [Load data](#load-data)
#' 4. [Analysis](#analysis)
#' 5. [Save data](#save-data)
#'
#' ## Objectives
#'
#' The objective of this script is to summarize the number of eQTL peaks per miRNA output from the first eQTL analysis of the 174 F2 MSUPRP pig miRNA expression profiles.
#'
#' ## Install libraries
#'
rm(list=ls())

#' Load eqtl function Rdata containing stb function, which will summarize the eQTL results
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")
ls()

library(limma)
library(edgeR)
library(gwaR)
library(plyr)
library(corrplot)

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts/")

#' ## Load data
#'
#' Load the gblup and eQTL output:
load("../1_gblup_results_summary.Rdata")

load("../2_gwa_results.Rdata")

# #' Load the dge object to obtain the mature miRNA annotation:
# load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")

#' Load the MSUPRP_miRNA gpdata object for the mapping information:
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")

#' Load the updated annotation file I created for my miR-eQTL miRNAs:
annotation<-read.csv("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.csv", header=TRUE, na.strings="NA", colClasses=c("character","character","character","numeric","numeric","numeric","numeric","character","character","character"))
annotation<-annotation[,1:10]
save(annotation, file="../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
#' Load the Z matrix:
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")
ls()

#' ## Analysis
#'
#' ### Summarize the heritability of the miRNAs, output from GBLUP:
#'
#' The average heritability of all the miRNA expression profiles:
mean(summary.rst.gblup$h2)
summary(summary.rst.gblup$h2)

#' The average heritability of the miRNAs found to be significantly heritable at FDR < 0.05:
#'
#' How many significantly heritable miRNAs found at FDR < 0.05? (Should match from eQTL scan output md file)
sum(summary.rst.gblup$qvalue<0.05)

#' Extract the significantly heritable miRNAs from the summary.rst.gblup dataset and calculate mean h2:
sigh2<-summary.rst.gblup[summary.rst.gblup$qvalue<0.05,]
dim(sigh2)
mean(sigh2$h2)
summary(sigh2$h2)

#' Define the minimum p-value that is not significant (based on q-value < 0.05) as the threshold for plotting significant points
summary(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh<-min(summary.rst.gblup$lrtpvalue[summary.rst.gblup$qvalue >= 0.05])
sigthresh

#' Plot h2 vs -log10(p-value) like before to determine trend in significance and h2:

#+ vol_sig_h2, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))

plot(summary.rst.gblup$h2, -log10(summary.rst.gblup$lrtpvalue),
    xlab = expression("Heritability"~(h^{2})),
    ylab = "-log10(p-value)",
    main = "Significance vs Heritability")
points(summary.rst.gblup$h2[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       -log10(summary.rst.gblup$lrtpvalue)[-log10(summary.rst.gblup$lrtpvalue)>(-log10(sigthresh))],
       pch=19,col="red")
abline(a = -log10(sigthresh), b = 0, lty = 5)

#' ---
#'
#' ### The Summary Table of GWAS Results:
#'
#' Assign the correct names to the different objects:
map <- MSUPRP_miRNA$map
colnames(map)
head(annotation)
colnames(annotation)

#' Annotation must contain columns "chr", "start", "end", and "strand"
colnames(annotation)<-c("Name","gene", "transc", "chr","start","end","width","strand","type","Alias")
annotation$mid.mir<-round(rowMeans(annotation[,c("start", "end")], na.rm=TRUE))
head(annotation)

#' Add the mid.mir to the map object for later use in manhattan plots (arrow where midpoint of miRNA precursor lies)
#' 
#' Need to remove the 2nd map position of miR-128 for this analysis:
#' 
annotation<-annotation[-4,]

#' Extract the chromosome and the position of the miRNA precursor and build a data.frame to add to the map object later:
mirpos<-data.frame(chr=annotation$chr,
	pos=annotation$mid.mir, row.names=annotation$Name)
str(mirpos)
head(mirpos)
dim(mirpos)

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

#' Create the summary table, using pergene = TRUE to get gene-wise eQTL peaks:
sum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=TRUE)
dim(sum.eqtl)
head(sum.eqtl)
#' How many unique miRNAs have eQTL?
length(unique(sum.eqtl$miRNA))
unique(as.character(sum.eqtl$miRNA))

#' miR-eQTL peaks at FDR<0.05:
sum.eqtl

sum.eqtl<-sum.eqtl[order(sum.eqtl$miRNA, sum.eqtl$chr.snp),]

#' #### Summary of GWAS results at FDR < 0.05
#'
#' Number of eQTL peaks per chromosome:
table(sum.eqtl$chr.snp)

#' Names of associated miRNAs:
table(sum.eqtl$miRNA)
#' Chromosomes of associated miRNAs:
table(sum.eqtl$chr.miR)

#' Names of associated peak markers:
table(as.character(sum.eqtl$SNP))

#' ---
#' 
#' ### Determining the ranges of associated SNPs per eQTL peak on SSC15 (for ISAG abstract):
#'
#' First, create the summary table at FDR 5% again, this time with pergene=F to identify all markers associated with each eQTL peak:
fullsum.eqtl<-eqtlsum(gwarst=rst.gwa,map=map,annot=annotation,threshold=0.05, pergene=FALSE)
dim(fullsum.eqtl)

#' Summarize the number of SNPs associated with each miRNA eQTL (some have multiple peaks)
numsnps<-by(fullsum.eqtl, as.character(fullsum.eqtl$miRNA), nrow)
numsnps<-ldply(numsnps, fun=NULL, id=names(numsnps))
colnames(numsnps) <- c("miRNA", "numsnps")
numsnps
sum(numsnps$numsnps)

#' --- 
#' 
#' Examine the correlation of the SNPs within a peak with the same q-value
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


#' ---
#' 
#' ### Extract peak range data from all miRNA eQTL peaks
#'
#' I can obtain information on the range of each peak based on miRNA (adapted from Deborah's function "peakrng"):
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

#' nmt = name transcript (in this case, miRNA); so, make a list of the miRNAs and loop through to get the peak range information for each miRNA with significant eQTL peaks
#' sumtb = output from the summary table function, with pergene = FALSE

sigmirnames <- unique(as.character(sum.eqtl$miRNA))
sigmirnames

mirpeaks<-data.frame(do.call(rbind, lapply(sigmirnames, peakrngmir, fullsum.eqtl)), sum.eqtl[,c("SNP","pos.snp","qvalue")])
mirpeaks

#' Separate the two miR-184 SNPs on SSC3:
fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$chr.snp == "3",]

map["DBWU0000430",]
rownames(map[map$pos=="9463123",])

mir184.1<-mirpeaks[8,]
mir184.1$range.peak<-0
mir184.1$max.pos<-mir184.1$min.pos
mir184.1$num.snp<-1
mir184.1

map["ASGA0016793",]
rownames(map[map$pos=="126505748",])

mir184.2<-mirpeaks[8,]
mir184.2$range.peak<-0
mir184.2$min.pos<-mir184.2$max.pos
mir184.2$num.snp<-1
mir184.2$SNP<-rownames(map[map$pos=="126505748",])
mir184.2$pos.snp<-map["ASGA0016793","pos"]
mir184.2$qvalue<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP == "ASGA0016793","qvalue"]
mir184.2


mirpeaks[8,]<-mir184.1

mirpeaks[(nrow(mirpeaks)+1),]<-mir184.2

mirpeaks<-mirpeaks[order(mirpeaks$miRNA, mirpeaks$chr.snp),]
rownames(mirpeaks)<- NULL

mirpeaks

#' Add this change to sum.eqtl object:
mir184.1
mir184.2

tmp<-sum.eqtl[8,1:7]
tmp<-cbind(tmp, 
    SNP=as.character(mir184.2$SNP),
    chr.snp=mir184.2$chr.snp, 
    pos.snp=mir184.2$pos.snp, 
    snp.sign=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "snp.sign"], 
    pvalue=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "pvalue"],
    qvalue=fullsum.eqtl[fullsum.eqtl$miRNA=="ssc-miR-184" & fullsum.eqtl$SNP=="ASGA0016793", "qvalue"])

tmp

sum.eqtl<-rbind(sum.eqtl, tmp)
sum.eqtl<-sum.eqtl[order(sum.eqtl$miRNA, sum.eqtl$chr.snp),]
rownames(sum.eqtl)<-NULL
sum.eqtl

sum.eqtl[sum.eqtl$miRNA=="ssc-miR-184",]
mirpeaks[mirpeaks$miRNA=="ssc-miR-184",]

rm(mir184.1)
rm(mir184.2)
rm(tmp)

#' ---
#' 
#' ### Creating Manhattan plots of the six miRNA with the highest numbers of associated SNP markers (for ISAG poster)
#'
#' First, convert the map positions to absolute map positions using Deborah's function "absmap"
#' 
#' Notice how the map object's positions are relative to chromosome:
head(map)
dim(map)

#' Add the map positions of the miRNA with eQTL
map.full<-rbind(map[,1:2], mirpos[sigmirnames,])
dim(map.full)
head(map.full)
if(nrow(map.full) - nrow(map) != length(sigmirnames)) stop ("miRNA map positions not added correctly")

#' Use the absmap function to convert the chromosomal map positions to absolute map positions:
absposmap<-absmap(map.full)
head(absposmap)
tail(absposmap)
#' Notice it didn't include the miRNA with no map position (ssc-miR-140-5p)
absposmap[which(names(absposmap) %in% sigmirnames)]


#' Divide by 1e6 to get absolute position in Mb (nicer x-axis for plots)
head(absposmap/1e6)

#' Use sigpval function (from DV) to calculate the significant p-value cutoff for plotting manhattan plots for each miRNA
#' 
#' Then, provide the vector of miRNA names and the absposmap object to the manhpt function (also from Deborah's func_eqtl.Rdata) and loop through the vector of miRNA names to create the Manhattan plots:

#' ## Visualize
#'
#' Correlation plots for miReQTL SNPs with same qvalue within peak:
#'
#+ corrplot_mirSNP, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in 1:length(names(corsnp))){
    corrplot(corsnp[[i]], method="number",mar=c(2,2,2,2), cl.cex=0.7, tl.cex=0.8, tl.col="black", title=names(pkmirqtl[i]))
}

#' Look again at the plots with many SNPs with similar qvalues:
#+ corrplot_mirSNP2, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in names(pkmirqtl)[c(10, 12, 15)]){
    corrplot.mixed(corsnp[[i]], tl.pos="lt", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", title=names(pkmirqtl[i]))
}

#' ---
sigmirnames
#+ man_plot_pval, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames){
	pvalcutoff<-sigpval(i, pvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.pval"], qvalues=rst.gwa[rst.gwa$miRNA==i,"gwa.qval"], fdr=0.05)
	pvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.pval"]
	names(pvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=pvals,map=map.full,annotation=NULL,pvalues=TRUE,cutoff=pvalcutoff, arrow=TRUE)

}

#' The Manhattan Plots of q-values should be put on equal y-axes for comparison on poster.
#' 
#' Three of the peaks are strong enough signals to be on a y-axis of 0-10 (-log10(qval))
#' 
#' The remaining 14 peaks are less strong, on a y-axis of 0-4
#' 
sigmirnames4<-sigmirnames[c(1:4,6,8,9,11,12,14:17)]
sigmirnames4
sigmirnames8<-sigmirnames[c(5,7,10)]
sigmirnames8
sigmirnames10<-sigmirnames[13]
sigmirnames10

#+ man_plot_qval_y4, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames4){
	qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
	names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

	manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,4))

}

#+ man_plot_qval_y8, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames8){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,8))

}

#+ man_plot_qval_y10, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in sigmirnames10){
    qvals<-rst.gwa[rst.gwa$miRNA==i, "gwa.qval"]
    names(qvals)<-rst.gwa[rst.gwa$miRNA==i, "SNPid"]

    manhpt(nm=i,abspos=absposmap/1e6,rst=qvals,map=map.full,annotation=NULL,pvalues=FALSE,fdr=0.05, arrow=TRUE, ylim=c(0,10))

}


#' ## Save data
save(sum.eqtl, fullsum.eqtl, absposmap, map.full, mirpeaks, file = "../3_eqtl_summary_tables_maps.Rdata")
write.table(sum.eqtl, file="../3_eqtl_summary.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(fullsum.eqtl, file="../3_eqtl_fullsummary.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)