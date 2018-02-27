#' **Script:** `5_mirna_eqtl_examine_peaks.R`
#' 
#' **Directory of Code:**  `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts`
#' 
#' **Date:**  `12/11/17`
#' 
#' **Input File Directory:**  
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' 2. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/4_dge_G_objects/`
#' 
#' **Input File(s):** 
#' 
#' 1. `2_gwa_results.Rdata`, `3_eqtl_summary_tables_maps.Rdata`
#' 
#' 2. `6_mirna_precursor_annot_ssc11.Rdata`, `3_msuprp_mirna_gpdata.Rdata`, `4_normalized_dge_voom.Rdata`, `5_Z_G_miRNA.Rdata`
#' 
#' **Output File Directory:** 
#' 
#' 1. `/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/`
#' 
#' **Output File(s):** 
#' 
#' 1. `7_exam_eqtl_peaks_fix_snps.Rdata`
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
#' For each of the significant miRNA eQTL peaks identifyed in the miRNA eQTL scan this code will:
#' 
#' 1. Compute the variance explained by the markers +/- 1Mb surrounding eQTL peaks significanly associated to the miRNA expression in question
#' 
#' 2. Test the significance of the variance accounted for by the markers in these regions associated to the miRNA expression
#' 
#' 3. Fit the top significant marker per miRNA eQTL peak (by p-value) as a covariate in the model
#' 
#' 4. Observe resulting peaks (if any) after fixing the top significant SNP (p-value) per eQTL peak
#'
#' ## Install libraries
rm(list=ls())

setwd("/mnt/research/ernstc_lab/miRNA_eQTL_Sscrofa11/5_gblup_gwa_eqtl/scripts")

library(regress)
library (limma)
library (edgeR)
library(gwaR)
library(parallel)
library(qvalue)
library(corrplot)


#' ## Load data
#' 
#' Load DV's eqtl functions:
load("/mnt/research/pigeqtl/analyses/eQTL/paper/output/funct_eqtl.Rdata")

#' Load required data:
load("../2_gwa_results.Rdata")
load("../3_eqtl_summary_tables_maps.Rdata")
load("../../4_dge_G_objects/6_mirna_precursor_annot_ssc11.Rdata")
load("../../4_dge_G_objects/3_msuprp_mirna_gpdata.Rdata")
load("../../4_dge_G_objects/4_normalized_dge_voom.Rdata")
load("../../4_dge_G_objects/5_Z_G_miRNA.Rdata")

ls()

#' ## Analysis
#' 
#' Create covariates for analysis:
X <- data.frame(sex=as.factor(MSUPRP_miRNA$covar$sex),
	 selcrit=as.factor(MSUPRP_miRNA$covar$growth_group))
rownames(X) <- MSUPRP_miRNA$covar$id

#' Change miRNAs and markers to characters in correct format:
sum.eqtl$miRNA<-gsub("-",".", as.character(sum.eqtl$miRNA))
sum.eqtl$SNP<-as.character(sum.eqtl$SNP)

fullsum.eqtl$miRNA<-gsub("-",".", as.character(fullsum.eqtl$miRNA))
fullsum.eqtl$SNP<-as.character(fullsum.eqtl$SNP)

rownames(map.full)<-gsub("-",".", rownames(map.full))
tail(map.full)

map.full$chr<-as.numeric(map.full$chr)
map.full<-map.full[order(map.full$chr),]

mirpeaks$miRNA<-gsub("-", ".", mirpeaks$miRNA)

#' Matrix of miRNAs expressions and covariates
data <- data.frame(t(v$E), X, Z[,as.character(unique(mirpeaks$SNP))])

colnames(wtcen)<-gsub("-",".", colnames(wtcen))
head(colnames(wtcen))

#' ### Compute variance accounted associated markers
#' 
#' Vector of genes to test
rsp <- gsub("-",".", unique(mirpeaks$miRNA))
length(rsp)

#' Design for GBLUP
design <- c(~sex + selcrit)

#' GBLUP
#' 
#+ results='hide'
system.time({
    rst.gblup <- lapply(rsp, function(x) gblup(rsp=x, data=data,
        design=design, G=G, vdata=NULL, wt=wtcen, pos=c(T,T)))
    names(rst.gblup) <- rsp
})

#' Create map without miRNA positions:

map.nomir<-map.full[-grep("ssc*", rownames(map.full)),]
grep("ssc*", rownames(map.nomir))
dim(map.nomir)
str(map.nomir)

#' Isolate SNP regions +/- 1Mb from each miReQTL peak for determining variance accounted for by SNP with in associated genomic region.
#' 
#' Notice the range in coverage of 2Mb window is between 10 - 69 SNPs
eqtl.2Mbregions<-lapply(1:nrow(mirpeaks), function(x) inrange(chr=mirpeaks$chr.snp[x], start=mirpeaks$pos.snp[x]-1E6, end=mirpeaks$pos.snp[x]+1E6, map=map.nomir, single="pos", range=NULL))

#' 
#' Test Peak
system.time({
lrt.peak<-lapply(1:length(mirpeaks$miRNA), function(a) test.peak(gb=rst.gblup[[mirpeaks$miRNA[a]]], x=t(Z), peak_pos=rownames(eqtl.2Mbregions[[a]])))
})

# Eliminate NULL results
length(lrt.peak)
lrt.peak <- lrt.peak[unlist(lapply(lrt.peak, length)) == 3]
length(lrt.peak)

#' Compute qvalues
pval <- unlist(lapply(lrt.peak, function(x) x$pvalue))
qval <- qvalue(pval, lambda=0)$qvalue
sum(qval<0.05)

#' Merge the results of the LRT for each eQTL peak
varpeak <- do.call (rbind, lapply(1:length(lrt.peak), function(x)
		cbind(t(lrt.peak[[x]]$vars[1]), h2=lrt.peak[[x]]$vars[1, 3],
		lrt.peak[[x]]$llik, pvalue=lrt.peak[[x]]$pvalue)))

rownames(varpeak) <- paste(mirpeaks$miRNA, "row", rownames(mirpeaks), sep=".")

varpeak <- cbind(varpeak, qvalue=qval)

#' ---
#' 
#' Summary variance explaned by +/- 1Mb peak window
summary(varpeak$h2)

varpeak

#' Investigate the miRNAs that yield NAs in GWA:
rst.gblup$"ssc.miR.184"
rst.gblup$"ssc.miR.6782.3p"
rst.gblup$"ssc.miR.7135.3p"
rst.gblup$"ssc.miR.9810.3p"

#' ### GWA fixed top SNP per eQTL peak
#' 
#' Genome Wide Association fixing significant SNPs per peak for each miRNA expression containing one or more eQTL

#+ results='hide'
Z <- t(Z)
system.time({
	rst.gwa.pkfx <- lapply(1:length(mirpeaks$miRNA), function(a) run.gwa(rsp=mirpeaks$miRNA[a], data=data,
		design=as.formula(paste("~sex + selcrit +", as.character(mirpeaks[a, "SNP"]))), G=G, vdata=NULL,
		wt=wtcen, x=Z[!rownames(Z) %in% mirpeaks[a,"SNP"],], LRT=F, threshold = 0.05, returnz = T,
        pos=c(T,T)))
names(rst.gwa.pkfx)<-1:length(rst.gwa.pkfx)
})

#' Calculate pvalues from GWA Zscores
gwa.pv <- lapply(rst.gwa.pkfx, getpvalue, log.p=F, is.z=T)

#' Merge the results of the GWA into a matrix
rst.gwa.pkfx <- do.call(cbind, rst.gwa.pkfx)

#' Gene-wise Multiple Test Correction (FDR) for GWA pvalues (compute qvalues)
system.time({
	gwa.qv <- do.call(cbind, lapply(gwa.pv, function(x) qvalue(x)$qvalues))
})

#' Merge the pvalues of the GWA into a matrix
gwa.pv <- do.call(cbind, gwa.pv)

#' Merge results of eQTL analysis into one R object
nms.gwa <- mirpeaks$miRNA
names(nms.gwa) <- colnames(gwa.pv) 
exam.eqtl <- list(varpeak=varpeak, gwa=rst.gwa.pkfx, gwa.pval=gwa.pv, gwa.qval=gwa.qv, nms.gwa=nms.gwa)

#' Number of genes still containing a significant peak
threshold <- 0.05
qval <- exam.eqtl$gwa.qval
sig <- qval < threshold
colSums(sig)
gw <- colSums(sig)[colSums(sig)!=0]
length(gw)
summary(gw)
gw
#' Remove NAs (ssc.miR.184 (column 11), ssc-miR-7135-3p (column 16), ssc-miR-874 (column 17), and ssc-miR-9810-3p (column 22))
gw<-gw[!is.na(gw)]
gw


#' Reduce qval matrix, retain only miRNAs with a significant association remaining
qval <- qval[,names(gw)]

nmt <- exam.eqtl$nms.gwa[colnames(qval)]
dim(qval)

#' Create pval object to extract correct peak SNP:
pval2<-exam.eqtl$gwa.pval[,colnames(qval)]

#' Separate nmt into groups based on y-axis appropriate for Manhattan plots:
nmt
nmt4<-nmt[c(1,2)]

nmt6<-nmt[c(5,6)]

nmt8<-nmt[c(3,4)]

nmt10<-nmt[7]

#' ---
#' Investigate why I'm getting NAs in the rst.gwa.pkfx:
SNP11<-c(as.character(mirpeaks[11,"SNP"]), rownames(rst.gwa.pkfx)[is.na(rst.gwa.pkfx[,11])])
corZmiR11<-cor(t(Z[SNP11,]))

SNP16<-c(as.character(mirpeaks[16,"SNP"]),rownames(rst.gwa.pkfx)[is.na(rst.gwa.pkfx[,16])])
corZmiR16<-cor(t(Z[SNP16,]))

SNP17<-c(as.character(mirpeaks[17,"SNP"]), rownames(rst.gwa.pkfx)[is.na(rst.gwa.pkfx[,17])])
corZmiR17<-cor(t(Z[SNP17,]))

SNP22<-c(as.character(mirpeaks[22,"SNP"]), rownames(rst.gwa.pkfx)[is.na(rst.gwa.pkfx[,22])])
corZmiR22<-cor(t(Z[SNP22,]))

#' It's because of the LD! Notice the extremely strong correlation between the SNPs in each case.
corZmiR11

corZmiR16

corZmiR17

corZmiR22

#' Identify the correlation of the SNPs in the four complete eQTL peaks yielding NAs:
#' 
#' miRNA-184 eQTL peak contains 46 SNPs:
snp.184<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.184","SNP"]
#' Only want the significant SNP on SSC7; this needs to be removed.
length(snp.184)
head(map.full[snp.184,])
snp.184<-snp.184[4:length(snp.184)]
length(snp.184)
head(snp.184)
#' Calculate correlation:
corsnp.184<-cor(t(Z[snp.184,]))

#' miRNA-6782-3p eQTL peak contains 4 SNP
snp.6782<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.6782.3p","SNP"]
length(snp.6782)

corsnp.6782<-cor(t(Z[snp.6782,]))

#' miRNA-7135-3p eQTL peak contains 14 SNP:
snp.7135<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.7135.3p","SNP"]
length(snp.7135)

#' Calculate correlation:
corsnp.7135<-cor(t(Z[snp.7135,]))

#' miRNA-9810-3p eQTL peak contains 2 SNP:
snp.9810<-fullsum.eqtl[fullsum.eqtl$miRNA=="ssc.miR.9810.3p","SNP"]
length(snp.9810)

#' Calculate correlation:
corsnp.9810<-cor(t(Z[snp.9810,]))

#' #### These correlations will be visualized using corrplot package later on.

#' ---
#' 
#' Create tables to summarize results of fixing peak eQTL:
annotation$miRNA<-gsub("-",".", annotation$miRNA)
annotation<-annotation[-4,]
rownames(annotation)<-annotation$miRNA

stb.nmiR<-function(nm,qval,pval,map,annotation,Z,threshold=0.01,gene.name="genes",pergene=TRUE){
    idx <- qval < threshold
    snp <- names(qval)[idx]
    snp.effect <- ifelse(Z[idx] < 0, x<-"-", x<-"+")
    gene <- rep(nm, sum(idx))

    rst <- data.frame(miR=gene, chr.miR=annotation[gene,"chr0"],
        start.miR=annotation[gene,"start"]/1e6, 
        end.miR=annotation[gene,"end"]/1e6,
        strand=annotation[gene,"strand"], 
        SNP=snp, 
        chr.snp=map[snp,"chr"], 
        pos.snp=map[snp,"pos"], 
        snp.effect=snp.effect,
        pvalue=pval[idx],
        qvalue=qval[idx], 
        row.names=NULL)

    if(pergene){
        id<-unique(rst[,"miR"])
        x<-list()

        for (i in id){
            a<-rst[rst[,"miR"]==i,]

            if (length(unique(a[,"chr.snp"]))==1){
                a<-a[order(a[,"pvalue"])[1],]

                } else {

                b<-a[order(a[,"pvalue"]),]
                a<-list()

                    for (j in unique(b$chr.snp)){
                        a[[j]]<-b[b[,"chr.snp"]==j,][1,]
                    }

                a<-do.call(rbind,a)
                }

        x[[i]]<-a

        }

    rst<-do.call(rbind,x)
    }
    rownames(rst)<-NULL
    return(rst)
}

#' Summary Table for all miRNA-marker associations
sumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[[x]],
    qval=qval[,x], pval=pval2[,x], map=map.full, annotation=annotation, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=F))
names(sumtb.exam) <- names(nmt)
length(sumtb.exam)
sumtb.exam

#' Summary table for all miRNA eQTL peaks (pergene=T)
rsumtb.exam <- lapply(names(nmt), function(x) stb.nmiR(nm=nmt[x],
    qval=qval[,x], pval=pval2[,x], map=map.full, annotation=annotation, Z=exam.eqtl$gwa[,x],
    threshold=0.05, pergene=T))
names(rsumtb.exam) <- names(nmt)
length(rsumtb.exam)
rsumtb.exam

#' Identify which SNP has the minimum p-value for reporting
lapply(1:length(sumtb.exam), function(x) sumtb.exam[[x]][sumtb.exam[[x]]$pvalue==min(sumtb.exam[[x]]$pvalue),])
lapply(1:length(rsumtb.exam), function(x) rsumtb.exam[[x]][rsumtb.exam[[x]]$pvalue==min(rsumtb.exam[[x]]$pvalue),])

#' Create data.frame of eQTL peaks for diferentiationg between local and distant regulators:

peakrngmir<-function(nmtb, sumtb){

# Positions of SNPs in eQTL peak
    nmtb <- nmtb
    map.CI <- sumtb[sumtb$miR==nmtb,c("SNP","chr.snp","pos.snp","qvalue")]

# Number of associated markers by chromosome
    chr <- table(map.CI$chr.snp)
    cat("Number of markers per chromosomal peak for",nmtb,":")
    print(chr)
    idx <- as.numeric(names(chr))

    min.pos <- unlist(lapply(idx, function(x) min(map.CI[map.CI$chr.snp == x,"pos.snp"])))
    max.pos <- unlist(lapply(idx, function(x) max(map.CI[map.CI$chr.snp == x,"pos.snp"])))

    start.miR <- rep(sumtb[sumtb$miR == nmtb,"start.miR"][1], length(idx))
    end.miR <- rep(sumtb[sumtb$miR == nmtb,"end.miR"][1], length(idx))

# Identify the position of marker extreams for each peak
    peaks <- data.frame(miRNA=rep(nmtb,length(idx)),
                chr.miR=rep(sumtb[sumtb$miR == nmtb,"chr.miR"][1], length(idx)),
                start.miR=start.miR, end.miR=end.miR, range.miR=end.miR-start.miR,
                chr.snp=idx, range.peak=max.pos-min.pos,
                min.pos=min.pos, max.pos=max.pos, num.snp=as.vector(chr))

    return(peaks)
}

# data.frame(peakrngmir(nmtb=nmt[[3]], sumtb=sumtb.exam[[3]]),rsumtb.exam[[3]][,c("SNP","pos.snp")])

#+ results='hide'
peaks.exam <- lapply(names(nmt), function(x) data.frame(peakrngmir(nmtb=nmt[x],sumtb=sumtb.exam[[x]]),rsumtb.exam[[x]][,c("SNP","pos.snp")]))
names(peaks.exam) <- names(nmt)

#' Compare the previous peaks with those obtained after fixing the sig assoc snps
for (i in names(nmt)){
cat(nmt[i], '\n')
cat('\n','All the peaks for this miRNA:','\n')
print(mirpeaks[mirpeaks$miRNA == nmt[i],])
cat('\n','Fixed this top snp:','\n')
print(mirpeaks[names(nmt),][i,])
cat('\n','Peak after fixing top SNP:','\n')
print(peaks.exam[[i]])
cat('\n','\n','\n')
}

#' ---
#' 
#' #### Summary
#' 

#' These GWA resulted in retained miRNA eQTL peaks after fixing peak SNP:
mirpeaks[names(nmt),]

#' These GWA resulted in elimated miRNA eQTL peaks after fixing peak SNP:
mirpeaks[!names(exam.eqtl$nms.gwa) %in% names(nmt),]

#' These GWA resulted in NAs after fixing peak SNP, due to highly correlated SNPs in peaks:
mirpeaks[c(11, 16, 17, 22),]

#' Note two SNPs of interest, where fixing one of multiple peak SNPs eliminated ALL peaks, while eliminating other peak SNPs did not remove all peaks:
#' 
#' 1. ASGA0034057 on SSC7 for miR-184 (four peaks)
#' 
#' 2. ALGA0016550 on SSC2 for miR-874 (two peaks)
#' 
#' ---
#' 
#' Additional investigation into distribution of p-values for original GWA and GWA with fixed peak SNP:

rst.gwa$miRNA<-gsub("-", ".", rst.gwa$miRNA)

lapply(unique(mirpeaks$miRNA), function(x) hist(rst.gwa[rst.gwa$miRNA==x, "gwa.pval"], breaks=20))

lapply(unique(mirpeaks$miRNA), function(x) summary(rst.gwa[rst.gwa$miRNA==x, "gwa.pval"]))

lapply(1:length(mirpeaks$miRNA), function(x) hist(exam.eqtl$gwa.pval[,x]))

lapply(1:length(mirpeaks$miRNA), function(x) summary(exam.eqtl$gwa.pval[,x]))

#' ## Visualize
#' 
#' ### Correlation of SNPs giving NAs in the gwa after fixing the peak SNP:
#' 
#' miR-184
#+ cor_plot_miR184, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.184, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.7, tl.col="black", title="miR-184 peak SNP ASGA0034057")

#' miR-6782-3p
#+ cor_plot_miR6782, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.6782, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.7, tl.col="black", title="miR-6782-3p peak SNP DIAS0000707")

#' miR-7135-3p
#+ cor_plot_miR7135, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.7135, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", tl.srt=45, title="miR-7135-3p peak SNP ALGA0124095")

#' miR-9810-3p
#+ cor_plot_miR9810, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
corrplot(corsnp.9810, method="color", type="upper", mar=c(2,2,2,2), tl.cex=0.8, tl.col="black", tl.srt=45, title="miR-9810-3p peak SNP ALGA0030853")


#' ### Manhattan Plots
par(oma=c(2,2,2,2))
#+ man_plot_qval_y4, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in names(nmt4)){
    manhpt(nm=nmt[i], abspos=absposmap/1e6, rst=qval[,i],
	map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,4))
    points(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), col="red")
    text(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), labels=sum.eqtl[i,"SNP"], pos=4)
    }

#+ man_plot_qval_y6, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in names(nmt6)){
    manhpt(nm=nmt[i], abspos=absposmap/1e6, rst=qval[,i],
    map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8))
    points(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), col="red")
    text(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), labels=sum.eqtl[i,"SNP"], pos=3)
    }

#+ man_plot_qval_y8, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
for(i in names(nmt8)){
    manhpt(nm=nmt[i], abspos=absposmap/1e6, rst=qval[,i],
    map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8))
    points(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), col="red")
    text(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), labels=sum.eqtl[i,"SNP"], pos=2)
    }

for(i in names(nmt10)){
    manhpt(nm=nmt[i], abspos=absposmap/1e6, rst=qval[,i],
    map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10))
    points(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), col="red")
    text(x=absposmap[sum.eqtl[i,"SNP"]]/1e6, y=(-log10(sum.eqtl[i, "qvalue"])), labels=sum.eqtl[i,"SNP"], pos=4)    
}

#' Checking the other peaks, where fixing the top SNP eliminated very large peaks:
#' 
#' miR-874
#+ man_plot_miR874, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
manhpt(nm="ssc.miR.874", abspos=absposmap/1e6, rst=exam.eqtl$gwa.qval[,18], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,10))
points(x=absposmap[sum.eqtl[18,"SNP"]]/1e6, y=(-log10(sum.eqtl[18, "qvalue"])), col="red")
text(x=absposmap[sum.eqtl[18,"SNP"]]/1e6, y=(-log10(sum.eqtl[18, "qvalue"])), labels=sum.eqtl[18,"SNP"], pos=4)

#' miR-429
#+ man_plot_miR429, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
manhpt(nm="ssc.miR.429", abspos=absposmap/1e6, rst=exam.eqtl$gwa.qval[,14], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8))
points(x=absposmap[sum.eqtl[14,"SNP"]]/1e6, y=(-log10(sum.eqtl[14, "qvalue"])), col="red")
text(x=absposmap[sum.eqtl[14,"SNP"]]/1e6, y=(-log10(sum.eqtl[14, "qvalue"])), labels=sum.eqtl[14,"SNP"], pos=4)

#' miR-184
#+ man_plot_miR184, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
manhpt(nm="ssc.miR.184", abspos=absposmap/1e6, rst=exam.eqtl$gwa.qval[,11], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,8))
points(x=absposmap[sum.eqtl[11,"SNP"]]/1e6, y=(-log10(sum.eqtl[11, "qvalue"])), col="red")
text(x=absposmap[sum.eqtl[11,"SNP"]]/1e6, y=(-log10(sum.eqtl[11, "qvalue"])), labels=sum.eqtl[11,"SNP"], pos=4)

#' miR-9785-5p
#+ man_plot_miR9785-5p, dpi=300, dev='tiff', dev.args=list(tiff = list(compression = 'lzw'))
manhpt(nm="ssc.miR.9785-5p", abspos=absposmap/1e6, rst=exam.eqtl$gwa.qval[,21], map=map.full[rownames(qval),], fdr=threshold, ylim=c(0,4))
points(x=absposmap[sum.eqtl[21,"SNP"]]/1e6, y=(-log10(sum.eqtl[21, "qvalue"])), col="red")
text(x=absposmap[sum.eqtl[21,"SNP"]]/1e6, y=(-log10(sum.eqtl[21, "qvalue"])), labels=sum.eqtl[21,"SNP"], pos=2)


#' ## Save data
#' 
#' Save eQTL results summary:
save(varpeak, lrt.peak, exam.eqtl, sumtb.exam, rsumtb.exam, peaks.exam, file="../7_exam_eqtl_peaks_fix_snps.Rdata")
write.table(varpeak, file="../8_var_expl_peaks.txt", quote=FALSE, row.names=TRUE, col.names=TRUE)