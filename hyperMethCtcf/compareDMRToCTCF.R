# recommended: R-bundle/20180906-foss-2017a-R-3.5.1
setwd("/omics/odcf/analysis/OE0290_projects/Ependymoma/results/wgbs/methAnalysis")

library(genefilter)
library(rtracklayer)


# PFA vs RELA only (without chrX)
dmrs <- read.table("/omics/odcf/analysis/OE0290_projects/Ependymoma/results/wgbs/dmrIdentification/PFA_vs_RELA/metilene_PFA_vs_RELA_specific_qval.0.05.out")

colnames(dmrs)[1:3] <- c("chr", "start", "end")
dmrGR <- GRanges(dmrs[,1:3])

# CTCF PFA vs RELA
diffPeaks <- read.table("/omics/odcf/analysis/OE0290_projects/Ependymoma/results/chip_seq_sequencing/CTCF/peaks/analysisDiffBind/PFA_vs_RELA.CTCF.cuffDiff_result.131119.txt", header=1)


diffPeakGr <- GRanges(diffPeaks[,1:3])
# not required for mehtylation
seqlevels(diffPeakGr) <- sub('chr','',seqlevels(diffPeakGr))

# find common
ovlps <- findOverlaps(diffPeakGr, dmrGR)

selPeaks <- diffPeaks[unique(queryHits(ovlps)), ]
selPeaksGr <- diffPeakGr[unique(queryHits(ovlps))]


colnames(dmrs) <- c("chr", "start", "end", "qval", "methDiff", "numCpG", "methPFA", "methRELA")
diffPeaksComb <- cbind(diffPeaks[queryHits(ovlps), ], dmrs[subjectHits(ovlps),] )

cols = rep("black", nrow(diffPeaksComb))
cols[ (diffPeaksComb$Fold < 0)  & (diffPeaksComb$V5 > 0) ] = "red"


write.table(diffPeaksComb, "DMRs_PFA_vs_RELA_ovlp_diff_peaks_CTCF.131119.txt", sep="\t",row.names=F, quote=F)

summary( (diffPeaksComb$Fold < 0)  & (diffPeaksComb$methDiff > 0) )
summary( (diffPeaksComb$Fold < 0)  & (diffPeaksComb$methDiff < 0) )

## Figure for manuscript

pdf("DMRs_CTCF_diff_peaks.PFA.061119.pdf", width=8, height=6)
plot( diffPeaksComb$Fold, diffPeaksComb$V5, pch=19, col=cols, 
    main = "Overlapping PFA-specific CTCF diff peaks and DMRs", 
    xlab = "log2FC CTCF", ylab = "meth diff"   )
abline(h = 0,lty=2)
abline(v = 0,lty=2)
dev.off()



fullDF <- read.table("/omics/odcf/analysis/OE0290_projects/Ependymoma/results/chip_seq_sequencing/CTCF/peaks/analysisDiffBind/PFA_vs_RELA.CTCF.cuffDiff_result_full.131119.txt", sep="\t", header=1)

# for variance figure result is more relaxed: p.val is 0.1

fullDF <- fullDF[ fullDF$p.value < 0.1,]


fullPeakGr <- GRanges(fullDF[,1:3])
seqlevels(fullPeakGr) <- sub('chr','',seqlevels(fullPeakGr))

ovlps <- findOverlaps(fullPeakGr, dmrGR)

diffPeaksComb <- cbind(fullDF[queryHits(ovlps), ], dmrs[subjectHits(ovlps),] )

sel <- diffPeaksComb[diffPeaksComb$V5 > 0,] # PFA hyper
sel <- diffPeaksComb[diffPeaksComb$V5 < 0,] # RELA hyper

# PFA
summary( (diffPeaksComb$V5 > 0) & (diffPeaksComb$Fold < 0) ) 
# RELA 
summary( (diffPeaksComb$V5 < 0) & (diffPeaksComb$Fold > 0) )



cols = rep("black", nrow(fullDF))
cols[ rownames(fullDF) %in% rownames(sel) ] <- "orange"


# MAIN FIGURES

pdf("CTCF_diff_peaks_pval.260320.pdf", width = 8, height=6 )

plot( fullDF$Fold, -log10(fullDF$p.value) , pch=19 , cex=1,
      main = "PFA vs RELA CTCF diff peaks\nhypermeth effect", col = cols,
      xlab = "log2FC", ylab= "-log10(p-val)" ) 
abline(v = 0, lty=2)
legend("topright", # title="", # cex=
       legend=c("PFA hypermeth"), fill = c("red") )
dev.off()

require(ggplot2)

pdf("CTCF_diff_peaks_pval.271021.pdf", width = 8, height=6 )
ggplot(fullDF, aes(Fold, -log10(p.value)))  + 
 geom_point(colour = cols, size = 1) + 
geom_vline(xintercept=0, linetype="dashed", color = "black") +
labs( title="PFA CTCF loss: hypermethylation effect" ) +
theme_classic() +
    theme(text = element_text(size=20), axis.text.y=element_text(size=18),
        axis.text.x = element_text(size=18))
dev.off()

# supplementary figure
# full selection
diffPeaksComb2 <- diffPeaksComb[   , c(1:3,9, 15:19)]

# filtering based on log2FC
diffPeaksComb2 <- diffPeaksComb[ (abs(diffPeaksComb$Fold) > 0.5)  , c(1:3,9, 15:19)]

# target : exclude non-focus
sel <- ( (diffPeaksComb$V5 > 0) & (diffPeaksComb$Fold < 0) ) | ( (diffPeaksComb$V5 < 0) & (diffPeaksComb$Fold > 0) ) 
diffPeaksComb2 <- diffPeaksComb[ sel  , c(1:3,9, 15:19)]


pdf("CTCF_variance_between_peaks_DMRs.010620.v2.pdf", width = 8, height=6 )
ggplot(diffPeaksComb2, aes(Fold, V5) ) + 
 geom_point(colour = "black", size = 1,alpha=0.5) + 
#geom_vline(xintercept=0, linetype="dashed", color = "grey") +
labs( title="PFA vs RELA", x="CTCF diff. peaks Fold Change (log2)",
        y="WGBS DMRs Variance" ) +
    
 theme_minimal() + geom_vline(xintercept=0) + geom_hline(yintercept=0 ) +  
    theme(text = element_text(size=20), axis.text.y=element_text(size=18),
        axis.text.x = element_text(size=18))
dev.off()



library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)

# NOTE: require diff. peaks GR with "chr" remaining!


pfaLossPeaksGr <- diffPeakGr[diffPeaks$Fold < 0 ]

# decrease result based on log2FC
#pfaLossPeaksGr <- diffPeakGr[diffPeaks$Fold < -2 ]
# top300
pfaLossPeaksGr <- diffPeakGr[1:300,]

# main focus tumor samples
sIds <- c("9EP1", "9EP9", "7EP18_neu" ,"11EP22", "7EP41", "9EP38" )


resList = list()

#wSize = 50 
wSize = 250
#wSize = 500

for (sId in sIds) {

    print(sId)
    sMeth <- import(paste0("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/wgbs/igvData/",sId,".meth.wig"), format="Wig")

    mat = normalizeToMatrix(sMeth, pfaLossPeaksGr , value_column = "score", mean_mode = "absolute",extend = 5000, w = wSize, target_ratio = 0.5, background = NA, smooth = TRUE)

    resList[[sId]] <- mat

}


# create 6 objects

obj1 <- EnrichedHeatmap(resList[[ sIds[1] ]], col = meth_col_fun, name = "methylation", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched( ylim=c(0.3,0.8), yaxis=FALSE,  gp = gpar(col = 1:5) ) ) ,  column_title = sIds[1] )

obj2 <- EnrichedHeatmap(resList[[ sIds[2] ]], col = meth_col_fun, name = "methylation 2", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched( ylim=c(0.3,0.8), yaxis=FALSE,  gp = gpar(col = 1:5))) ,    column_title = sIds[2], show_heatmap_legend = FALSE ) 
 
obj3 <- EnrichedHeatmap(resList[[ sIds[3] ]], col = meth_col_fun, name = "methylation 3", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched( ylim=c(0.3,0.8), yaxis=FALSE,  gp = gpar(col = 1:5))) ,    column_title = sIds[3], show_heatmap_legend = FALSE ) 

obj4 <- EnrichedHeatmap(resList[[ sIds[4] ]], col = meth_col_fun, name = "methylation 4", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0.3,0.8), yaxis=FALSE, gp = gpar(col = 1:5))) ,    column_title = sIds[4], show_heatmap_legend = FALSE ) 

obj5 <- EnrichedHeatmap(resList[[ sIds[5] ]], col = meth_col_fun, name = "methylation 5", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0.3,0.8), yaxis=FALSE,gp = gpar(col = 1:5))) ,    column_title = sIds[5], show_heatmap_legend = FALSE ) 
 
obj6 <- EnrichedHeatmap(resList[[ sIds[6] ]], col = meth_col_fun, name = "methylation 6", axis_name_rot = 90,  top_annotation = HeatmapAnnotation(lines = anno_enriched(ylim=c(0.3,0.8), gp = gpar(col = 1:5))) ,    column_title = sIds[6], show_heatmap_legend = FALSE ) 
 
# 6 samples togeth
ht_list = obj1 + obj2 + obj3 + obj4 + obj5 + obj6

# main figure selection
resName = paste0("methStatus.PFA_CTCF_loss.top300.wSize_",wSize,".190821.pdf")


pdf(resName, width=16,height=12)
draw(ht_list)
dev.off()

resName = paste0("methStatus.PFA_CTCF_loss.wSize_",wSize,".310520.svg")
svg(resName, width=16,height=12)
draw(ht_list)
dev.off()


### Connect to HIC data 
### approach 1: check if within loop of EAG

selPeaks <- read.table("DMRs_PFA_vs_RELA_ovlp_diff_peaks_CTCF.131119.txt",sep="\t", header =T)
selPeaksGr <- GRanges(selPeaks[,1:3])


# PFA specific main
loopsDf <- read.table("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_full.pval_0.01.enh_contacts_closest.cor_analysis.txt", sep="\t", header=1)

# RELA specific main
loopsDf <- read.table("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.RELA_combined_full.pval_0.01.enh_contacts_closest.cor_analysis.txt", sep="\t", header=1)


loopsGR <- GRanges(data.frame(chr = loopsDf$chr1, start=loopsDf$fragmentMid1, end=loopsDf$fragmentMid2 ) )
#seqlevels(loopsGR) <- sub('chr','',seqlevels(loopsGR))

resWithinLoops = findOverlaps(selPeaksGr, loopsGR, type="within")

# check
#selPeaks[queryHits(resWithinLoops),]

resDf <- cbind(loopsDf[subjectHits(resWithinLoops),], selPeaks[queryHits(resWithinLoops),] )
resDf <- resDf[resDf$Cor > 0.5,] # group-specific only

write.table(resDf, "FitHiC_res5000.PFA_combined_full.pval_0.01.loops_with_CTCF_DMR_in_between.041219.txt",quote=F, sep="\t")

write.table(resDf, "FitHiC_res5000.RELA_combined_full.pval_0.01.loops_with_CTCF_DMR_in_between.050320.txt",quote=F, sep="\t")


## approach 2: work on SE regions adjustment 
## requies selPeaks
## here only PFA is shown, RELA demonstarted no match

pfaLoss <- selPeaks[ (selPeaks$methDiff > 0) & (selPeaks$Fold < 0), ]

targSites <- GRanges(pfaLoss[,1:3] ) 

# find connection to motifs

dmrGr <- GRanges(seqnames = paste0("chr",pfaLoss$chr), ranges = IRanges(start=pfaLoss$start.1, end = pfaLoss$end.1)) 

checkRegions <- pintersect(targSites, dmrGr)
motifGr <- import.bed("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/targ_CTCF_motifs.hg19.full.bed.gz")
seqlevels(motifGr) <- paste0("chr",seqlevels(motifGr))
ovlpMeasure <- countOverlaps(checkRegions, motifGr)

pfaLoss$ctcfLoci <- as.character(checkRegions)
pfaLoss$ctcfMotifs <- ovlpMeasure


# load SE regions
seGr <- import.bed("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/chip_seq_sequencing/enhancer/combined_approach/Ependymoma_SE.bed")
names(seGr) <- paste0("SE",1:length(seGr))

# filter loops to have 
loops <- read.delim("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/enhancersAnalysis/fithic_res/FitHiC_res5000.PFA_combined_full.qval_0.05.txt")

loopsDf <- data.frame(chr=loops$chr1, start=loops$fragmentMid1,end=loops$fragmentMid2, bias1=loops$bias1,bias2=loops$bias2)
loopsGr <- GRanges(loopsDf)

ctcfInLoop <- findOverlaps(targSites,loopsGr,type="within")

targLoops <- loopsDf[ subjectHits(ctcfInLoop), ]

# for annotation
refCtcf <- pfaLoss[ queryHits(ctcfInLoop), ]
targLoopsExt <- cbind(targLoops, refCtcf)

# format correctly
targLoopsFormated <- data.frame(chr1=targLoops$chr, pos1=targLoops$start, chr2=targLoops$chr, pos2=targLoops$end)



# find loops that have SE on one side  and gene on other 
library(InTAD)

geneGr <- import.gff("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/star_index_gencode19/gencode.v19.annotation.genes_only.gtf")
names(geneGr) <- geneGr$gene_id

# these are empty matrices since correlation is ignored
cntMtx <- matrix(1, length(geneGr), 5)
rownames(cntMtx) <- names(geneGr)
colnames(cntMtx) <- paste0("s",1:5)
seMtx <- matrix(1, length(seGr), 5)
rownames(seMtx) <- names(seGr)
colnames(seMtx) <- colnames(cntMtx)

inTadObj <- newSigInTAD(seMtx, seGr, cntMtx, geneGr)

rownames(targLoops) <- paste0("loop",1:nrow(targLoops))
rownames(targLoopsExt) <- rownames(targLoops)

inTadObj <- combineWithLoops(inTadObj,targLoopsFormated,fragmentLength=5000)

require(dplyr)
resDf <- bind_rows(inTadObj@loopConnections)
resDf2 <- cbind(resDf,targLoops[resDf$Loop,] )

resDf2 <-cbind(resDf2, targLoopsExt[resDf$Loop,c("ctcfLoci", "ctcfMotifs")] )


resDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/ctcfLossHypomethSE/"

extDist = 0
resName <- paste0(resDir, "PFA_gene_SE_loop_connected.with_motifs.txt")
write.table(resDf2, resName, sep="\t", quote=F,row.names=F)

# this image without motifs
#save.image(paste0(resDir,"PFA_gene_SE_loop_initial.300321.RData"))
#load(paste0(resDir,"PFA_gene_SE_loop_initial.300321.RData"))





