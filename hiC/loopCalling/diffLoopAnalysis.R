# tested on R 4.0.1
# source : https://rpubs.com/caleblareau/diffloop_vignette

# test
#bed_dir <- system.file("extdata", "esc_jurkat", package="diffloopdata")

library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)

# main check
bed_dir <- "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/diffAnalysis/diffloop/PFA_vs_RELA"

# arima capture : DMSO vs Aza
bed_dir <- "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/CaptureHiC/loops_arcplots/diff_loop_input"

resDir = "/omics/odcf/analysis/OE0290_projects/Ependymoma/results/HiC/hicPro/CaptureHiC/loops_arcplots/"

resName = "diff_loops"
full <- loopsMake(bed_dir)


# tumors
grps <- c("RELA","RELA","PFA","RELA","PFA","PFA")
# de-methylation
grps <- c("Aza", "Aza", "Aza", "DMSO", "DMSO")


full <- updateLDGroups(full, grps)
head(full,4)

# filtering
# this version is more optimal
full <- subsetLoops(full, full@rowData$loopWidth >= 5000) 

p1 <- loopDistancePlot(full)

#samples <- c("11EP22", "4EP53", "7EP18", "7EP41", "9EP1", "9EP9")

# default
pcp2 <-pcaPlot(full)

# adjusted
pcp2 <-pcaPlot(full) +# geom_text_repel(aes(label=samples)) +
    scale_x_continuous(limits = c(-120, 120)) + ggtitle("PC Plot with Size Factor Correction") + theme(legend.position="none")

qcName <- paste0(bed_dir, ".QC_check.pdf")
pdf(qcName, width=8,height=6)
p1
pcp2
dev.off()


km_res <- quickAssoc(full)

# manual check
fullRes <- km_res@rowData
ids <- as.numeric(rownames(fullRes[fullRes$PValue < 0.05,]))
km_res[ids,]

showRegion <- function(res, id) {
    cons <- res@interactions[id,]
    paste(res@anchors[cons[1]],res@anchors[cons[2]])
}

# main result extraction
getFullRes <- function(res, minPval) {
    rObj <- res@rowData
    ids <- as.numeric(rownames(rObj[rObj$PValue < minPval,]))
    resDf <- NULL
    for ( id in ids) {
        cons <-  res@interactions[id,]
        contact <- data.frame( pos1 = as.character(res@anchors[cons[1]]), 
                        pos2 = as.character(res@anchors[cons[2]]) )
        resDf <- rbind(resDf,contact)
    }
    
    resDf <- cbind(resDf, rObj[ids,])
    resDf
    
}


res <- getFullRes(km_res,0.1)
res2 <- res[ res$PValue < 0.05,]

resName <- paste0(bed_dir, ".diff_loops.min_pval_0.01.txt")
write.table(res,resName,quote=F,sep="\t")



resName2 <- paste0(resDir, "diff_loops.",resName, ".min_pval_0.05.txt")
write.table(res2, resName2, quote=F, sep="\t")


## arima results check

# Plot difference

resName4 <- paste0(bed_dir, ".diff_loops.min_pval_0.05.volcano_plot.pdf")

require(ggplot2)
cols = rep("orange",nrow(res))
cols[res$logFC < 0 ] <- "gold4"

pdf(resName4, width = 8, height=6 )
ggplot(res, aes(logFC, -log10(PValue)))  +
 #geom_point(colour = cols, size = 3,alpha=0.6) +
 geom_point(colour = cols,size = 0.5) +
geom_vline(xintercept=0, linetype="dashed", color = "black") +
labs( title="5_Aza vs DMSO:  3627 vs 27016 " ) +
theme_classic() +
    theme(text = element_text(size=20), axis.text.y=element_text(size=18),
        axis.text.x = element_text(size=18))
dev.off()



# requires res2 sorted + normMtx

mtx <- km_res@counts
normMtx <- t( t(km_res@counts) / km_res@colData$sizeFactor )

require(ComplexHeatmap)
require(viridis)
selCols <- magma(20)


nLoops = 500
d1 <- rownames(res2)[1:nLoops]

targMtx <- normMtx[as.numeric(d1), ]

resName5 <- paste0(resDir, "diff_loops.", resName, ".heatmap_top",nLoops,"_diff_loops.pdf")


pdf(resName5, width = 4, height = 9)
ht1 <- Heatmap(targMtx, #top_annotation = ha,
               col = selCols,
               #top_annotation_height = unit(1, "cm"),
               cluster_columns  = FALSE,
               # use control order
               cluster_rows = FALSE,
               # perform clustering
               #cluster_rows = TRUE,
               #clustering_method_rows = "ward.D2",
               column_title = paste0("Top ",nLoops," diff.loops"),
               show_column_names = TRUE,
               show_row_names = FALSE,
               name = "Loop signal" )

ht1  
dev.off()









