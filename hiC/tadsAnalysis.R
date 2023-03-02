#install.packages("UpSetR")
library(UpSetR)
library(rtracklayer)
options(stringsAsFactors=FALSE)
setwd("/Users/okonechn/work/ependymoma/hiC/tumorTADs/")

# use Upset

tadStatus <- read.table("/Users/okonechn/work/ependymoma_src/results/HiC/hicPro/commonTADs/inputForUpset.260919.txt",header = T)
colnames(tadStatus) <- gsub("X","",colnames(tadStatus))

cols <- rep("black",57)
# for 6 samples
cols[c(2,7,11,16)] <- "orange"
cols[c(3,4,8,14)] <- "red"


pdf("EPN_HiC_comon_TADs.050220.pdf", width=8,height=6)
upset(tadStatus, #nsets=6, 
      sets = colnames(tadStatus)[2:ncol(tadStatus)], 
      #scale.intersections="log2",
      nintersects = 20,
      # mb.ratio = c(0.8, 0.2), # vertical
      main.bar.color = cols,
      keep.order = TRUE,
      cutoff = 500,
      sets.bar.color = "#56B4E9", 
      order.by = "freq")
dev.off()

# TADs with AstroCb

tadStatus <- read.table("inputForUpset_with_AstroCb.150420.txt",header = T)
colnames(tadStatus) <- gsub("X","",colnames(tadStatus))

selPFA <- as.integer(rowSums(tadStatus[,2:4]) > 1)
selPFB <- as.integer(rowSums(tadStatus[,5:7]) > 1)
selCB <- as.integer(rowSums(tadStatus[,8:9]) > 1)

common <- sum(selPFA & selPFB & selCB)
sel12 <- sum(selPFA & selPFB)
sel13 <- sum(selPFA & selCB)
sel23 <- sum(selPFB & selCB)

pdf("common_TADs.PFA_RELA_AstroCB.160420.pdf",width=6,height = 6)
grid.newpage()
draw.triple.venn(area1 = sum(selPFA), area2 = sum(selPFB), area3 = sum(selCB),
                 n12 = sel12, n23 = sel23, n13 = sel13, n123 = common, 
                 category = c("PFA", "RELA", "AstroCb"), lty = "blank", 
                 fill = c("orange", "red", "blue"))
#require(gridExtra)
#grid.arrange(gTree(children=venn.plot), top="Common TADs")
dev.off()


## TAds similarity

library(corrplot)

corData <- read.table("TADs_common_proprotion.with_AstroCb.dist_50Kbp.txt",sep="\t",header=1)
colnames(corData) <- gsub("X","",colnames(corData))
# 03.08.20 exclude BT214
corData <- corData[c(1:6,8:9),c(1:6,8:9) ]


#M <- cor(mtcars)
col3 <- colorRampPalette(c("red", "white", "green4" ))
labCols <- c(rep("orange",3),rep("red",3),rep("blue",2))
pdf("TADs_common_proprotion.030820.pdf",width=6,height=6)
corrplot(as.matrix(corData), method = "circle",
         tl.col = labCols, cl.lim = c(0,1),col = col3(200))
dev.off()



library(gplots)
pdf("heatmap2_cor_plot.pdf",width=8,height=8)
heatmap.2(t(as.matrix(corData)),
          scale = "none", dendrogram = "row",
          hclustfun =function(x) hclust(x,method = "ward.D2"))
dev.off()



