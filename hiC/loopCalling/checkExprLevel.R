library(rtracklayer)
options(stringsAsFactors = FALSE)
# from local
setwd("/Users/okonechn/work/ependymoma/hiC/clustering/")

epnCounts <- read.table("~/work/ependymoma_src/results/rna_sequencing/DESeq/epn_rna_rpkm_ensmbl.package_test.gn.txt")

mrCounts <- log2(epnCounts + 1)
colnames(mrCounts) <- gsub( "X","",colnames(mrCounts))

annDf <- read.delim("~/work/ependymoma/hiC/2018_12_Ependymoma_Tumors_HiC_QC_Dixon.txt", sep="\t")
rownames(annDf) <- annDf$ID.1

annDf <- annDf[ rownames(annDf) %in% colnames(mrCounts), ]

dataDir = "~/work/ependymoma_src/results/HiC/hicPro/tumors/"
fDir = "/hicpro_result/fithic2_res/loops5KB/fithicoutput_low_bound_10K_res/"

# inspired by 
# https://stats.stackexchange.com/questions/160368/bootstrapping-a-t-test-in-r

# NOTE: t-test might not be suitable due to distribution!!

# other sources
# https://stats.stackexchange.com/questions/175321/independent-samples-t-test-with-unequal-sample-sizes
# https://stats.stackexchange.com/questions/92542/how-to-perform-a-bootstrap-test-to-compare-the-means-of-two-samples/92553#92553


bootstrap.t.test <- function(testExpr, targExpr, B= 1000) {

    t.vect <- vector(length=B)
    p.vect <- vector(length=B)

    for(i in 1:B){
        #print(i)
        boot.targ <- sample(targExpr, size=length(testExpr), replace=T)
        ttest <- t.test(testExpr, boot.targ, alternative="greater")
        t.vect[i] <- ttest$statistic
        p.vect[i] <- ttest$p.value
    }
    
    mean(p.vect)
}

set.seed(1315)
fullRes <- NULL
statRes <- NULL

for (sId in rownames(annDf)) {
    print(sId)
    fPath = paste0(dataDir,sId,fDir, "FitHiC.spline_pass2.res5000.significances.qval_0.05.gene_contacts_closest.txt")
    contacts <- read.table(fPath, header=TRUE, sep="\t")

    genes1 = unlist(strsplit( contacts$target1, ","))
    genes2 = unlist(strsplit( contacts$target2, ","))

    targetGenes <- unique( c(genes1, genes2))
    #print(summary(targetGenes %in% rownames(mrCounts)))

    # group-specifc enhancers - genes
    fPath2 = paste0(dataDir,sId,fDir, "FitHiC.spline_pass2.res5000.significances.qval_0.05.enh_contacts_closest.txt")
    enhContacts <- read.table(fPath2, header=TRUE, sep="\t")
    
    enhGenes1 <- unlist(strsplit( enhContacts$target1, ","))
    enhGenes2 <- unlist(strsplit( enhContacts$target2, ","))
    
    enhGenes <- unique(c(enhGenes1,enhGenes2))
    enhGenes <- enhGenes[ enhGenes %in% rownames(mrCounts) ]
    
    # filter out enh genes
    targetGenes <- setdiff(targetGenes,enhGenes)
    
    otherGenes <- rownames(mrCounts)[ !( rownames(mrCounts) %in% c(targetGenes,enhGenes) ) ]

    t1 <- data.frame( expr = mrCounts[targetGenes, sId], type= "Genes", sample=sId)
    rownames(t1) <- targetGenes

    t2 <- data.frame( expr = mrCounts[enhGenes, sId], type= "Enhancers", sample=sId)
    rownames(t2) <- enhGenes
    
    t3 <- data.frame( expr = mrCounts[otherGenes, sId], type= "NoLoops", sample=sId)
    rownames(t3) <- otherGenes
    
    res1 <- bootstrap.t.test(t1$expr, t3$expr)
    res2 <- bootstrap.t.test(t2$expr, t3$expr)
    statRes <- rbind(statRes,c( nrow(t1), res1, nrow(t2), res2,nrow(t3) ))
    
    fullRes <- rbind(fullRes, t1, t2, t3)
    print(paste(sId,":", nrow(t1),nrow(t2),nrow(t3)))

}

# save stat.res
colnames(statRes) <- c("TargGenes","Pval_genes", "TargEnh", "Pval_enh",  "Other")
rownames(statRes) <- rownames(annDf)
write.table(statRes, "genes_in_loops.bootstrap_ttest_pvals.txt",sep="\t",quote=F)

# control order
fullRes$sample <- factor(fullRes$sample,levels = c("7EP18","9EP1","9EP9", "7EP41","11EP22","4EP53"))


my_comparisons <- list( c("Genes", "NoLoops"), c("Enhancers", "NoLoops") )
    
library(ggplot2)
pdf("genes_in_loops_closest_expr_boxplot_extended.group_spec_enh.220520.pdf", height=5, width=10)
ggplot(fullRes, aes(x=sample, y=expr, fill=type)) +
    ylim(0,9.5) +
    scale_fill_manual(values=c( "yellow","cadetblue3", "grey")) +
    theme(axis.text.x = element_text( color=c(rep("orange",3),rep("red",3)),size=10,face="bold" ) ) +
    labs(title = "Gene expression: loops effect",x = "Sample", y="log2 RPKM") +
    geom_boxplot(outlier.shape = NA)
dev.off()

fullRes$expr <- round(fullRes$expr,2)
write.table(fullRes, "genes_in_loops_closest_expr_boxplot_extended.group_spec_enh.output.txt")

