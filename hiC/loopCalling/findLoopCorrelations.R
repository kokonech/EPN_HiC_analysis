library(rtracklayer)
options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)


#######################
# find the correlations

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("ERROR! Input file is not provided.", call.=FALSE)
} else {
  # default output file
    checkFile = args[1]
    print(paste("Input file:", checkFile))
}


resName = gsub(".txt", ".cor_analysis.txt", checkFile)

#checkFile = "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/HiC/hicPro/PFA_EPD210FH/hicpro_result_detailed/fithic_res/loops5KB/fithicoutput_low_bound_10K_res/fithic.spline_pass2.significances.qval_0.05.EPD210FH_enh_contacts.txt"

checkDf <- read.table(checkFile, sep="\t", header =1)

countsDf <- read.table( "/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/rna_sequencing/DESeq/epn_rna_rpkm_ensmbl.package_test.gn.txt", sep="\t", header=1)

enhDf <- read.table("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/results/associateGeneEnhancer/subgroup_specific_stitched_data.package_test.txt", sep="\t", header=1)


#ex 1
# eId <- "chr1:1826712-1827210"
# gName<- "NADK"
#ex 2
# eId = "chr19:13800695-13801193"
# gName = "CACNA1A"

cNames <- paste0(enhDf$chr,":",(enhDf$start + 1),"-",enhDf$end)
enhDf2 <- enhDf[,4:ncol(enhDf)]
rownames(enhDf2) <- cNames

#print("Sample names fit")
#summary(colnames(enhDf2) == colnames(countsDf))
print("Performing analysis...")

id1 = ncol(checkDf) - 1
id2 = ncol(checkDf) 
print(paste("Focus columns:", colnames(checkDf)[id1] , colnames(checkDf)[id2] ))


resFull <- NULL
for (i in 1:nrow(checkDf) ) {
    val1 = checkDf[i,id1]
    val2 = checkDf[i,id2]
    if (val1 %in% rownames(enhDf2)) {
        eId = val1
        gName = strsplit(val2, ",")[[1]][1]
    } else {
        eId = val2
        gName = strsplit(val1, ",")[[1]][1]
    } 
    #print(i)
    enhVals <- enhDf2[eId,]
    exprVals <- countsDf[gName, ]
    if ( sum(exprVals) == 0) {
        next
    }
    res <- cor.test(as.numeric(enhVals), as.numeric(log2(exprVals + 1)), method="spearman")
    # v1
    #resFull <- rbind(resFull, c( eId, gName, res$estimate, res$p.value ) )
    # v2
    resFull <- rbind(resFull, cbind( checkDf[i,1:7], eId, gName, res$estimate, res$p.value  ) )
}

# v1 : with filter
#colnames(resFull) <- c("Enh", "Gene", "Cor", "P-value")
#corFull <- unique(corFull)

# v2 : no filter
colnames(resFull)[8:11] <- c("Enh", "Gene", "Cor", "CorPval")

print("Positive correlation > 0.5")
print(summary(resFull$Cor > 0.5))

# test result
# "binInfo/PFA_connected_correlations.5kbp.180618.txt"

write.table(resFull, resName, sep="\t", quote=F, row.names=F)

print("Done!")


