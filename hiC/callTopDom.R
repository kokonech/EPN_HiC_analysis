source("/icgc/dkfzlsdf/analysis/dktk/Ependymoma/scripts/hiC/TopDom_v0.0.2.R")


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Input file is not provided.", call.=FALSE)
} else {
    inFile = args[1]
    print(paste("Input file:", inFile))
}

resFile= gsub("TopDom_input.txt", "TopDom_result", inFile)

TopDom(matrix.file= inFile, window.size=10,outFile=resFile)


print("Done!")





