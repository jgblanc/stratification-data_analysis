# This script concatenates results from computing the error in FGr

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_results.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 4)

for (i in 2:length(args)) {

  # Get results file name
  filename = args[i]

  # Extract parameters
  tmp <- strsplit(filename, "/")[[1]]
  test <- tmp[4]
  gwas <- tmp[5]
  contrast <- strsplit(strsplit(tmp[7], "_")[[1]][2], ".txt")[[1]][1]

  # Read in results
  dfTmp <- fread(filename)
  varE <- dfTmp[1,1]

  # Save results in table
  x <- c(varE, test, gwas,contrast)
  dfOut <- rbind(dfOut, x)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
colnames(dfOut) <- c("varE", "subdataset", "gwas",  "contrast")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




