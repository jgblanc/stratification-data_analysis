# This script concatenates results from computing the error in FGr

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_results.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 8)

for (i in 2:length(args)) {

  # Get results file name
  filename = args[i]

  # Extract parameters
  tmp <- strsplit(filename, "/")[[1]]
  rep <- tmp[4]
  nsnp <- as.numeric(strsplit(tmp[5], "L-")[[1]][2])
  phenotype <- tmp[6]
  covar <- strsplit(tmp[8], "_")[[1]][1]
  contrast <- strsplit(strsplit(tmp[8], "_")[[1]][2], "-PC")[[1]][1]
  pcNum <- strsplit(strsplit(strsplit(tmp[8], "_")[[1]][2], "-PC")[[1]][2], "all.results")[[1]][1]

  # Read in results
  dfTmp <- fread(filename)
  q <- dfTmp[1,1]
  pval <- dfTmp[1,2]

  # Save results in table
  x <- c(q,pval, rep, nsnp,contrast, phenotype, covar, pcNum)
  dfOut <- rbind(dfOut, x)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
colnames(dfOut) <- c("q", "pval","rep", "L", "contrast", "phenotype","covar", "pc")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




