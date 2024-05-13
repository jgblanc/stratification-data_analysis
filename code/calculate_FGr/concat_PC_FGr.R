# This script concatenates results from computing the error in FGr

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_results.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 10)
colnames(dfOut) <- c("PC","r2",   "rho",  "lc", "up", "test", "gwas", "contrast", "chrtype_gwas", "chrtype_test")
#print(dfOut)

for (i in 2:length(args)) {

  # Get results file name
  filename = args[i]

  # Extract parameters
  tmp <- strsplit(filename, "/")[[1]]
  test <- tmp[4]
  gwas <- tmp[5]
  contrast <- strsplit(tmp[8], "_PCs_FGr.txt")[[1]][1]
  chrtype_gwas <- tmp[6]
  chrtype_test <- tmp[7]

  # Read in results
  dfTmp <- fread(filename,header=TRUE)
  dfTmp$test <- test
  dfTmp$gwas <- gwas
  dfTmp$contrast <- contrast
  dfTmp$chrtype_gwas <- chrtype_gwas
  dfTmp$chrtype_test <- chrtype_test
  #print(dfTmp)

  # Save results in table
  dfOut <- rbind(dfOut, dfTmp)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
#colnames(dfOut) <- c("PC",	"r2",	"rho", 	"lc", "up", "subdataset", "gwas", "constrast", "chrtype_gwas", "chrtype_test")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




