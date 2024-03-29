# This script concatenates results from selection scan into plotting format

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_results.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 12)

for (i in 2:length(args)) {

  # Get results file name
  filename = args[i]

  # Extract phenotype
  phenotype <- strsplit(filename, "/")[[1]][6]

  # Extract which covariate was used
  tmp <- strsplit(filename, "/")[[1]][8]
  covar <- strsplit(tmp, "_")[[1]][1]
  contrast <- strsplit(strsplit(tmp, "_")[[1]][2], ".all")[[1]][1]
  if (covar == "FGr-LOCO") {
    covar <- TRUE
  } else if (covar == "no-FGr") {
    covar <- FALSE
  } else {
    covar <- NA
  }

  # Extract number of PCs used
  pc <- regmatches(tmp, gregexpr("([0-9]|[1-9][0-9]|1[0-9]{2}|2[0-4][0-9]|25[0-5])", tmp))[[1]]

  # Test replicate
  tmp2 <- strsplit(filename, "/")[[1]][4]
  test <- tmp2
  testSize <- as.numeric(strsplit(tmp2, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)[[1]][2])

  # Gwas size
  tmp3 <- strsplit(filename, "/")[[1]][5]
  gwas <- tmp3
  gwasSize <- as.numeric(strsplit(tmp3, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)[[1]][2])

  # Read in results
  dfTmp <- fread(filename)
  q <- dfTmp[1,1]
  p <- dfTmp[1,2]
  q_NS <- dfTmp[1,3]
  p_NS <- dfTmp[1,4]

  # Save results in table
  x <- c(q, p, q_NS, p_NS, phenotype, covar, pc, contrast, test, testSize,gwas, gwasSize)
  dfOut <- rbind(dfOut, x)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
colnames(dfOut) <- c("q", "p", "q_NoSign", "p_NoSign", "phenotype", "covar", "pc", "contrast", "test", "testSize","gwas", "gwasSize")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")




