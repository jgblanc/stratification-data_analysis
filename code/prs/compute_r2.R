## This script takes in PRSs and Phenotypes and compues r and r^2


args=commandArgs(TRUE)

if(length(args)<3){stop("Provide <phenotype> <prs> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))


pheno_file = args[1]
prs_file = args[2]
out_file = args[3]

# Read in phenotypes
phenos <- fread(pheno_file)

# Read in PRS
prs <- fread(prs_file)

# Join data
df <- inner_join(prs, pheno, by = c("#IID"= "IID"))

# Compute r
r <- cor(df$SCORE1_SUM, df$Height)

# Compute r^2
r2 <- r^2

# Set up output table
dfOut <- matrix(NA, ncol = 2, nrow = 1)
colnames(dfOut) <- c("r", "r2")

# Save output
fwrite(dfOut, out_file, row.names = F, col.names = F, quote = F, sep = "\t")

