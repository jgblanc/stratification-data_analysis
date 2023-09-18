# Script to format summary statistics to match vilma
## Plink2 uses the the alt allele as the "coded" allele so need effect size for alt
## Vilma uses CHR:POS_REF_ALT so need to convert ss IDs

## fastGWA ALT is the effect allele

args=commandArgs(TRUE)

if(length(args)<3){stop("Provide <fastGWA summary stats> <plink2 pvar file> <SNP IDs from vilma ID matrix> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ss_file = args[1]
snp_file = args[2]
out_file = args[3]


# Read inn summary statistics
dfSS <- fread(ss_file)

# Read in SNP IDs used in Vilma
dfSNP <- fread(snp_file)

# Combine dfs
df <- inner_join(dfSS, dfSNP, by = c("CHR" = "CHROM", "POS" = "BP", "A1" = "A1", "A2" = "A2") )

# Save output
fwrite(df, out_file,col.names=T,row.names=F,quote=F,sep="\t")
