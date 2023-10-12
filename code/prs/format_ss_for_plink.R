# Script to readjust SNPs IDs to only include CHR:POS

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide <vilma format> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ss=args[1]
outfile = args[2]


# Read in first data frame
df <- fread(ss)

# Fix IDs
#df <- df %>% separate(ID, c("ID", "tmp1", "tmp2"), "_") %>%
#  select(ID,A1, posterior_ukbb)

colnames(df) <- c("#CHROM", "ID", "POS", "REF", "ALT", "r", "posterior_ukbb","block")
df <- df %>% separate(ID, c("ID", "tmp1", "tmp2"), "_") %>%
  select(ID, ALT, posterior_ukbb)


# Save output
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

