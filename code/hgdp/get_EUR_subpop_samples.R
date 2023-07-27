## Get EUR samples
## This script subsets the list of HGDP individuals to only those in europe

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript get_EUR_subpop_samples.R <population meta data> <output_prefix>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

popfile = args[1]
outpre = args[2]


# Read in metadata
pops <- fread(popfile)

# Subset to Europeans
eur <- pops %>% filter(region == "EUROPE")

# Get list of subpops and save list of sample IDs
subpops <- unique(eur$population)
for (i in 1:length(subpops)) {

  # Filter IDs
  df <- eur %>% filter(population == subpops[i]) %>% select(sample)
  df$SEX <- rep(NA, nrow(df))
  colnames(df) <- c("#IID", "SEX")

  # Write output
  outfile <- paste0(outpre, subpops[i], ".txt")
  print(outfile)
  fwrite(df, outfile,row.names = F, col.names = T, quote = F, sep = "\t")

}


