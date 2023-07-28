## Get pairwise EUR_AF_contrasts
# This script get chromosome specific allele frequency contrasts (r) for all pairwise HGDP EUR pops

args=commandArgs(TRUE)

if(length(args)<1){stop("Rscript get_pairwise_EUR_AF_contrasts.R <number of inputs> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

spnum = as.numeric(args[1])
out_pre = args[2+spnum]

# Get all pairwise combos
pw_comb <- t(combn(seq(1, spnum),2))

for (i in 1:nrow(pw_comb)) {

  # Get pairwise IDs
  id1 <- pw_comb[i,1]
  id2 <- pw_comb[i,2]

  # Read in dataframes
  df1 <- fread(args[id1+1])
  df2 <- fread(args[id2+1])

  # Compute r
  r <- df1$ALT_FREQS - df1$ALT_FREQS

  # Set up output dataframe
  out <- df1[,1:4]
  out$r <- r

  # Get subpop names
  nm1 <- strsplit(strsplit(args[id1+1], "/")[[1]][5], "_")[[1]][1]
  nm2 <- strsplit(strsplit(args[id2+1], "/")[[1]][5], "_")[[1]][1]

  # Get chromosome number
  chr <- strsplit(strsplit(strsplit(x, "/")[[1]][5], "_")[[1]][2], ".afreq")[[1]][1]


  # Save output
  outfile <- paste0(out_pre, nm1, "-", nm2, "_", chr, ".rvec")
  print(outfile)
  fwrite(df, outfile,row.names = F, col.names = T, quote = F, sep = "\t")

}




