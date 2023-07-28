## Get pairwise EUR_AF_contrasts
# This script get chromosome specific allele frequency contrasts (r) for all pairwise HGDP EUR pops

args=commandArgs(TRUE)

if(length(args)<1){stop("Rscript get_pairwise_EUR_AF_contrasts.R <input>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

infiles = args[1]
print(infiles)




