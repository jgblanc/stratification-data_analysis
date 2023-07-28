## Get pairwise EUR_AF_contrasts
# This script get chromosome specific allele frequency contrasts (r) for all pairwise HGDP EUR pops

args=commandArgs(TRUE)

if(length(args)<1){stop("Rscript get_pairwise_EUR_AF_contrasts.R <number of inputs> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

spnum = args[1]
print(spnum)

for (i in 2:(spnum+1)) {

  print(args[i])

}



