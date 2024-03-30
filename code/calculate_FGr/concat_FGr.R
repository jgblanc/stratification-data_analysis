## Concatenate projections
## This script reads in all the projected test vectors for each chromosome and adds them all together, leaving the focal chromosome out, and scales the final covariate

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript conccat_FGr.R <prefix to Tm chromosomes> <chromosome number> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

outfile = args[1]
print(length(args))

# Add Tm for each chromosome to each other
df <- fread(args[2])
for (i in 3:length(args)) {

  new <- fread(args[i])
  df$FGr <- df$FGr + new$FGr

}

df$FGr <- scale(df$FGr)

# Save output
fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







