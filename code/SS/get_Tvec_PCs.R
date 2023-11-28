## Get North/East Test Vector
## This script uses UKBB North or East birth locations and outputs test vector for a subset of individuals

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript get_Tvec.R <metadata> <id file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

pc_file = args[1]
prefix_out = args[2]

# Read in metadata
dfPCA <- fread(pc_file)


# Subset columns
for (i in 1:5) {

  # Ger correct column
  n <- as.numeric(args[2+i])
  selected_columns <- c(1, 2, n+2)
  df_selected <- dfPCA[, selected_columns, with = FALSE]
  name <- paste0("pc", n)
  colnames(df_selected)[3] <- name

  # Save output
  outfile <- paste0(prefix_out, name, ".txt")
  print(outfile)
  fwrite(df_selected,outfile, row.names = F, col.names = T, quote = F, sep = "\t")
}


