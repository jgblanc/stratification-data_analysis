# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

infile = args[1]
nsnp = strsplit(args[2], "L-")[[1]][2]
outfile = args[3]

if (nsnp != "all") {
  nsnp = as.numeric(nsnp)
  df <- fread(infile)[,2:3]
  df <- df %>% sample_n(nsnp)
  fwrite(df, outfile ,row.names=F,quote=F,sep="\t", col.names = T)
} else {
  df <- fread(infile)[,2:3]
  fwrite(df, outfile ,row.names=F,quote=F,sep="\t", col.names = T)
}



