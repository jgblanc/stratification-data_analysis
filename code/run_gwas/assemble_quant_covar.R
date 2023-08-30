# Script to format quantitative covariates

args=commandArgs(TRUE)

if(length(args)<5){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

age = args[1]
fgr = args[2]
pc = args[3]
pcNum = as.numeric(args[4])
outfile = args[5]

# Read in dataframes
dfAge <- fread(age)
dfPc <- fread(pc)
dfFgr <- fread(fgr)

# Subset to correct number of PCs
dfPc <- dfPc[,1:(pcNum + 2)]
print(ncol(dfPc))

# Join files
df <- inner_join(dfAge, dfPc, dfFgr, by = c("#FID", "IID"))

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

