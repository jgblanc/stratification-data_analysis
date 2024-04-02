# Script to format quantitative covariates when not including FGr

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

age = args[1]
outfile = args[2]

# Read in dataframes
df <- fread(age)


# Drop rows with NA
df <- drop_na(df)

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

