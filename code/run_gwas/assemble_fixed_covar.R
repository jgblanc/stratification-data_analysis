# Script to format categorical covariates

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

sex = args[1]
age = args[2]
batch = args[3]
outfile = args[4]

# Read in phenotypes
dfSex <- fread(sex)
head(dfSex)
dfAge <- fread(age)
head(dfAge)
dfBatch <- fread(batch)
head(dfBatch)

# Join files
df <- inner_join(dfSex, dfAge, dfBatch, by = c("#FID", "IID"))

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

