# Script to format quantitative covariates

args=commandArgs(TRUE)

if(length(args)<5){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

age = args[1]
fgr = args[2]
outfile = args[3]

# Read in dataframes
dfAge <- fread(age)
dfFgr <- fread(fgr)
print(head(dfFgr))
colnames(dfFgr) <- c("#FID" ,"IID","FGr")
dfFgr <- dfFgr %>% select("#FID" ,"IID", "FGr")
print(head(dfFgr))


# Join files
df <- inner_join(dfAge, dfFgr, by = c("#FID", "IID"))

# Drop rows with NA
df <- drop_na(df)
print(head(df))

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

