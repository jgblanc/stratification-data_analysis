# Script formats standing height phenotype file

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide path to phenotype file")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

pheno_file = args[1]
out_file = args[2]

# Read in phenotypes
df <- fread(pheno_file)
#df <-fread("~/Downloads/ethnic_grouping_22006.txt")

# Select ID and initial visit heights and rename columns
df <- df[,1:2]
colnames(df) <- c("IID", "Pheno")

# Remove individuals with no initial value for height
df <- df %>% drop_na(Pheno)
df$FID <- df$IID
df <- df[,c("FID", "IID", "Pheno")]
colnames(df) <- c("#FID", "IID", "Pheno")

# Get name of phenotype
name <- tools::file_path_sans_ext("ABCD.csv")
colnames(df)[3] <- name

# Save file
fwrite(df, out_file,col.names=T,row.names=F,quote=F,sep="\t")
