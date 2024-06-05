# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sex = args[1]
batch = args[2]
age = args[3]
genotyped = args[4]
phenotyped = args[5]
withdraw = args[6]
outGWAS = args[7]

## Read in all dataframes and join them
df <- fread(sex,  colClasses = 'character')[,1:2]
df <- inner_join(df, fread(batch,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(age,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(genotyped,  colClasses = 'character')[,1:2])

# Get withdrawn individuals
df_wd <- fread(withdraw, header = F)
colnames(df_wd) <- "IID"
print(nrow(df_wd))


## Select gwas panel
df_GWAS <- df %>% filter(!IID %in% df_wd$IID) %>% select("#FID", "IID")
fwrite(df_GWAS, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


