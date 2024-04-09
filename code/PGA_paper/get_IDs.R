# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sex = args[1]
batch = args[2]
north = args[3]
east = args[4]
age = args[5]
genotyped = args[6]
outGWAS = args[7]
outTest = args[8]
gwasSize = as.numeric(args[9])
testSize = as.numeric(args[10])
wbs= args[11]
white = args[12]
phenotyped = args[13]

## Read in all dataframes and join them
df <- fread(sex,  colClasses = 'character')[,1:2]
df <- inner_join(df, fread(batch,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(age,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(genotyped,  colClasses = 'character')[,1:2])

# Get GWAS
df_white <- inner_join(df,fread(white,  colClasses = 'character')[,1:2])
df_white <- inner_join(df_white,fread(phenotyped,  colClasses = 'character')[,1:2])

# Get test
df_wbs <- inner_join(df, fread(wbs,  colClasses = 'character')[,1:2])
df_wbs <- inner_join(df_wbs,fread(north,  colClasses = 'character')[,1:2])
df_wbs <- inner_join(df_wbs,fread(east,  colClasses = 'character')[,1:2])


## Select test panel
df_test <- df_wbs %>% sample_n(testSize) %>% select("#FID", "IID")
fwrite(df_test, outTest ,row.names=F,quote=F,sep="\t", col.names = T)

## Select gwas panel
df_GWAS <- df_white %>% filter(!IID %in% df_test$IID) %>% sample_n(gwasSize) %>% select("#FID", "IID")
fwrite(df_GWAS, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


