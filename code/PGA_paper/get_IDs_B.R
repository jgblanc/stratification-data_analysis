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
gwasSizeWBS = as.numeric(args[9])
gwasSizeOther = as.numeric(args[10])
testSize = as.numeric(args[11])
wbs= args[12]
white = args[13]
all = args[14]

## Read in all dataframes and join them
df <- fread(sex,  colClasses = 'character')[,1:2]
df <- inner_join(df, fread(batch,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(age,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(genotyped,  colClasses = 'character')[,1:2])
df_wbs <- inner_join(df, fread(wbs,  colClasses = 'character')[,1:2])
df_white <- inner_join(df,fread(white,  colClasses = 'character')[,1:2])
df_all <- inner_join(df,fread(all,  colClasses = 'character')[,1:2])
df_other <- df_all %>% filter(!IID %in% df_white$IID)

## Select test panel
df_test <- df_wbs %>% sample_n(testSize) %>% select("#FID", "IID")
fwrite(df_test, outTest ,row.names=F,quote=F,sep="\t", col.names = T)

## Select gwas panel
df_GWAS1 <- df_wbs %>% filter(!IID %in% df_test$IID) %>% sample_n(gwasSizeWBS) %>% select("#FID", "IID")
df_GWAS2 <- df_other %>% filter(!IID %in% df_test$IID) %>% sample_n(gwasSizeOther) %>% select("#FID", "IID")
df_GWAS <- rbind(df_GWAS1, df_GWAS2)
fwrite(df_GWAS, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


