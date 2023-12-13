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


## Read in all dataframes and join them
df <- fread(sex)[,1:2]
df <- inner_join(df, fread(batch)[,1:2])
df <- inner_join(df, fread(age)[,1:2])
df <- inner_join(df, fread(genotyped)[,1:2])

## Get WBS for test panel
dfWBS <- inner_join(df, fread(wbs)[,1:2])
dfWBS <- inner_join(dfWBS, fread(north)[,1:2])
dfWBS <- inner_join(dfWBS, fread(east)[,1:2])

## Select test panel
df_test <- dfWBS %>% sample_n(testSize) %>% select("#FID", "IID")
fwrite(df_test, outTest ,row.names=F,quote=F,sep="\t", col.names = T)

## Select gwas panel
df_GWAS_other <- df %>% filter(!IID %in% dfWBS$IID) %>% select("#FID", "IID")
nOther <- nrow(df_GWAS_other)
print(paste0("The number of other ancestries is ", nOther))
df_GWAS_WBS <- dfWBS  %>% filter(!IID %in% df_test$IID) %>% sample_n(gwasSize - nOther) %>% select("#FID", "IID")
df_GWAS <- rbind(df_GWAS_other, df_GWAS_WBS)
print(paste0("The total number in the GWAS is ", nrow(df_GWAS)))

fwrite(df_GWAS, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


