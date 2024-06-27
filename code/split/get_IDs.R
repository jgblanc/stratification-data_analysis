# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<10){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

sex = args[1]
batch = args[2]
age = args[3]
genotyped = args[4]
outGWAS = args[5]
outTest = args[6]
wbs= args[7]
african = args[8]
phenotyped = args[9]
withdraw = args[10]

## Read in all dataframes and join them
df <- fread(sex,  colClasses = 'character')[,1:2]
df <- inner_join(df, fread(batch,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(age,  colClasses = 'character')[,1:2])
df <- inner_join(df, fread(genotyped,  colClasses = 'character')[,1:2])

# Get WBS
df_wbs <- inner_join(df,fread(wbs,  colClasses = 'character')[,1:2])
df_african <- inner_join(df,fread(african,  colClasses = 'character')[,1:2])

# Get test
df_african_test <- df_african %>% sample_n(size = 1000)
df_african_test$POP <- "AFR"
df_wbs_test <- df_wbs %>% sample_n(size = 1000)
df_wbs_test$POP <- "WBS"
df_test <- rbind(df_african_test, df_wbs_test) %>% select("#FID", "IID", african)

# Get GWAS
df_african_gwas <- df_african %>% filter(!IID %in% df_test$IID) %>% sample_n(size = 5000)
df_african_gwas$POP <- "AFR"
df_wbs_gwas <- df_wbs %>% filter(!IID %in% df_test$IID) %>% sample_n(size = 5000)
df_wbs_gwas$POP <- "WBS"
df_gwas <- rbind(df_african_gwas, df_wbs_gwas) %>% select("#FID", "IID", african)


## Select test panel
fwrite(df_test, outTest ,row.names=F,quote=F,sep="\t", col.names = T)

## Select gwas panel
fwrite(df_gwas, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


