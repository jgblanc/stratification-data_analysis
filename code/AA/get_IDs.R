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


## Read in all dataframes and join them
df <- fread(sex)[,1:2]
df <- inner_join(df, fread(batch)[,1:2])
df <- inner_join(df, fread(age)[,1:2])
df <- inner_join(df, fread(genotyped)[,1:2])

## Get WBS for test panel
dfWBS <- inner_join(df, fread(wbs)[,1:2])
dfWBS <- inner_join(dfWBS, fread(north)[,1:2])
dfWBS <- inner_join(dfWBS, fread(east)[,1:2])
print("The WBS total size is ", nrow(dfWBS))

## Get White
dfWhite <- inner_join(df, fread(white)[,1:2])
dfWhite <- dfWhite %>% filter(!IID %in% dfWBS$IID)
print("The White total size is ", nrow(dfWhite))

## Other ancestries
dfOther <- df %>% filter(!IID %in% dfWBS$IID) %>% filter(!IID %in% dfWhiteS$IID)  %>% select("#FID", "IID")
print("The other total size is ", nrow(dfOther))

## Select test panel
df_test <- dfWBS %>% sample_n(testSize) %>% select("#FID", "IID")
fwrite(df_test, outTest ,row.names=F,quote=F,sep="\t", col.names = T)

## Select gwas panel

# Get 50K "white"
dfGWAS_white <- dfWhite %>% sample_n(50000) %>% select("#FID", "IID")

# Calculate the number of other ancestries we need
nOther <- gwasSize - 50000
print(paste0("The number of other ancestries is ", nOther))

# Select other
dfGWAS_other <- dfOther %>% sample_n(nOther) %>% select("#FID", "IID")

# Combine GWAS panel
df_GWAS <- rbind(df_GWAS_white, df_GWAS_other)

# Save output
fwrite(df_GWAS, outGWAS ,row.names=F,quote=F,sep="\t", col.names = T)


