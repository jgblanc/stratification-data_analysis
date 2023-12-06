## Simulate non-genetic phenotype along sharp gradient in WBS UKBB

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript simulate_sharp.R <north coordinates> <outfile> <total shift>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_cords = args[1]
east_cords = args[2]
gwas_ids = args[3]
outfile = args[4]
shift = as.numeric(args[5])
print(paste0("The total shift in phenotype is ", shift))

# Read in IDs and coordinates
ids <- fread(gwas_ids)
north <- fread(north_cords)
east <- fread(east_cords)
df <- inner_join(ids, north)
df <- inner_join(df, east)

# Draw random environment
df$env <- rnorm(n = nrow(df), mean = 0, sd = 1)

# Shift the environment as a function of cordinates
df <- df %>% mutate(pheno = case_when(
  (PlaceOfBirthNorthCord_129 <  190381 & PlaceOfBirthNorthCord_129 >  170381
   & PlaceOfBirthEastCord_130 < 540034 & PlaceOfBirthEastCord_130 > 520034) ~ env + shift ,
  TRUE ~  env))

# Prepare output
dfOut <- df %>% select("#FID", "IID", "pheno")
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








