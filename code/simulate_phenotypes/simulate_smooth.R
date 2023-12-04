## Simulate non-genetic phenotype along north south gradient in WBS UKBB

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simulate_smooth.R <north coordinates> <outfile> <total shift>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_cords = args[1]
gwas_ids = args[2]
outfile = args[3]
shift = as.numeric(args[4])
print(paste0("The total shift in phenotype is ", shift))

# Read in IDs and coordinates
ids <- fread(gwas_ids)
north <- fread(north_cords)
df <- inner_join(ids, north)


# Draw random environment
df$env <- rnorm(n = nrow(df), mean = 0, sd = 1)

# Shift the environment as a function of cordinates
maxN <- max(df$PlaceOfBirthNorthCord_129)
tmp <- shift / maxN
df$pheno <- df$env + (tmp * df$PlaceOfBirthNorthCord_129)

# Prepare output
dfOut <- df %>% select("#FID", "IID", "pheno")
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








