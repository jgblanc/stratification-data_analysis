## Get North Test Vector
## This script uses UKBB North birth locations and outputs test vector for a subset of individuals

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_Tvec_North.R <metadata> <id file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_file = args[1]
id_file = args[2]
outfile_north = args[3]

# Read in metadata
dfNorth <- fread(north_file)

# Read in TP IDs
dfID <- fread(id_file)

# Combine files
dfNorth <- inner_join(dfNorth, dfID)
colnames(dfNorth)[3] <- "smooth"
print(head(dfNorth))


# Set up outfiles
dfNorth <- dfNorth %>% select("#FID","IID", "smooth")


# Save output
fwrite(dfNorth,outfile_north, row.names = F, col.names = T, quote = F, sep = "\t")
