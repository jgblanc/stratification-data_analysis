## Get North/East Test Vector
## This script uses UKBB North or East birth locations and outputs test vector for a subset of individuals

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript get_Tvec.R <metadata> <id file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_file = args[1]
east_file = args[2]
id_file = args[3]
outfile_north = args[4]
outfile_east = args[5]

# Read in metadata
dfNorth <- fread(north_file)
dfEast <- fread(east_file)

# Read in TP IDs
dfID <- fread(id_file)

# Combine files
dfNorth <- inner_join(dfNorth, dfID)
colnames(dfNorth)[3] <- "north"
print(head(dfNorth))

dfEast <- inner_join(dfEast, dfID)
colnames(dfEast)[3] <- "east"
print(head(dfEast))

# Set up outfiles
dfNorth <- dfNorth %>% select("IID", "north")
dfEast <- df %>% select("IID", "east")


# Save output
fwrite(dfNorth,outfile_north, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfEast,outfile_east, row.names = F, col.names = T, quote = F, sep = "\t")
