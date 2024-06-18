## Get Lat/Long Test Vector
## This script uses HGDP metadata and a psam to format test vectors for latitide and longitude

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript get_LatLong_Tvec.R <metadata> <psam> <outfile Lat> <outfile Long>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

md_file = args[1]
psam_file = args[2]
outfile_lat = args[3]
outfile_long = args[4]

# Read in metadata
dfMeta <- fread(md_file)

# Read in Psam
dfPSam <- fread(psam_file)

# Combine files
df <- inner_join(dfPSam, dfMeta, by = c("#IID" = "sample"))

# Only keep Eurasian samples
df <- df %>% filter(region %in% c("CENTRAL_SOUTH_ASIA", "EUROPE",  "MIDDLE_EAST",  "EAST_ASIA"))

# Set up outfiles
dfLat <- df %>% select("#IID", "latitude")
dfLong <- df %>% select("#IID", "longitude")

# Save output
fwrite(dfLat,outfile_lat, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfLong,outfile_long, row.names = F, col.names = T, quote = F, sep = "\t")
