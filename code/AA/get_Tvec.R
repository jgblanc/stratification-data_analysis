## Get North/East Test Vector
## This script uses UKBB North or East birth locations and outputs test vector for a subset of individuals

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript get_Tvec.R <metadata> <id file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_file = args[1]
east_file = args[2]
id_file = args[3]
outfile_north = args[4]
outfile_east = args[5]
outfile_sharp = args[6]

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
dfNorth <- dfNorth %>% select("#FID","IID", "north")
dfEast <- dfEast %>% select("#FID","IID", "east")

# Get sharp Tvec
df <- inner_join(dfNorth, dfEast)
df <- df %>% mutate(sharp = case_when(
  (north <  190381 & north >  170381
   & east < 540034 & east > 520034) ~ 1 ,
  TRUE ~  0))

df$sharp <- df$sharp - mean(df$sharp)

# Set up outfiles
dfSharp <- df %>% select("#FID","IID", "sharp")

# Save output
fwrite(dfNorth,outfile_north, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfEast,outfile_east, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfShar,outfile_sharp, row.names = F, col.names = T, quote = F, sep = "\t")
