## Get Sharp Test Vector
## This script uses UKBB North and East birth locations and outputs test vector for a subset of individuals

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript get_Tvec_North.R <metadata> <id file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

north_file = args[1]
east_file = args[2]
id_file = args[3]
outfile = args[4]

# Read in metadata
dfNorth <- fread(north_file)
dfEast <- fread(east_file)

# Read in TP IDs
dfID <- fread(id_file)

# Combine files
dfNorth <- inner_join(dfNorth, dfID)
df <- inner_join(dfNorth, dfEast)

# Mark if in coordinate range
df <- df %>% mutate(sharp = case_when(
  (PlaceOfBirthNorthCord_129 <  190381 & PlaceOfBirthNorthCord_129 >  170381
   & PlaceOfBirthEastCord_130 < 540034 & PlaceOfBirthEastCord_130 > 520034) ~ 1 ,
  TRUE ~  0))

# Set up outfiles
df <- df %>% select("#FID","IID", "sharp")


# Save output
fwrite(df,outfile, row.names = F, col.names = T, quote = F, sep = "\t")
