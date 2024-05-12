## Overlapping SNPs
## This script takes a .afreq file from UKBB, subsets it to MAF > 1%, and intersects it with SDS snps (already limited to MAF > 5% in UK10K)
## It then flips the SDS value if the derived allele is the reference - this ensures that the SDS value is for the ALT allele which matches the rest of the pipeline


args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript flip_overlapping_snps.R <ukbb.freq> <SDS> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ukbb_file_wbs = args[1]
ukbb_file_all = args[2]
sds_file = args[3]
outpath = args[4]

# Read in UKBB freq file
ukbb_wbs <- fread(ukbb_file_wbs)
ukbb_wbs <- ukbb_wbs %>% separate(ID, c("tmp", "POS"), remove = FALSE) %>% filter(ukbb_wbs$ALT_FREQS > 0.01 & ukbb_wbs$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")
ukbb_all <- fread(ukbb_file_all)
ukbb_all <- ukbb_all %>% separate(ID, c("tmp", "POS"), remove = FALSE) %>% filter(ukbb_all$ALT_FREQS > 0.01 & ukbb_all$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")
ukbb <- full_join(ukbb_wbs, ukbb_all)

# Read in contrasts
sds <- fread(sds_file)
sds$POS <- as.character(sds$POS)

# Join files by chromosome and position
df <- inner_join(ukbb, sds, by = c("#CHROM" = "CHR", "POS" = "POS", "REF" = "REF", "ALT" = "ALT"))

# Change allele frequency to alternate allele freq by taking 1 -
df_flipped <- df_filter %>% mutate(SDS = case_when(DA == REF ~ (-1 * SDS), DA == ALT ~ SDS))

# Select correct columns
df_out <- df_flipped %>% select("#CHROM", "ID.x", "REF", "ALT", "SDS")
colnames(df_out) <- c("#CHROM", "ID", "REF", "ALT", "r")

# Save output
fwrite(df_out,outfile_r, row.names = F, col.names = T, quote = F, sep = "\t")


