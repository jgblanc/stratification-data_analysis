## Overlapping SNPs
## This script takes a .afreq file from UKBB, subsets it to MAF > 1%, and intersects it with SDS snps (already limited to MAF > 5% in UK10K)
## It then flips the SDS value if the derived allele is the reference - this ensures that the SDS value is for the ALT allele which matches the rest of the pipeline


args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript flip_overlapping_snps.R <ukbb.freq> <SDS> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ukbb_file = args[1]
sds_file = args[2]
outfile_snps = args[3]

# Read in UKBB freq file
ukbb <- fread(ukbb_file)
ukbb <- ukbb %>% separate(ID, c("chr", "POS"), remove = FALSE) %>% filter(ukbb$ALT_FREQS > 0.01 & ukbb$ALT_FREQS < 0.99)

# Read in SDS
sds <- fread(sds_file)
sds$POS <- as.character(sds$POS)

# Join files by chromosome and position
df <- inner_join(ukbb, sds, by = c("#CHROM" = "CHR", "POS" = "POS"))

# Get rid of rows where any alleles don't match
df_filter <- df %>% filter((DA == REF & AA == ALT) | (DA == ALT | AA == REF))

# If derived and alternate allele don't match, flip SDS
df_flipped <- df_filter %>% mutate(SDS = case_when(DA == REF ~ (-1 * SDS), DA == ALT ~ SDS))

# Select correct columns
df_out <- df_flipped %>% select("#CHROM", "ID.x", "REF", "ALT", "SDS")
colnames(df_out) <- c("#CHROM", "ID", "REF", "ALT", "r")

# Format for overlapping snps
df_snps <- df_out %>% select(ID)

# Save output
fwrite(df_snps,outfile_snps, row.names = F, col.names = T, quote = F, sep = "\t")
