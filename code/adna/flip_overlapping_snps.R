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
chr_num = args[5]
print(chr_num)


# Read in UKBB freq file
ukbb_wbs <- fread(ukbb_file_wbs)
ukbb_wbs <- ukbb_wbs %>% separate(ID, c("tmp", "POS"), remove = FALSE) %>% filter(ukbb_wbs$ALT_FREQS > 0.01 & ukbb_wbs$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")
ukbb_all <- fread(ukbb_file_all)
ukbb_all <- ukbb_all %>% separate(ID, c("tmp", "POS"), remove = FALSE) %>% filter(ukbb_all$ALT_FREQS > 0.01 & ukbb_all$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")
ukbb <- full_join(ukbb_wbs, ukbb_all)
colnames(ukbb)[1] <- "CHROM"
ukbb <- ukbb %>% separate(ID, c("tmp", "POS"), remove = FALSE)
ukbb$POS <- as.numeric(ukbb$POS)

# Read in contrasts
sds <- fread(sds_file)
head(str(sds))
sds <- sds %>% filter(CHROM == as.integer(chr_num))

# Join files by chromosome and position
df <- inner_join(ukbb, sds)
head(df)

# Change allele frequency to alternate allele freq by taking 1 -
df_flipped <- ss %>% mutate(EUROPEAN_NEOLITHIC = 1 - EUROPEAN_NEOLITHIC, BRONZE_AGE = 1 - BRONZE_AGE, HISTORICAL = 1 - HISTORICAL,
                            ANATOLIA_NEOLITHIC = 1 -ANATOLIA_NEOLITHIC, MESOLITHIC = 1 - MESOLITHIC, STEPPE = 1 - STEPPE)

# Calculate statistic
df_EN <- df_flipped %>% mutate(expected = (0.84*ANATOLIA_NEOLITHIC) + (0.16*MESOLITHIC)) %>% mutate(ss = EUROPEAN_NEOLITHIC - expected)
df_BA <- df_flipped %>% mutate(expected = (0.52*STEPPE) + (0.48*EUROPEAN_NEOLITHIC)) %>% mutate(ss = BRONZE_AGE - expected)
df_H <- df_flipped %>% mutate(expected = (0.15*EUROPEAN_NEOLITHIC) + (0.85*BRONZE_AGE)) %>% mutate(ss = HISTORICAL - expected)

# Select correct columns and save output

df_EN <- df_EN %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_EN) <- c("CHROM", "ID", "REF", "ALT", "r")
fwrite(df_EN, paste0(outpath,"EUROPEAN_NEOLITHIC_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")

df_BA <- df_BA %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_BA) <- c("CHROM", "ID", "REF", "ALT", "r")
fwrite(df_BA, paste0(outpath,"BRONZE_AGE_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")

df_H <- df_H %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_H) <- c("CHROM", "ID", "REF", "ALT", "r")
fwrite(df_H, paste0(outpath,"HISTORICAL_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")



