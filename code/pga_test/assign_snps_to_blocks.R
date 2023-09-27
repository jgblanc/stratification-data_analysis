# Script to assign snps to LD blocks and harmonize SNP directions
# .rvec is altFreq - altFreq
# A1 = ALT

args=commandArgs(TRUE)

if(length(args)<4){stop("Provide <ld block> <vilma format> <rvec> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ldFile =args[1]
ssFile = args[2]
rFile = args[3]
outfile = args[4]


# Read in LD block file
ld <- fread(ldFile)

# Read in summary stats
ss <- fread(ssFile)
dfss <- ss %>% separate(ID, c("ID","REF", "ALT"), "_")

# Read in r file
rvec <- fread(rFile)

# Combine r and ss
df <- inner_join(rvec, dfss, by = c("ID", "REF", "ALT"))
df <- df %>% separate(ID, c("chr", "POS"), remove = FALSE)

# Save header
header <- c("#CHROM", "ID", "POS", "REF", "ALT", "r", "posterior_ukbb","block")
fwrite(header, outfile_header, row.names = F, col.names = T, quote = F, sep = "\t")

# Assign SNPs to blocks

assign_SNP_to_block <- function(CHR, BP, block = ld) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Add block info - takes a while
df_blocks <- df %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[3]), as.numeric(params[4]))))

# Format output
out <- df_blocks %>% select("#CHROM", "ID", "POS", "REF", "ALT", "r", "posterior_ukbb","block")


# Save output
fwrite(out, outfile_betas, row.names = F, col.names = T, quote = F, sep = "\t")





