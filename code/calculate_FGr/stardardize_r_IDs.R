## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript standardized_r_IDs.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

gwas_prefix = args[1]
out_prefix = args[2]
r_prefix = args[3]
snp_prefix = args[4]
outfile = args[5]
gwas_IDs = args[6]
ldFile = args[7]


# Read in GWAS individuals
dfGWAS_IDs <- fread(gwas_IDs)[,1:2]
m <- nrow(dfGWAS_IDs)
print(paste0("m is ", m))

# Read in and format r for all chromosome
r_file <- paste0(r_prefix, "1.rvec")
r <- fread(r_file)
for (i in 2:22) {
  r_file <- paste0(r_prefix, i, ".rvec")
  tmp <- fread(r_file)
  r <- rbind(r, tmp)
}

# Read in and format all SNPs
snp_file <- paste0(snp_prefix, "1.txt")
dfSnps <- fread(snps_file)
for (i in 2:22) {
  snp_file <- paste0(snp_prefix, i, ".txt")
  tmp <- fread(snp_file)
  dfSnps <- rbind(dfSnps, tmp)
}

# Combine R and SNPs file
colnames(dfSnps) <- "ID"
r <- inner_join(r, dfSnps)
r <- r %>% dplyr::select("ID", "ALT", "r")
colnames(r) <- c("ID", "A1", "BETA")

# Separate ID into CHR and BP
r_blocks <- r %>% separate("ID", into = c("CHR", "BP"), sep = ":", remove = FALSE)
print(paste0("r has", nrow(r), " rows"))

# Read in LD block file
ld <- fread(ldFile)

# Assign SNPs to blocks
assign_SNP_to_block <- function(CHR, BP, block = ld) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  print(block_chr)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Add block info - takes a while
r_blocks <- r_blocks %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[2]), as.numeric(params[3])))) %>%
  drop_na()
print(paste0("Now r blocks has", nrow(r_blocks), " rows"))


# Compute GWAS genotype counts of all SNPs by chromosome
outlist <- list()
for (i in 1:22) {

  # Save snps on chromosome
  selected_snps <- r_blocks %>% filter(chr == i) %>% select("ID")
  r_chr <- r_blocks %>% filter(chr == i)
  snps_file =  paste0(out_prefix, "_chr", i,"_allSNPs.txt")
  fwrite(selected_snps, snps_file, row.names = F, col.names = T, quote = F, sep = "\t")

  # Get genotype counts
  outfile_count <- paste0(out_prefix,"_chr", i, "G_count")
  plinkfile <- paste0(gwas_prefix,i, "_v3")
  cmd_count <- paste("sh code/calculate_FGr/compute_GWAS_count_ID.sh", plinkfile, outfile_count, snps_file, gwas_IDs, sep = " ")
  print(cmd_count)
  system(cmd_count)

  # Calculate variance of GWAS panel genotypes from counts
  count_plink <- fread(paste0(out_prefix, "G_count.gcount"))
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)
  length_mc_genos <- length_mc_genos * (1/(m-1))

  # Store variance
  r_chr$geno_var <- length_mc_genos
  outlist[[i]] <- r_chr
}

# Scale all R values
## Unpack list
dfOut <- do.call("rbind", outlist)
print(dfOut)
# Divide r by SD and scale
dfOut$BETA <- dfOut$BETA * (1/sqrt(dfOut$geno_var))
dfOut$BETA <- dfOut$BETA - mean(dfOut$BETA)
dfOut$BETA <- dfOut$BETA / sqrt(sum(dfOut$BETA^2))
dfOut$BETA <- dfOut$BETA * (1/sqrt(dfOut$geno_var))
dfOut <- dfOut %>% select("ID", "A1", "BETA", "block")
colnames(dfOut) <- c("ID", "ALT", "r", "block")
print(dfOut)

# Save output
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")

