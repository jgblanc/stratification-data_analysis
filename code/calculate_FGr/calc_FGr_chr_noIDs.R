## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript calc_FGr_chr.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gwas_prefix = args[1]
out_prefix = args[2]
r_file = args[3]
snps_file = args[4]
outfile = args[5]


# Read in GWAS individuals
gwasID <- fread(paste0(gwas_prefix, ".psam"))
m <- nrow(gwasID)

# Read in and format r
r <- fread(r_file)
dfSnps <- fread(snps_file)
colnames(dfSnps) <- "ID"
r <- inner_join(r, dfSnps)
print(nrow(r))
r <- r %>% dplyr::select("ID", "ALT", "r")
colnames(r) <- c("ID", "A1", "BETA")
print(nrow(r))

# Compute GWAS genotype counts
outfile_count <- paste0(out_prefix, "G_count")
cmd_count <- paste("sh code/calculate_FGr/compute_GWAS_count.sh", gwas_prefix, outfile_count, snps_file, sep = " ")
print(cmd_count)
system(cmd_count)

# Calculate variance of GWAS panel genotypes from counts
count_plink <- fread(paste0(out_prefix, "G_count.gcount"))
nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
mean_gc <- counts / nOBS
length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)
length_mc_genos <- length_mc_genos * (1/(m-1))

#  Re-write .linear file with correct betas
r$BETA <- r$BETA * (1/length_mc_genos)
r[is.na(r)] <- 0
r[is.infinite(r$BETA),3] <- 0

# Save r to use as scoring weights
fwrite(r, paste0(out_prefix, ".xt_temp.glm.linear"), sep = "\t")

# Compute Gr
cmd_b <- paste("sh code/calculate_FGr/GWAS_score.sh",
               gwas_prefix,
               paste0(out_prefix, ".xt_temp.glm.linear"),
               paste0(out_prefix,".gxt_tmp"), snps_file, sep = " ")
print(cmd_b)
system(cmd_b)

# Read in FGr
dfFGr = fread(paste0(out_prefix, ".gxt_tmp.sscore"))
print(head(dfFGr))
print(nrow(dfFGr))
FGr = as.matrix(dfFGr$BETA_SUM)

# Format output
gwasID <- dfFGr %>% select("#IID")
colnames(gwasID) <- "#FID"
gwasID$IID <- gwasID$`#FID`
gwasID$FGr <- FGr
print(head(gwasID))

# Save output
fwrite(gwasID, outfile, row.names = F, col.names = T, quote = F, sep = "\t")

# Remove tmp files
cmd <- paste("rm", paste0(out_prefix, ".xt_temp.glm.linear*"),  paste0(out_prefix,".gxt_tmp*" ), paste0(out_prefix,"G_count*" ), sep = " ")
system(cmd)
