## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_FGr_chr.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

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
colnames(gwasID) <- c("FID", "IID",  "Sex")

# Read in and format r
r <- fread(r_file)
r <- r %>% dplyr::select("ID", "ALT", "r")
colnames(r) <- c("ID", "A1", "BETA")
print(head(r))

# Save r to use as scoring weights
fwrite(r, paste0(out_prefix, "xt_temp.glm.linear"), sep = "\t")

# Compute Gr
cmd_b <- paste("sh code/calculate_FGr/GWAS_score.sh",
               gwas_prefix,
               paste0(out_prefix, "xt_temp.glm.linear"),
               paste0(out_prefix,"gxt_tmp" ), snps_file, sep = " ")
print(cmd_b)
system(cmd_b)

# Read in FGr
FGr = fread(paste0(out_prefix, "gxt_tmp.sscore"))
FGr = as.matrix(FGr$BETA_SUM)

# Format output
gwasID$FGr <- FGr

# Save output
fwrite(gwasID, outfile, row.names = F, col.names = T, quote = F, sep = "\t")
