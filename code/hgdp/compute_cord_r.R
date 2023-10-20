## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript project_Tvec_chr.R <test panel prefix> <gwas panel prefix> <test vec file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

test_prefix = args[1]
tvec_file = args[2]
out_prefix = args[3]
overlap_snps = args[4]
outfile_lat = args[5]
outfile_long = args[6]

####################
## Functions #######
####################

# Compute G %*% t(X) %*% T
compute_b <- function(path_to_test, path_to_testvec, test_type, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  pheno_file <- paste0(path_to_testvec, "_", test_type, ".txt")
  cmd_XT <- paste("sh code/calculate_Tm/compute_XT.sh", path_to_test, pheno_file, test_type, outfile_XT, overlap_snps, sep = " ")
  system(cmd_XT)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "xt_temp." , test_type ,".glm.linear"))
  count_plink <- fread(paste0(outpath, "xt_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  betas_plink_norm <- beta_plink$BETA * length_mc_genos

  #  Re-write .linear file with correct betas
  beta_plink$BETA <- betas_plink_norm
  beta_reformat <- beta_plink %>% dplyr::select(ID, A1, BETA)
  beta_reformat[is.na(beta_reformat)] <- 0

  return(beta_reformat)
}


#####################
##     Main       ###
#####################


# Gather parameters
testID <- fread(paste0(test_prefix, ".psam"))
head(testID)
colnames(testID) <- c("IID",  "Sex")
n <- nrow(testID)

# Compute b for latitude
r_lat = compute_b(path_to_gwas = gwas_prefix, path_to_testvec = tvec_file, test_type = "latitude", outpath = out_prefix)
print(head(r_lat))

# Compute b for latitude
b_long = compute_b(path_to_test = test_prefix, path_to_gwas = gwas_prefix, path_to_testvec = tvec_file, test_type = "longitude", outpath = out_prefix)
b_long = as.data.frame(b_long)
colnames(b_long) <- "longitude"

b <- cbind(b_lat, b_long)


fwrite(b, outfile, row.names = F, col.names = T, quote = F, sep = "\t")
