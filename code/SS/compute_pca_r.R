## Compute r for each chromosome

args=commandArgs(TRUE)

if(length(args)<8){stop("Rscript project_Tvec_chr.R <test panel prefix> <gwas panel prefix> <test vec file> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

test_prefix = args[1]
out_prefix = args[2]
overlap_snps = args[3]
id_file = args[4]
tvec_file_1 = args[5]
tvec_file_2 = args[6]
tvec_file_3 = args[7]
tvec_file_4 = args[8]
tvec_file_5 = args[9]
outfile_1 = args[10]
outfile_2 = args[11]
outfile_3 = args[12]
outfile_4 = args[13]
outfile_5 = args[14]

####################
## Functions #######
####################

# Compute G %*% t(X) %*% T
compute_b <- function(path_to_test, testvec_file, test_type, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  pheno_file <- testvec_file
  cmd_XT <- paste("sh code/calculate_FGr/compute_XT_IDs.sh", path_to_test, pheno_file, test_type, outfile_XT, overlap_snps, id_file, sep = " ")
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
  beta_reformat <- beta_plink %>% dplyr::select("#CHROM","ID", "REF", "ALT",  "BETA")
  beta_reformat[is.na(beta_reformat)] <- 0

  return(beta_reformat)
}


#####################
##     Main       ###
#####################


# Gather parameters
testID <- fread(id_file)
head(testID)
n <- nrow(testID)

# Compute b for pcs
tmp <- fread(tvec_file_1)
test_type <- colnames(tmp)[3]
r_1 = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_1, test_type = test_type, outpath = out_prefix)
r_1$BETA  <- r_1$BETA * (1/n)
colnames(r_1)[5] <- "r"
print(head(r_1))
fwrite(r_1, outfile_1, row.names = F, col.names = T, quote = F, sep = "\t")

tmp <- fread(tvec_file_2)
test_type <- colnames(tmp)[3]
r_2 = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_2, test_type = test_type, outpath = out_prefix)
r_2$BETA  <- r_2$BETA * (1/n)
colnames(r_2)[5] <- "r"
print(head(r_2))
fwrite(r_2, outfile_2, row.names = F, col.names = T, quote = F, sep = "\t")

tmp <- fread(tvec_file_3)
test_type <- colnames(tmp)[3]
r_3 = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_3, test_type = test_type, outpath = out_prefix)
r_3$BETA  <- r_3$BETA * (1/n)
colnames(r_3)[5] <- "r"
print(head(r_3))
fwrite(r_3, outfile_3, row.names = F, col.names = T, quote = F, sep = "\t")

tmp <- fread(tvec_file_4)
test_type <- colnames(tmp)[3]
r_4 = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_4, test_type = test_type, outpath = out_prefix)
r_4$BETA  <- r_4$BETA * (1/n)
colnames(r_4)[5] <- "r"
print(head(r_4))
fwrite(r_4, outfile_4, row.names = F, col.names = T, quote = F, sep = "\t")


tmp <- fread(tvec_file_5)
test_type <- colnames(tmp)[3]
r_5 = compute_b(path_to_test = test_prefix, testvec_file = tvec_file_5, test_type = test_type, outpath = out_prefix)
r_5$BETA  <- r_5$BETA * (1/n)
colnames(r_5)[5] <- "r"
print(head(r_5))
fwrite(r_5, outfile_5, row.names = F, col.names = T, quote = F, sep = "\t")
