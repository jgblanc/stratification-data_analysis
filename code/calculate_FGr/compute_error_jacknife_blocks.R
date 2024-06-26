## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
outfile = args[2]
snp_prefix = args[3]

# Get list of all chromosomes
chrs <- seq(1,22)

# Add Tm for each chromosome to each other
data <- fread(paste0(Tm_prefix, "_1.txt"))
data <- data[,3:ncol(data)]

for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  filename <- paste0(Tm_prefix,"_", i, ".txt")
  tmp <- fread(filename)
  tmp <- tmp[,3:ncol(tmp)]
  data <- cbind(data, tmp)

}
data <- as.data.frame(data)
print(dim(data))

# Find total number of SNPs
snp_nums <- fread(paste0(snp_prefix, "_1_SNPcount.txt"))
for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  snp_nums <- rbind(snp_nums, fread(paste0(snp_prefix,"_", i, "_SNPcount.txt")))

}
print(dim(snp_nums))
print(snp_nums)

# Set parameters
L <- sum(snp_nums$nSNP)
M <- nrow(data)
print(paste0("L is ", L))
print(paste0("M is ", M))

# Compute D
FGr_hat <- apply(data, 1, sum) * (1/L)
D <- t(FGr_hat) %*% FGr_hat * (L^2)

# Expected D
expD <- (M-1)

# Compute LOCO
nblocks <- ncol(data)
allFGrs <- matrix(NA, nrow = nrow(data), ncol = nblocks)
allDs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  print(paste0("block num ",i))
  # Drop ith  column
  loco <- data[,-i]

  # Drop block nsnps
  block_snps <- snp_nums[i, 2]
  locoL <- as.integer(L - block_snps)
  print(locoL)

  # Compute loco FGR
  FGr_loco <- apply(loco, 1, sum,na.rm=TRUE) * (1/locoL) # Fix to get SNPs without block
  print(str(FGr_loco))
  allFGrs[,i] <- FGr_loco

  # Compute Loco D
  D_loco <- t(FGr_loco) %*% FGr_loco * (locoL^2) # Fix to get SNPs without block
  allDs[i] <- D_loco

}
print("Computed LOCO gamma")
print(allDs)
print(mean(allDs))
# Calculate variance of D
varD <-  ((nblocks -1)/nblocks) * sum((allDs - mean(allDs, na.rm=TRUE))^2)


# Calculate variance of all entries of FGr
Fbar <- apply(allFGrs, 1, mean)
sigmas <- rep(0, nrow(data))
for (i in 1:nrow(data)) {

  sigmas[i] <- ((nblocks - 1)/nblocks) * sum((allFGrs[i,] - Fbar[i])^2)

}
jkVar <- mean(sigmas, na.rm = TRUE)


# Test D for significance
pval <- pnorm(abs(D -expD) ,mean =0, sd = sqrt(varD), lower.tail = FALSE) * 2

# Find Error
varFGr <- var(FGr_hat, na.rm = TRUE)
error <- jkVar / varFGr

# Find Z
Z <- D / ((M-1) * L)


# Make output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 8))
colnames(dfOut) <- c("D","ExpD", "varD", "pvalD","Z", "jkVar", "varFGr", "error")
dfOut[1,] <- c(D,expD,varD,pval,Z, jkVar,varFGr,error)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







