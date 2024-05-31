## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
outfile = args[2]


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
print(head(data))

# Find total number of SNPs
snp_nums <- rep(NA,22)
snp_nums[1] <- nrow(fread(paste0(snp_prefix, "_1_allSNPs.txt")))
for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  snp_nums[i] <- paste0(snp_prefix,"_", i, "_allSNPs.txt")

}

# Set parameters
L <- sum(snp_nums)
M <- nrow(data)
print(paste0("L is ", L))
print(paste0("M is ", M))

# Compute D
FGr_hat <- apply(data_block, 1, sum) * (1/L)
D <- t(FGr_hat) %*% FGr_hat

# Expected D
expD <- (M-1) / (L)


# Compute LOCO gamma
nblocks <- ncol(data)
nblocks = 5
allFGrs <- matrix(NA, nrow = nrow(data), ncol = nblocks)
allDs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  print(paste0("block num ",i))
  # Drop ith  column
  loco <- data[,-i]

  # Compute loco FGR
  FGr_loco <- apply(loco, 1, sum,na.rm=TRUE) * (1/L)
  allFGrs[,i] <- FGr_loco

  # Compute Loco D
  D_loco <- t(FGr_loco) %*% FGr_loco
  allDs <- D_loco

}
print("Computed LOCO gamma")

# Calculate variance of D
varD <-  ((nblocks -1)/nblocks) * sum((allDs - mean(allDs))^2)


# Calculate variance of all entries of FGr
Fbar <- apply(allFGrs, 1, mean)
sigmas <- rep(0, nrow(data))
for (i in 1:nrow(data)) {

  sigmas[i] <- ((nblocks - 1)/nblocks) * sum((allFGrs[i,] - Fbar[i])^2)

}
jkVar <- mean(sigmas, na.rm = TRUE)


# Test D for significance
pval <- pnorm(abs(D -expD) ,mean =0, sd = sqrt(varD), lower.tail = FALSE)

# Find Error
varFGr <- var(FGr_hat, na.rm = TRUE)
error <- jkVar / varFGr

# Make output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
colnames(dfOut) <- c("D", "varD", "pvalD", "jkVar", "varFGr", "error")
dfOut[1,] <- c(D,varD,pval,jkVar,varFGr,error)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







