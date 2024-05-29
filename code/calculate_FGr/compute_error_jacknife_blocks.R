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

# Compute LOCO gamma
nblocks <- ncol(data)
gammas <- matrix(NA, nrow = nrow(data), ncol = nblocks)
for (i in 1:nblocks) {

  print(paste0("block num ",i))
  # Drop ith  column
  loco <- data[,-i]

  # Compute mean across entries
  means <- apply(loco, 1, mean,na.rm=TRUE)
  gammas[,i] <- means

}
print("Computed LOCO gamma")

# Calculate error
Fbar <- apply(gammas, 1, mean)
sigmas2 <- rep(0, nrow(data))
for (i in 1:nrow(data)) {

  sigmas2[i] <- ((nblocks - 1)/nblocks) * sum((gammas[i,] - Fbar[i])^2)

}
print(mean(sigmas2, na.rm = TRUE))

FGr_hat <- apply(data, 1, mean)

# Get test stat
f <- mean(FGr_hat)
f2 <- f^2

# Get variance
jkVar <- mean(sigmas2, na.rm = TRUE)

# Test for significance
pval <- pchisq(f2, df =1, lower.tail = F)

# Find proportion
varFGr <- var(FGr_hat, na.rm = TRUE)
error <- jkVar / varFGr

# Make output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
colnames(dfOut) <- c("f", "f2", "jkVar", "pval", "varFGr", "error")
dfOut[1,] <- c(f,f2,jkVar,pval,varFGr,error)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







