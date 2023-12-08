## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
outfile = args[2]


# Get list of all but the focal chromosome
chrs <- seq(1,22)

# Add Tm for each chromosome to each other
data <- fread(paste0(Tm_prefix, "_1.txt"))[,4]
colnames(data) <- "1"

for (i in 2:22) {

  # Read in new chromosome
  filename <- paste0(Tm_prefix,"_", i, ".txt")
  data$tmp <- fread(filename)[,4]
  colnames(data)[i] <- as.character(i)

}
data <- as.data.frame(data)

# Compute LOCO gamma
gammas <- matrix(NA, nrow = nrow(data), ncol = 22)
for (i in 1:22) {

  # Drop ith  column
  loco <- data[,-i]

  # Compute mean across entries
  means <- apply(loco, 1, mean)
  gammas[,i] <- means

}

# Calculate error
Fbar <- apply(gammas, 1, mean)
sigmas2 <- rep(0, nrow(data))
for (i in 1:nrow(data)) {

  sigmas2[i] <- (21/22) * sum((gammas[i,] - Fbar[i])^2)

}
FGr_hat <- apply(data, 1, sum)

# Find proportion
error <- mean(sigmas2) / var(FGr_hat)

# Make output table
dfOut <- as.data.frame(error)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")







