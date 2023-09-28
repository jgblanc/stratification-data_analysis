# Script to run selection test with jacknife error

args=commandArgs(TRUE)

if(length(args)<2){stop("Provide <formatted betas> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

infile=args[1]
outfile = args[2]


# Read in
df <- fread(infile)
colnames(df) <- c("#CHROM", "ID", "POS", "REF", "ALT", "r", "posterior_ukbb","block")


# Function to calculate \hat{q}
calc_q <- function(df) {

  B <- df$posterior_ukbb
  r <- df$r

  # compute q
  q <- t(B) %*%  r

  return(q)
}

# Set up output file
num_blocks <- length(unique(df$block))
jacknives <- rep(0, num_blocks)

# Calculate \hat{q} with LOCO
for (i in 1:num_blocks) {

  if (i %in% seq(0,2000, 10)) {
    print(i)
  }

  # Get rid of block
  block_num <- unique(df$block)[i]
  df_LOCO <- df %>% filter(block != block_num)

  # Compute q
  jacknives[i] <- calc_q(df_LOCO)

}

# Compute mean
qBar <- mean(jacknives)

# Compute sigma squared
sigma2 <- ((num_blocks -  1)/num_blocks) * sum((jacknives - qBar)^2 )

# Compute full q
q <- calc_q(df)

# Compute p-value
pval <- pnorm(abs(q), mean = 0, sd = sqrt(sigma2),lower.tail = FALSE) * 2

# Set up output table
out <- as.data.frame(matrix(nrow = 1, ncol =2))
colnames(out) <- c("q", "pval")
out[1,] <- c(q, pval)

# Save output
fwrite(out, outfile,col.names=T,row.names=F,quote=F,sep="\t")
