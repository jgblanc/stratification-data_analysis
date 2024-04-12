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
colnames(df) <- c("CHR", "ID", "POS", "REF", "ALT", "r", "BETA","block", "P")


# Function to calculate \hat{q}
calc_q <- function(df) {

  B <- df$BETA
  r <- df$r

  # compute q
  q <- t(B) %*%  r

  return(q)
}

main <- function(df) {

  # Compute full q
  q <- calc_q(df)

  flips <- rep(0, 1000)
  for (i in 1:1000) {

    dftmp <- df
    dftmp$Beta <- dftmp$Beta * sample(c(-1,1), nrow(dftmp), replace = T,prob = c(0.5, 0.5))

    flips[i] <- calc_q(dftmp)
  }

  # Compute p-value
  pval <- pnorm(abs(q), mean = 0, sd = sd(flips),lower.tail = FALSE) * 2

  # Divide q by var
  q <- calc_q(df) / sd(flips)

  x <- c(q, pval)

  return(x)
}


# Set up output table
out <- as.data.frame(matrix(nrow = 1, ncol =4))
colnames(out) <- c("q", "pval", "q_NoSign", "pval_NoSign")
out[1,1:2] <- main(df)
df$BETA <- sign(df$BETA)
out[1,3:4] <- main(df)


# Save output
fwrite(out, outfile,col.names=T,row.names=F,quote=F,sep="\t")
