# Script to combine all chromosome from vilma output and readjust SNPs IDs to only include CHR:POS

args=commandArgs(TRUE)

if(length(args)<4){stop("Provide <list of chromosome numbers> <first prefix> <second prefix> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

chrlist = args[1]
p1 = args[2]
p2 = args[3]
outfile = args[4]

print(chrlist)
print(p1)
print(p2)


