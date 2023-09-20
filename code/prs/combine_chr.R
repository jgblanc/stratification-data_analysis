# Script to combine all chromosome from vilma output and readjust SNPs IDs to only include CHR:POS

args=commandArgs(TRUE)

if(length(args)<5){stop("Provide <number of total chrosomesome> <first prefix> <second prefix> <output> <list of chromosome numbers> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

chr_num = as.numeric(args[1])
p1 = args[2]
p2 = args[3]
outfile = args[4]


print(chr_num)
print(p1)
print(p2)

# Read in first data frame
fn <- paste0(p1, "_", args[5], p2)
print(fn)
df <- fread(fn)

# Fix IDs
df <- df %>% separate(ID, c("ID", "tmp1", "tmp2"), "_") %>%
  select(ID,A1, A2, posterior_ukbb, posterior_variance_ukbb)


for (i in 2:chr_num) {

  # Read in data frame
  fn <- paste0(p1, "_", args[4+i], p2)
  print(fn)
  tmp <- fread(fn)

  # Fix IDs
  tmp <- tmp %>% separate(ID, c("ID", "tmp1", "tmp2"), "_") %>%
    select(ID,A1, A2, posterior_ukbb, posterior_variance_ukbb)

  # Add to existing df
  df <- rbind(df, tmp)

}

# Save output
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")

