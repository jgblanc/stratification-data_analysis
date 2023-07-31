#!/bin/bash
pfile_path=$1
outfile=$2

plink2 \
  --pfile $pfile_path \
  --geno-counts \
  --out $outfile
