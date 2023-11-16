#!/bin/bash
pfile_path=$1
outfile=$2
overlap_snps=$3
IDs=$4

plink2 \
    --pfile $pfile_path \
    --keep $IDs \
    --extract $overlap_snps \
    --geno-counts \
    --out $outfile
