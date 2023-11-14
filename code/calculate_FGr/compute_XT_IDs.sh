#!/bin/bash
pfile_path=$1
pheno_path=$2
test_type=$3
outfile=$4
overlap_snps=$5
ids=$6

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --glm omit-ref allow-no-covars \
  --pheno iid-only $pheno_path \
  --pheno-name $test_type \
  --geno-counts \
  --keep $ids \
  --out $outfile
