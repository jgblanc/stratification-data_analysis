#!/bin/sh

chrList=$1
snpList=$2
nparts=$3
indx=$4
outPrfx=$5

./bin/gcta64 --mpfile $chrList --extract $snpList --make-grm-part $nparts $indx --thread-num 10 --out $outPrfx
