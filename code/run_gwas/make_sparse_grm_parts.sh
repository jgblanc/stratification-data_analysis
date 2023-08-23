#!/bin/sh

chrList=$1
snpList=$2
nparts=$3
outPrfx=$4

for i in {1..$nparts}
do
    echo $i
    ./bin/gcta64 --mpfile $chrList --extract $snpList --make-grm-part $nparts $i --thread-num 10 --out $outPrfx
done
