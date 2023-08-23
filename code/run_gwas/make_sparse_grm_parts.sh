#!/bin/sh

chrList=$1
snpList=$2
nparts=$3
outPrfx=$4

#for i in $( seq 2 $nparts )
for i in {3..3}
do
    echo $i
    ./bin/gcta64 --mpfile $chrList --extract $snpList --make-grm-part $nparts $i --thread-num 10 --out $outPrfx
done
