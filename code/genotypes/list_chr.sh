#!/bin/bash
prefix=$1
outfile=$2

touch $outfile
for i in {1..22}
do
    echo "${prefix}${i}_v3" >> $outfile
done
