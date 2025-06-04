#!/bin/bash

indir=$1

if [ -z $indir ]; then
    echo "[i] No directory is input. Exit!"
    exit
fi

#for file in `ls $indir/*.root`; do
for file in `ls $indir/*.hepmc`; do
    file_size_kb=`du -k "$file" | cut -f1`
    if [ $file_size_kb -lt 5000 ]; then
	rm -v $file
	#echo $file $file_size_kb
    fi
done
