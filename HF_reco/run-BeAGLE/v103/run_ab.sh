#!/bin/bash

# this script needs to be run within eic-shell

indir=$1

if [ -z $indir ]; then
    echo "[e] no input directory"
    exit
fi

files=`ls $indir/*.hepmc`
for file in $files; do
	echo $file
	outfile=$(basename $file)
	outfile=ab.$outfile
	outfile=`echo $outfile | sed "s/.hepmc//g"`
	echo $outfile
	abconv --plot-off -p ip6_eau_110x10 $file -o $outfile
done

