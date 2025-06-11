#!/bin/bash

# this script needs to be run within eic-shell

indir=$1
odir=$2

if [ -z $indir ]; then
    echo "[e] no input directory"
    exit
fi

if [ -z $odir ]; then
    echo "[e] no output directory"
    exit
fi

files=`ls $indir/*.hepmc`
for file in $files; do
	echo $file
	outfile=$(basename $file)
	outfile=ab.$outfile
	outfile=`echo $outfile | sed "s/.hepmc//g"`
	echo $outfile
	#abconv --plot-off -p ip6_eau_110x10 $file -o $outfile
	abconv --plot-off -p ip6_eau_41x5 $file -o $outfile
	rm $outfile.hist.root
	mv $outfile.hepmc $odir
done

