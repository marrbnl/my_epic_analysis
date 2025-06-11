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
	#outfile=ab.$outfile
	outfile=`echo $outfile | sed "s/.hepmc//g"`
	outfile=${indir}/${outfile}.hepmc3.tree.root
	echo $file $outfile
	~/Util/hepmc3ascii2root/hepmc3ascii2root $file $outfile
done

