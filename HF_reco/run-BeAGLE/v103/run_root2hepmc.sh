#!/bin/bash

# NOTE: for this script to work, the following needs to be run
# setenv EIC_LEVEL dev
# source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_cshrc.csh
# they are included in .cshrc

indir=$1

if [ -z $indir ]; then
    echo "[e] no input directory"
    exit
fi

files=`ls $indir/*.filterD0.root`
for file in $files; do
    echo $file
    outfile=$(basename $file)
    #outfile=ab.$outfile
    filenumber=`echo $outfile | sed "s/root/hepmc/g"`
    echo $filenumber
    echo "TreeToHepMC(\"${file}\")" | /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/env/dev/PACKAGES/eic-smear/build/eic-smear |& grep -vi exception | grep -v skipped
done

