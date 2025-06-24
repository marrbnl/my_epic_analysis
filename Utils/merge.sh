#!/bin/bash

# NOTE: for this script to work, the following needs to be run
# setenv EIC_LEVEL dev
# source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_cshrc.csh
# they are included in .cshrc

indir=$1
odir=$2
nfiles=$3
prefix=$4

if [ -z $indir ]; then
    echo "[e] no input directory"
    exit
fi

if [ -z $odir ]; then
    echo "[e] no output directory"
    exit
fi

if [ -z $nfiles ]; then
    echo "[e] please input # of files to merge"
    exit
fi

if [ -z $prefix ]; then
    echo "[e] please input prefix for merged files"
    exit
fi

filelist=files.all.txt
cp blank $filelist

ls $indir/*.root > $filelist

split -d -l $nfiles --suffix-length=4 $filelist sublist.

for list in `ls sublist.*`; do
    num=`echo $list | cut -d "." -f 2`
    echo $list $num
    sources=`cat $list`
    outfile=$prefix.$num.root
    hadd $outfile $sources
    mv $outfile $odir
    rm $list
    #cat $list | hadd $prefix
done

