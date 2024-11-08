#!/bin/bash

submit=$1

if [ -z $submit ]; then
    echo "[i] Set it to submit = 0 for local test"
    submit=0
fi

##################################################
### Local test ###
##################################################
if [ $submit -eq 0 ]; then
    echo "[i] Running locally"
    root -l -q 'ana_d0.C("file.list","test.root")'
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    #configs=(Geo202406_Truth_default Geo202406_Truth_d0z0 Geo202406_Truth_d0z0_usePt Geo202406_Real_default Geo202406_Real_d0z0 Geo202406_Real_d0z0_usePt)
    configs=(Geo202409_Real_default)

    for config in "${configs[@]}"; do
	echo "[i] Running config = $config"

	filelist=input_files/file.${config}.list
	rm $filelist
	find /eic/u/rongrong/gpfs02/D0/${config}/reco*.root > $filelist
	output=DIS.D0.${config}.root
	root -l -q ana_d0.C\(\"${filelist}\",\"${output}\"\)
    done
fi

