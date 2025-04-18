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
    echo "[i] Running local tests"
    export BEAGLESYS="/gpfs02/eic/mdbaker/BeAGLEdev/BeAGLE"
    $BEAGLESYS/BeAGLE_1.03.01 < eAu.inp > eA.log
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    odir=/gpfs/mnt/gpfs02/eic/rongrong/BeAGLE/D0/eAu_10x100_Q2min1
    logdir=$odir/outForPythiaMode/log
    if [ ! -d $logdir ]; then
	mkdir -pv $logdir
    fi
    rm -rf $logdir/*

    executable=run_eAu.sh
    cp -v ${executable} $odir/.
    cp -v eAu.inp $odir/.
    cp -v S3ALL003 $odir/.
    cp -v nuclear.bin $odir/.
    cp -v make_tree.C $odir/.
    cp -v HF_filter.C $odir/.
    cp -v bins.h $odir/.
    chmod a+x $odir/${executable}
    
    echo $odir/$executable
    #Initialising Condor File
    condor_file=CondorFile_gen
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv  =  True" >> ${condor_file}
    #echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/jug_xl:nightly"' >> ${condor_file}

    for ifile in `seq 0 999`; do
	LogFile=${logdir}/eAujob.${ifile}.out
	ErrFile=${logdir}/eAujob.${ifile}.err
	    
	echo "" >> CondorFile_gen
	echo "Output    = ${LogFile}" >> ${condor_file}
	echo "Error     = ${ErrFile}" >> ${condor_file}
	echo "Arguments = ${ifile} ${odir}" >> ${condor_file}
	echo "Queue" >> ${condor_file}
    done
    cp ${condor_file} $odir/.
fi

