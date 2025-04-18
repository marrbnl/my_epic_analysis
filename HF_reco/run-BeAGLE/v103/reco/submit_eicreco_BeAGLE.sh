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
    echo "[i] Running eicreco locally"
    source /opt/detector/epic-25.03.1/bin/thisepic.sh epic
    eicrecon test.edm4hep.root -Ppodio:output_file=eicrecon_output.test.root
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    config=Geo25.03.1
    
    echo "[i] Running config = $config"
    
    odir=/gpfs/mnt/gpfs02/eic/rongrong/BeAGLE/D0/eAu_10x100_Q2min1
    outdir=$odir/$config
    edmdir=$odir/edm4hep

    logdir=$outdir/log
    if [ ! -d $logdir ]; then
	mkdir -pv $logdir
    fi
    rm -rf $logdir/*

    executable=run_eicreco_BeAGLE.sh
    cp -v run_eicreco_BeAGLE.sh $odir/.
    chmod a+x $odir/run_eicreco_BeAGLE.sh

    echo $executable
    #Initialising Condor File
    condor_file=CondorFile_eicrecon_BeAGLE
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv  =  False" >> ${condor_file}
    echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:nightly"' >> ${condor_file}

    files=`ls $edmdir/sim_*.edm4hep.root`
    for file in $files; do
	ArgName=`basename ${file} | sed "s/.edm4hep.root//g" | sed "s/sim_//g"`
	LogFile=${logdir}/${ArgName}.out
	ErrFile=${logdir}/${ArgName}.err
	
	echo "" >> ${condor_file}
	echo "Output    = ${LogFile}" >> ${condor_file}
	echo "Error     = ${ErrFile}" >> ${condor_file}
	echo "Arguments = ${ArgName} ${config}" >> ${condor_file}
	echo "Queue" >> ${condor_file}
    done
    cp ${condor_file} $odir/.
fi

