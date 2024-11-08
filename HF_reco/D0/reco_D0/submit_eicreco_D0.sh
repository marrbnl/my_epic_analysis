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
    #source /opt/detector/epic-24.05.0/setup.sh
    source /opt/detector/epic-24.09.0/bin/thisepic.sh epic
    
    #eicrecon ~/gpfs02/vertexing/pythiaDIS_18x275_minQ210/edm4hep/sim_pythiaDIS_18x275_minQ2=10_100.edm4hep.root -Ppodio:output_file=eicrecon_output.test.root
    #eicrecon test.edm4hep.root -Pacts:LogLevel=debug -Ppodio:output_file=eicrecon_output.test.root
    #eicrecon test.edm4hep.root -Ptracking:CentralTrackVertices:LogLevel=debug -Ppodio:output_file=eicrecon_output.test.root
    eicrecon test.edm4hep.root -Ppodio:output_file=eicrecon_output.test.root
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    config=Geo202409_Real_default
    
    echo "[i] Running config = $config"
    
    odir=/gpfs/mnt/gpfs02/eic/rongrong/D0
    outdir=$odir/$config
    edmdir=$odir/edm4hep_ab

    logdir=$outdir/log
    if [ ! -d $logdir ]; then
	mkdir -pv $logdir
    fi
    rm -rf $logdir/*

    cp -v run_eicreco_D0.sh $odir/.
    chmod a+x $odir/run_eicreco_D0.sh
    executable=$odir/run_eicreco_D0.sh
    echo $executable
    #Initialising Condor File
    condor_file=CondorFile_eicrecon_D0
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${executable}" >> ${condor_file}
    echo "GetEnv  =  False" >> ${condor_file}
    echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/jug_xl:nightly"' >> ${condor_file}

    files=`ls $edmdir/*.root`
    for file in $files; do
	ArgName=`basename ${file} | sed "s/.edm4hep.root//g" | sed "s/sim_//g"`
	LogFile=${logdir}/${ArgName}.out
	ErrFile=${logdir}/${ArgName}.err
	
	echo "" >> CondorFile_eicrecon
	echo "Output    = ${LogFile}" >> ${condor_file}
	echo "Error     = ${ErrFile}" >> ${condor_file}
	echo "Arguments = ${ArgName} ${config}" >> ${condor_file}
	echo "Queue" >> ${condor_file}
    done
    cp ${condor_file} $odir/.
fi

