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
    source /opt/detector/epic-25.03.1/bin/thisepic.sh 
    npsim --compactFile $DETECTOR_PATH/$DETECTOR_CONFIG.xml --numberOfEvents 20 --inputFiles "~/gpfs02/BeAGLE/D0/eAu_10x100_Q2min1/eAu_0.hepmc3.tree.root" --outputFile "test.edm4hep.root"
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    odir=/gpfs/mnt/gpfs02/eic/rongrong/BeAGLE/D0/eAu_10x100_Q2min1
    logdir=$odir/edm4hep/log
    if [ ! -d $logdir ]; then
	mkdir -pv $logdir
    fi
    rm -rf $logdir/*

    executable=run_npsim_BeAGLE.sh
    cp -v ${executable} $odir/.
    chmod a+x $odir/${executable}
    
    echo $odir/$executable
    #Initialising Condor File
    condor_file=CondorFile_npsim_BeAGLE
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv  =  False" >> ${condor_file}
    echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:nightly"' >> ${condor_file}

    hepmcdir=$odir/ab_hepmc
    files=`ls $hepmcdir/*hepmc`
    for file in $files; do
	#for ifile in `seq 0 447`; do
	#for ifile in `seq 0 100`; do
	for ifile in `seq 0 0`; do
	    ArgName=`basename ${file} | sed "s/.hepmc//g"`
	    LogFile=${logdir}/${ArgName}.${ifile}.out
	    ErrFile=${logdir}/${ArgName}.${ifile}.err
	    SkipEvents=`expr $ifile \* 400`
	    
	    echo "" >> ${condor_file}
	    echo "Output    = ${LogFile}" >> ${condor_file}
	    echo "Error     = ${ErrFile}" >> ${condor_file}
	    echo "Arguments = ${ArgName} ${ifile} ${SkipEvents}" >> ${condor_file}
	    echo "Queue" >> ${condor_file}
	done
    done
    cp ${condor_file} $odir/.
fi

