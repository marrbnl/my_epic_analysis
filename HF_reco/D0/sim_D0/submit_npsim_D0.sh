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
    source /opt/detector/epic-24.09.0/bin/thisepic.sh 
    npsim --compactFile $DETECTOR_PATH/$DETECTOR_CONFIG.xml --numberOfEvents 20 --inputFiles "/eic/u/rongrong/gpfs02/D0/D0_ABCONV/pythia8.306-1.0/18x275/hiAcc/pythia8.306-1.0_D0_18x275_hiAcc_run1.hepmc3.tree.root" --outputFile "test.edm4hep.root"
    exit
fi

###############################################################
### Batch prodcution for real data using file list          ###
###############################################################

if [ $submit -eq 1 ]; then
    odir=/gpfs/mnt/gpfs02/eic/rongrong/D0
    outdir=$odir/edm4hep_ab
    hepmcdir=$odir/D0_ABCONV/pythia8.306-1.0/18x275/hiAcc

    logdir=$outdir/log
    if [ ! -d $logdir ]; then
	mkdir -pv $logdir
    fi
    rm -rf $logdir/*

    executable=run_npsim_reco_D0.sh
    cp -v ${executable} $odir/.
    chmod a+x $odir/${executable}
    
    echo $odir/$executable
    #Initialising Condor File
    condor_file=CondorFile_npsim_D0
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv  =  False" >> ${condor_file}
    echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/jug_xl:nightly"' >> ${condor_file}

    files=`ls $hepmcdir`
    for file in $files; do
	#for ifile in `seq 0 447`; do
	#for ifile in `seq 0 100`; do
	for ifile in `seq 101 447`; do
	    ArgName=`basename ${file} | sed "s/.hepmc3.tree.root//g"`
	    LogFile=${logdir}/${ArgName}.${ifile}.out
	    ErrFile=${logdir}/${ArgName}.${ifile}.err
	    SkipEvents=`expr $ifile \* 400`
	    
	    echo "" >> CondorFile_npsim
	    echo "Output    = ${LogFile}" >> ${condor_file}
	    echo "Error     = ${ErrFile}" >> ${condor_file}
	    echo "Arguments = ${ArgName} ${ifile} ${SkipEvents}" >> ${condor_file}
	    echo "Queue" >> ${condor_file}
	done
    done
    cp ${condor_file} $odir/.
fi

