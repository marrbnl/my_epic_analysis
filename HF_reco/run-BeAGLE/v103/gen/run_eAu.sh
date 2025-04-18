#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_eAu.sh jobnumber"
        echo "Exiting..."
        exit 1
fi

#Go into scratch directory
chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

#Working directory
WDIR=${2}

#Make subdirectory and move there
INPUT=$(( 0 + $1 ))
echo $INPUT
DIR=`printf "%04d" $INPUT`
mkdir $DIR
cd $DIR

#Environmental Variables
export BEAGLESYS="/gpfs/mnt/gpfs02/eic/mdbaker/BeAGLEdev/BeAGLE"

#Soft links to necessary files
ln -s ${WDIR}/eAu.inp
ln -s ${WDIR}/S3ALL003
ln -s ${WDIR}/nuclear.bin
ln -s ${WDIR}/make_tree.C
ln -s ${WDIR}/HF_filter.C
ln -s ${WDIR}/bins.h

#Run simulation
echo "start running in directory $PWD"

echo "Running Job Number $1"
$BEAGLESYS/BeAGLE_1.03.01 < eAu.inp > eA.log
echo "Completed Simulation!!!"

echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eA.txt\")"
echo "Done!!!"
echo ""

echo "Filter ROOT File..."
root -l -b -q "HF_filter.C(\"eA.root\")"
echo "Done!!!"
echo ""

#Move output files and cleanup
echo "Cleaning Up..."
#mv -v eA.txt  $WDIR/outForPythiaMode/eAu_${INPUT}.txt
#mv -v eA.root $WDIR/outForPythiaMode/eAu_${INPUT}.root
mv -v eA.filterD0.root $WDIR/outForPythiaMode/eAu_${INPUT}.filterD0.root
mv -v eA.log  $WDIR/log/eAu_${INPUT}.log
echo "DONE!!!"
