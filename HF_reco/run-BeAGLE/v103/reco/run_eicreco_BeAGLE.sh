#!/bin/bash

condorDir=/gpfs/mnt/gpfs02/eic/rongrong/BeAGLE/D0/eAu_10x100_Q2min1
cd ${condorDir}

#condorDir=/gpfs/mnt/gpfs02/eic/rongrong/data/pythiaDIS_18x275_minQ210/
#cd ${condorDir}

#cd /gpfs/mnt/gpfs02/eic/rongrong/data/pythiaDIS_18x275_minQ210/edm4hep
#cd /direct/eic+u/rongrong
#EICSHELL=./eic-shell

FILENAME=${1}
CONFIG=${2}

SIM_FILE=${condorDir}/edm4hep/sim_${FILENAME}.edm4hep.root

RECO_FILE=${condorDir}/${CONFIG}/reco_${FILENAME}.root

echo $SIM_FILE
echo $RECO_FILE


#LOG=${condorDir}/edm4hep/${CONFIG}.log
#DETECTOR_PATH=/opt/detector/epic-nightly/share/epic
#DETECTOR_CONFIG=epic
# ls /opt/detector
pwd 
#source /opt/detector/setup.sh

#${HOME}/bin/eic-shell

#source /opt/detector/epic-24.05.0/setup.sh epic
source /opt/detector/epic-25.03.1/bin/thisepic.sh epic

#source /gpfs/mnt/gpfs02/eic/rongrong/D0/EICrecon/bin/eicrecon-this.sh

#env


eicrecon ${SIM_FILE} -Ppodio:output_file=${RECO_FILE}

exit

