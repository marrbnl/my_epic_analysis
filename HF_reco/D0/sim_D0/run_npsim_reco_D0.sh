#!/bin/bash

condorDir=/gpfs/mnt/gpfs02/eic/rongrong/D0
cd ${condorDir}

CONFIG=${1}
IFILE=${2}
SKIP_NEVENTS=${3}

HEPMC_FILE=${condorDir}/D0_ABCONV/pythia8.306-1.0/18x275/hiAcc/${CONFIG}.hepmc3.tree.root
SIM_FILE=${condorDir}/edm4hep_ab/sim_${CONFIG}.${IFILE}.edm4hep.root
LOG=${condorDir}/edm4hep_ab/${CONFIG}.${IFILE}.log

pwd 

source /opt/detector/epic-24.09.0/bin/thisepic.sh epic

env

npsim --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml  \
      --inputFiles ${HEPMC_FILE}  \
      --outputFile ${SIM_FILE}  \
      --runType batch  \
      --numberOfEvents 400 \
      --skipNEvents ${SKIP_NEVENTS}
      
#eicrecon ${SIM_FILE} -Ppodio:output_file=${RECO_FILE}

exit

