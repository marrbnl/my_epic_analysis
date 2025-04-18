#!/bin/bash

condorDir=/gpfs/mnt/gpfs02/eic/rongrong/BeAGLE/D0/eAu_10x100_Q2min1
cd ${condorDir}

CONFIG=${1}
IFILE=${2}
SKIP_NEVENTS=${3}

HEPMC_FILE=${condorDir}/ab_hepmc/${CONFIG}.hepmc
SIM_FILE=${condorDir}/edm4hep/sim_${CONFIG}.${IFILE}.edm4hep.root
LOG=${condorDir}/edm4hep/${CONFIG}.${IFILE}.log

pwd

#${condorDir}/eic-shell << EOF
ls /opt/detector

source /opt/detector/epic-25.03.1/bin/thisepic.sh epic

env

npsim --compactFile ${DETECTOR_PATH}/${DETECTOR_CONFIG}.xml  \
      --inputFiles ${HEPMC_FILE}  \
      --outputFile ${SIM_FILE}  \
      --runType batch  \
      --numberOfEvents 400 \
      --skipNEvents ${SKIP_NEVENTS}

#eicrecon ${SIM_FILE} -Ppodio:output_file=${RECO_FILE}

exit

