Instructions for running pca_global_impactpoint.C
It has to run inside the eic-shell since it uses DD4Hep and Acts classes. You run as

eic-shell -c eic_xl -v 24.09-stable
source /opt/detector/epic-24.09.0/bin/thisepic.sh
root -l -b -q pca_global_impactpoint.C

The code uses an EICRecon ROOT file as input. 
