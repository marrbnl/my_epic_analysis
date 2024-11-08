set maindir = /gpfs/mnt/gpfs02/eic/hsingh1/khushi/eic/big_stat
cd $maindir
echo $maindir

set executable = `pwd`/run_npsim_reco.sh
echo $executable
#Initialising Condor File
echo "" > CondorFile
echo "Universe    = vanilla" >> CondorFile
echo "Executable  = ${executable}" >> CondorFile
echo "GetEnv  =  False" >> CondorFile
echo '+SingularityImage="/cvmfs/singularity.opensciencegrid.org/eicweb/jug_xl:nightly"' >> CondorFile

set inputdir = $maindir/pythiaDIS_18x275_minQ2=10/pythiaDIS_18x275_minQ2=10_*_1000.hepmc

foreach file(${inputdir}*)
	set ArgName = `basename ${file} | sed "s/_1000.hepmc//g"`
	set LogFile = ${maindir}/log_files/${ArgName}.out
	set ErrFile = ${maindir}/err_files/${ArgName}.err
	
	echo "" >> CondorFile
	echo "Output    = ${LogFile}" >> CondorFile
	echo "Error     = ${ErrFile}" >> CondorFile
	echo "Arguments = ${ArgName}" >> CondorFile
	echo "Queue" >> CondorFile
end

