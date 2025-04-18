## Instructions for producing D0 eAu events with BeAGLE and ePIC software

### Step 1: generate D0 eAu events with BeAGLE
- Relevant codes are in "gen" directory
- setup.sh: set needed environment variables
- submit_eAu.sh: master script for the entire routine
- eAu.inp: specifiy collision species, collision energy, y, Q2 ranges, etc
- S3ALL003: configure PYTHIA process
- Need to run with BeAGLE v103:
```
export BEAGLESYS="/gpfs/mnt/gpfs02/eic/mdbaker/BeAGLEdev/BeAGLE"
$BEAGLESYS/BeAGLE_1.03.01 < eAu.inp > eA.log
```
- make_tree.C: convert output txt file into a tree
- HF_filter.C: filter out D0 -> pi+k events, and save them to a new root file

### Step 2: prepare files for ePIC simulation
- Set up the environment
```
setenv EIC_LEVEL dev
source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_cshrc.csh
```
- run_root2hepmc.sh: convert BeAGLE root files into HEPMC format

### Step 3: simulate events with npsim
- Relevant codes are in "sim" directory
- submit_npsim_BeAGLE.sh: master script for the routine

### Step 4: reconstruct events with EICrecon
- Relevant codes are in "reco" directory
- submit_eicreco_BeAGLE.sh: master script for the routine
