#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_10.py angular2D_Profiled_Afb 10
#./nll NLLwFix test 10
#python fitFCToys_unCons_ToyMC.py 10
python fitFCToys_unCons_NLLMC_0.py 10
