#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_2.py angular2D_Profiled_Afb 2
#./nll NLLwFix test 2
#python fitFCToys_unCons_ToyMC.py 2
python fitFCToys_unCons_NLLMC.py 2
