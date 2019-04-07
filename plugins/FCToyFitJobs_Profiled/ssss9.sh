#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_9.py angular2D_Profiled_Afb 9
#./nll NLLwFix test 9
#python fitFCToys_unCons_ToyMC.py 9
python fitFCToys_unCons_NLLMC.py 9

