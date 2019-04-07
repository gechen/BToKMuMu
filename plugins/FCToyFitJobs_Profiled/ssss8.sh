#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_8.py angular2D_Profiled_Afb 8
#./nll NLLwFix test 8
#python fitFCToys_unCons_ToyMC.py 8
python fitFCToys_unCons_NLLMC.py 8

