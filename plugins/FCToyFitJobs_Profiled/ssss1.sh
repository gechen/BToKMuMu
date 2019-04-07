#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_1.py angular2D_Profiled_Afb 1
#./nll2 NLLwFix test 1
#python fitFCToys_unCons_ToyMC.py 1
python fitFCToys_unCons_NLLMC.py 1

