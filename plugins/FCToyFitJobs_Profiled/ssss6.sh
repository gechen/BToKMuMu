#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_6.py angular2D_Profiled_Afb 6
#./nll2 NLLwFix test 6
#python fitFCToys_unCons_ToyMC.py 6
python fitFCToys_unCons_NLLMC.py 6

