#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python scan_update.py angular2D_Profiled_Afb 0
#./nll2 NLLwFix test 0
python fitFCToys_unCons_NLLMC.py 0

