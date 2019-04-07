#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

#python fitFCToys_unCons_ToyMC.py 7
python scan_7.py angular2D_Profiled_Afb 7
