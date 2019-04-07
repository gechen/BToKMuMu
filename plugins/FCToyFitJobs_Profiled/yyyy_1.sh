#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

./the3 NLLwFix test 7
#./fit createNLLToys test 7
#python fitFCToys_unCons_NLLMC.py 7
