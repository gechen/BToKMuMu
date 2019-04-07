#!/bin/sh
date
hostname

.  ~/.bashrc

source ~/setrootenv6.sh

cd /afs/cern.ch/work/g/gechen/work/BToKMuMu/Ana/FCToyFitJobs_Profiled

./the1 NLLwFix test 7
#./nll createNLLToys test 7
#python fitFCToys_unCons_NLLMC_0.py 7
