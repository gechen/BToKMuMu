eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_1.root ./
./sel data cut1 BToKMuMu_1.root Data_2012A_bfAnti_cut1_1.root 
rm -rf BToKMuMu_1.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_2.root ./
./sel data cut1 BToKMuMu_2.root Data_2012A_bfAnti_cut1_2.root 
rm -rf BToKMuMu_2.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_3.root ./
./sel data cut1 BToKMuMu_3.root Data_2012A_bfAnti_cut1_3.root 
rm -rf BToKMuMu_3.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_4.root ./
./sel data cut1 BToKMuMu_4.root Data_2012A_bfAnti_cut1_4.root 
rm -rf BToKMuMu_4.root

hadd Data_2012A_bfAnti_cut1.root Data_2012A_bfAnti_cut1_*
rm -rf Data_2012A_bfAnti_cut1_*



