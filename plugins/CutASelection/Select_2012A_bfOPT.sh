eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_1.root ./
./sel data cutJpsi BToKMuMu_1.root Data_2012A_bfAnti_cutJpsi_1.root 
./sel data cutPsip BToKMuMu_1.root Data_2012A_bfAnti_cutPsip_1.root 
rm -rf BToKMuMu_1.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_2.root ./
./sel data cutJpsi BToKMuMu_2.root Data_2012A_bfAnti_cutJpsi_2.root 
./sel data cutPsip BToKMuMu_2.root Data_2012A_bfAnti_cutPsip_2.root 
rm -rf BToKMuMu_2.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_3.root ./
./sel data cutJpsi BToKMuMu_3.root Data_2012A_bfAnti_cutJpsi_3.root 
./sel data cutPsip BToKMuMu_3.root Data_2012A_bfAnti_cutPsip_3.root 
rm -rf BToKMuMu_3.root

eos cp /eos/cms/store/user/gechen/crab3_run/DoubleMuParked/BuToKMuMu_Data_2012A_8TeV_v8/BToKMuMu_4.root ./
./sel data cutJpsi BToKMuMu_4.root Data_2012A_bfAnti_cutJpsi_4.root 
./sel data cutPsip BToKMuMu_4.root Data_2012A_bfAnti_cutPsip_4.root 
rm -rf BToKMuMu_4.root

hadd Data_2012A_bfAnti_cutJpsi.root Data_2012A_bfAnti_cutJpsi_*
hadd Data_2012A_bfAnti_cutPsip.root Data_2012A_bfAnti_cutPsip_*
rm -rf Data_2012A_bfAnti_cutJpsi_*
rm -rf Data_2012A_bfAnti_cutPsip_*



