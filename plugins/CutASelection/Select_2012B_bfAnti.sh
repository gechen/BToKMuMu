
eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_1.root ./
./sel data cut1 BToKMuMu_Data_2012B_1.root Data_2012B_bfAnti_cut1_1.root 
rm -rf BToKMuMu_Data_2012B_1.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_2.root ./
./sel data cut1 BToKMuMu_Data_2012B_2.root Data_2012B_bfAnti_cut1_2.root 
rm -rf BToKMuMu_Data_2012B_2.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_3.root ./
./sel data cut1 BToKMuMu_Data_2012B_3.root Data_2012B_bfAnti_cut1_3.root 
rm -rf BToKMuMu_Data_2012B_3.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_4.root ./
./sel data cut1 BToKMuMu_Data_2012B_4.root Data_2012B_bfAnti_cut1_4.root 
rm -rf BToKMuMu_Data_2012B_4.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_5.root ./
./sel data cut1 BToKMuMu_Data_2012B_5.root Data_2012B_bfAnti_cut1_5.root 
rm -rf BToKMuMu_Data_2012B_5.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_6.root ./
./sel data cut1 BToKMuMu_Data_2012B_6.root Data_2012B_bfAnti_cut1_6.root 
rm -rf BToKMuMu_Data_2012B_6.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_7.root ./
./sel data cut1 BToKMuMu_Data_2012B_7.root Data_2012B_bfAnti_cut1_7.root 
rm -rf BToKMuMu_Data_2012B_7.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_8.root ./
./sel data cut1 BToKMuMu_Data_2012B_8.root Data_2012B_bfAnti_cut1_8.root 
rm -rf BToKMuMu_Data_2012B_8.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_9.root ./
./sel data cut1 BToKMuMu_Data_2012B_9.root Data_2012B_bfAnti_cut1_9.root 
rm -rf BToKMuMu_Data_2012B_9.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_10.root ./
./sel data cut1 BToKMuMu_Data_2012B_10.root Data_2012B_bfAnti_cut1_10.root 
rm -rf BToKMuMu_Data_2012B_10.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_11.root ./
./sel data cut1 BToKMuMu_Data_2012B_11.root Data_2012B_bfAnti_cut1_11.root 
rm -rf BToKMuMu_Data_2012B_11.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_12.root ./
./sel data cut1 BToKMuMu_Data_2012B_12.root Data_2012B_bfAnti_cut1_12.root 
rm -rf BToKMuMu_Data_2012B_12.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_13.root ./
./sel data cut1 BToKMuMu_Data_2012B_13.root Data_2012B_bfAnti_cut1_13.root 
rm -rf BToKMuMu_Data_2012B_13.root

eos cp /eos/cms/store/user/gechen/crab3_run/MuOniaParked/BuToKMuMu_Data_2012B_8TeV_v8/BToKMuMu_Data_2012B_14.root ./
./sel data cut1 BToKMuMu_Data_2012B_14.root Data_2012B_bfAnti_cut1_14.root 
rm -rf BToKMuMu_Data_2012B_14.root

hadd Data_2012B_bfAnti_cut1.root Data_2012B_bfAnti_cut1_*
rm -rf Data_2012B_bfAnti_cut1_*


