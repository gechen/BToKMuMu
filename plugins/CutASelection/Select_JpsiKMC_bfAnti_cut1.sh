
eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_0.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_0.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_0.root
rm -rf BuToJPSiK_MC_8TeV_v4_0.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_1.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_1.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_1.root
rm -rf BuToJPSiK_MC_8TeV_v4_1.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_2.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_2.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_2.root
rm -rf BuToJPSiK_MC_8TeV_v4_2.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_3.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_3.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_3.root
rm -rf BuToJPSiK_MC_8TeV_v4_3.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_4.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_4.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_4.root
rm -rf BuToJPSiK_MC_8TeV_v4_4.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_5.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_5.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_5.root
rm -rf BuToJPSiK_MC_8TeV_v4_5.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_6.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_6.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_6.root
rm -rf BuToJPSiK_MC_8TeV_v4_6.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_7.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_7.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_7.root
rm -rf BuToJPSiK_MC_8TeV_v4_7.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_8.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_8.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_8.root
rm -rf BuToJPSiK_MC_8TeV_v4_8.root

eos cp /eos/cms/store/user/gechen/BuToJPSiK_MC_8TeV_v4/BuToJPSiK_MC_8TeV_v4_9.root ./
./sel mc.lite cut1 ./BuToJPSiK_MC_8TeV_v4_9.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_9.root
rm -rf BuToJPSiK_MC_8TeV_v4_9.root

hadd BToJpsiK_MC_8TeV_v4_bfAnti_cut1.root BToJpsiK_MC_8TeV_v4_bfAnti_cut1_*
rm -rf BToJpsiK_MC_8TeV_v4_bfAnti_cut1_*


