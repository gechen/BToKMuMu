
eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_1.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_1.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_1.root
rm -rf BuToPSi2SK_MC_8TeV_v4_1.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_2.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_2.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_2.root
rm -rf BuToPSi2SK_MC_8TeV_v4_2.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_3.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_3.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_3.root
rm -rf BuToPSi2SK_MC_8TeV_v4_3.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_4.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_4.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_4.root
rm -rf BuToPSi2SK_MC_8TeV_v4_4.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_5.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_5.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_5.root
rm -rf BuToPSi2SK_MC_8TeV_v4_5.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_6.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_6.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_6.root
rm -rf BuToPSi2SK_MC_8TeV_v4_6.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_7.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_7.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_7.root
rm -rf BuToPSi2SK_MC_8TeV_v4_7.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_8.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_8.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_8.root
rm -rf BuToPSi2SK_MC_8TeV_v4_8.root

eos cp /eos/cms/store/user/gechen/BuToPsi2SK_MC_8TeV_v4/BuToPSi2SK_MC_8TeV_v4_9.root ./
./sel mc.lite cutPsip ./BuToPSi2SK_MC_8TeV_v4_9.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_9.root
rm -rf BuToPSi2SK_MC_8TeV_v4_9.root

hadd BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip.root BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_*
rm -rf BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip_*



