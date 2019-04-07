
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 999 > AllCuts_v3/a999.log 

#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 0 > AllCuts_v3/a0.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 1 > AllCuts_v3/a1.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 2 > AllCuts_v3/a2.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 3 >& AllCuts_v3/a3.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 4 > AllCuts_v3/a4.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 5 >& AllCuts_v3/a5.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 6 > AllCuts_v3/a6.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 7 > AllCuts_v3/a7.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 8 > AllCuts_v3/a8.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 9 > AllCuts_v3/a9.log 
#./fit accXrecoEffLimited ../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3.root 10 > AllCuts_v3/a10.log 
python scan_0.py accXrecoEffLimited
python scan_1.py accXrecoEffLimited
python scan_2.py accXrecoEffLimited
python scan_4.py accXrecoEffLimited
python scan_6.py accXrecoEffLimited
python scan_7.py accXrecoEffLimited
python scan_8.py accXrecoEffLimited
python scan_9.py accXrecoEffLimited
python scan_10.py accXrecoEffLimited
