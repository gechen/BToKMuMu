

bsub -R "pool>3000" -q 2nd -J bin1 < ssss1.sh 

./fit createFCToys test 0 

./fit PlotAfbFh_f test 7 angular2D_NLL
python scan_7.py angular2D_NLL 7


python scan_8.py angular2D_Profiled_Afb 8
./fit PlotFCN_Profiled test 8 angular2D_Profiled_Afb
python scan_8.py angular2D_Profiled_Afb_refit


