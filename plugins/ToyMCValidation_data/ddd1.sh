

python scan.py ToyVali_7
rm -rf AllCuts_v3/2D_7_* 
rm -rf plots/*.pdf

python scan.py ToyVali_8
rm -rf AllCuts_v3/2D_8_* 
rm -rf plots/*.pdf

python scan.py ToyVali_9
rm -rf AllCuts_v3/2D_9_* 
rm -rf plots/*.pdf

python scan.py ToyVali_10
rm -rf AllCuts_v3/2D_10_* 
rm -rf plots/*.pdf

source fcn_data.sh 
source ReFit.sh
