
#./fit createValiToys test 0 >& log.log
#echo -e "bin 0 done\n"
#./fit createValiToys test 1 >& log.log
#echo -e "bin 1 done\n"
#./fit createValiToys test 2 >& log.log
#echo -e "bin 2 done\n"
#./fit createValiToys test 4 >& log.log
#echo -e "bin 4 done\n"
#./fit createValiToys test 6 >& log.log
#echo -e "bin 6 done\n"
#./fit createValiToys test 7 >& log.log
#echo -e "bin 7 done\n"
#./fit createValiToys test 8 >& log.log
#echo -e "bin 8 done\n"
#./fit createValiToys test 9 >& log.log
#echo -e "bin 9 done\n"
./fit createValiToys test 10 >& log.log
echo -e "bin 10 done\n"
python fitFCToys_unCons.py 10
