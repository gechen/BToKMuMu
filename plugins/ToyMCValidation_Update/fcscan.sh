
python fitFCToys_1.py 0 >& log.log
echo -e "bin 0 done\n"
python fitFCToys_1.py 1 >& log.log
echo -e "bin 1 done\n"
python fitFCToys_1.py 2 >& log.log
echo -e "bin 2 done\n"
python fitFCToys_1.py 4 >& log.log
echo -e "bin 4 done\n"
python fitFCToys_1.py 6 >& log.log
echo -e "bin 6 done\n"
python fitFCToys_1.py 7 >& log.log
echo -e "bin 7 done\n"
python fitFCToys_1.py 8 >& log.log
echo -e "bin 8 done\n"

./fit FCScan test 0 >& log.log
echo -e "bin 0 done\n"
./fit FCScan test 1 >& log.log
echo -e "bin 1 done\n"
./fit FCScan test 2 >& log.log
echo -e "bin 2 done\n"
./fit FCScan test 4 >& log.log
echo -e "bin 4 done\n"
./fit FCScan test 6 >& log.log
echo -e "bin 6 done\n"
./fit FCScan test 7 >& log.log
echo -e "bin 7 done\n"
./fit FCScan test 8 >& log.log
echo -e "bin 8 done\n"



