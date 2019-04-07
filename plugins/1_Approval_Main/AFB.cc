
#include <TString.h>
#include <fstream>

#include<iostream>

using namespace std;



char LHCb_files[19][200] ={ "q2_0_to_1.dat",
							"q2_1_to_2.dat",
							"q2_2_to_3.dat",
							"q2_3_to_4.dat",
							"q2_4_to_5.dat",
							"q2_5_to_6.dat",
							"q2_6_to_7.dat",
							"q2_7_to_8.dat",
							"q2_11_to_11.dat",
							"q2_11_to_12.dat",
							"q2_15_to_16.dat",
							"q2_16_to_17.dat",
							"q2_17_to_18.dat",
							"q2_18_to_19.dat",
							"q2_19_to_20.dat",
							"q2_20_to_21.dat",
							"q2_21_to_22.dat",
							"q2_1_to_6.dat",
							"q2_15_to_22.dat"};

void AFB() {
	for (int iibin = 0; iibin < 19; iibin++) {
		ifstream inputfile(TString::Format("./LHCb_data/%s",LHCb_files[iibin]));
		string s;
		string ss1, ss2;
		int number = 0;
		int number_max = 0;
		float afb_max = -1, afb_min = 1, fh_max = -1, fh_min = 3;
		float a[2000], b[2000], c[2000];
		while (getline( inputfile, s ))
		{
			if (ss2 != "AFB") {
				inputfile >> ss1 >> ss2;
				continue;
			} else {
				inputfile >> a[number] >> b[number] >> c[number];
				number++;
			}
		}
		float pvalue = 0;
		for (int i = 0; i<number; i++) {
			if (pvalue < c[i]) {
				pvalue = c[i];
				number_max = i;
			}
			if (c[i] > 0.6) {
				if (afb_max < a[i]) afb_max = a[i];
				if (afb_min > a[i]) afb_min = a[i];
				if (fh_max < b[i]) fh_max = b[i];
				if (fh_min > b[i]) fh_min = b[i];
			}
		}
		cout << LHCb_files[iibin] << " " << number_max << ": " << a[number_max] << " " << b[number_max] << " " << c[number_max] << endl;
		cout << "afb_min = " << afb_min << "; afb_max = " << afb_max << endl;
		cout << "fh_min  = " << fh_min << "; fh_max  = " << fh_max << endl;
	}
}
