#define MCB2KMuMu_cxx
#include "MCB2KMuMu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>

void MCB2KMuMu:: ClearEvent() {
	Triggers            = -999;
	DimuonPt            = -999;
	Mumumass            = -999;
	Mumumasserr         = -999;
	TkChg               = -999;
	Trkdcasigbs         = -999;
	Bmass               = -999;
	Bchg                = -999;
	Bvtxcl              = -999;
	Blxysig             = -999;
	Bcosalphabs         = -999;
	Bcosalphabs2D       = -999;
	Bctau               = -999;
	Q2                  = -999;
	CosThetaL           = -999;
	
}
void MCB2KMuMu::Loop()
{
	TDatime t_begin_ , t_now_ ;
	t_begin_.Set();
	printf("\n ---------- Begin Job ---------- \n \n");
	t_begin_.Print();
	cout<<endl;

//   In a ROOT session, you can do:
//      Root > .L MCB2KMuMu.C
//      Root > MCB2KMuMu t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<nentries<<endl;
	
//	TFile *newfile = new TFile("./MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT+Anti3_NoGEN.root","recreate");

//	TFile *newfile = new TFile("./MC_Kst0MuMu_8TeV_v4_cut0+HLT+Q+B1+resonance-1+OPT.root","recreate");
//	TFile *newfile = new TFile("./MC_JpsiX_8TeV_v4_cut0+HLT+Q+B2+resonance-1.root","recreate");
//	TFile *newfile = new TFile("./MC_Psi2SK_8TeV_v4_cut0+HLT+Q+B2+resonance-1.root","recreate");
//	TFile *newfile = new TFile("./MC_JpsiK_8TeV_v4_cut0+HLT+Q+B2+resonance-1.root","recreate");
//	TFile *newfile = new TFile("./MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1_NoGEN.root","recreate");
	
//	TFile *newfile = new TFile("./Data_2012A_8TeV_v5_cut0+HLT+B2+resonance-1+OPT+Anti3.root","recreate");
//	TFile *newfile = new TFile("./Data_2012A_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT+Anti3.root","recreate");
//	TFile *newfile = new TFile("./Data_2012A_8TeV_v8_cut0+HLT+Q+B2+resonance-1+OPT+Anti3.root","recreate");
//	TFile *newfile = new TFile("./Data_2012B_8TeV_v8_cut0+HLT+Q+B2+resonance-1+OPT+Anti3.root","recreate");
//	TFile *newfile = new TFile("./Data_2012B_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT+Anti3.root","recreate");
//	TFile *newfile = new TFile("./Data_2012C_8TeV_v8_cut0+HLT+Q+B2+resonance-1+OPT+Anti3_2.root","recreate");
//	TFile *newfile = new TFile("./Data_2012C_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT+Anti3_2.root","recreate");
//	TFile *newfile = new TFile("./Data_2012D_8TeV_v8_cut0+HLT+Q+B2+resonance-1+OPT+Anti3_2.root","recreate");
//	TFile *newfile = new TFile("./Data_2012D_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT+Anti3_4.root","recreate");

	TTree *tree_ = fChain->CloneTree(0);

	Long64_t a=0, gen=0;
	int q0=0, q1=0, q2=0, q3=0, q4=0;
	int r1=0, r2=0;
	int p1=0, p2=0, p3=0, p4=0, p5=0;
	int A0=0, A1=0, A2=0, A3=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		GetEntry(jentry);
		if (Q2 == -999) {
			ClearEvent();
		//	tree_->Fill();
		//	gen++;
			continue;
		}
		else {
			if (Triggers != 1) {
			//	ClearEvent();
			//	tree_->Fill();
				continue;
			}
			q0+=1;
/////////////////////////// quality cuts   	////////////////////////	
		if (B_4vec->Pt() < 8.0) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		q1+=1;
		if (fabs(B_4vec->Eta()) > 2.2) {
		//	ClearEvent();
		//	tree_->Fill();	
			continue;
		}
		q2+=1;
/*		if (Bctau <1e-5) {
	//		ClearEvent();
	//		tree_->Fill();
			continue;
		}
*/		q3+=1;
	//	if (Bmass > 5.56 || Bmass <5.0) {
		if (Bmass > 5.60 || Bmass <5.10) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		q4+=1;
/////////////////////////// resonance cuts   	////////////////////////	
		if (Mumumass > 3.096916-5. *Mumumasserr && Mumumass < 3.096916+3. *Mumumasserr) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		r1+=1;
		if (Mumumass > 3.686109-3. *Mumumasserr && Mumumass < 3.686109+3. *Mumumasserr) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		r2+=1;
/////////////////////////// O.P.T. cuts   	////////////////////////	
		if (Bvtxcl < 0.08) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		p1+=1;
		if (Bcosalphabs2D < 0.9997) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		p2+=1;
		if (Trkdcasigbs < 2.7) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		p3+=1;
		if (Blxysig < 7.0) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		p4+=1;
		if (Tk_4vec->Pt() < 1.3) {
		//	ClearEvent();
		//	tree_->Fill();
			continue;
		}
		p5+=1;
	//	if (Trkdcasigbs == -999) continue;
		
/////////////////////////// anti-radiation	////////////////////////	
/*	// LHCb	
 		if (Bmass < 5.23) {
			if (Mumumass > 2.997 &&  Mumumass < 3.15) continue;
			if (Mumumass > 3.566 &&  Mumumass < 3.75) continue;

		}
*/		
/*	// Anti1	
		if (Bmass < 5.23) {
			A0+=1;
			if (fabs(Bmass-Mumumass-2.182) <0.140) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
			A1+=1;
			if (fabs(Bmass-Mumumass-1.593) <0.060) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
		}
*/
/*		if (Bmass < 5.23) {
			A0+=1;
			if (fabs(Bmass-Mumumass-2.182) <0.140) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
			A1+=1;
			if (fabs(Bmass-Mumumass-1.593) <0.070) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
			A2+=1;
		} else if (Bmass > 5.33) {
			if (fabs(Bmass-Mumumass-2.182) <0.115) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
			if (fabs(Bmass-Mumumass-1.593) <0.055) {
	//		ClearEvent();
	//		tree_->Fill();
				continue;
			}
		} else {
			A3+=1;
		}
*/		
	// Anti3	
		if (fabs(Bmass-Mumumass-2.182) <0.125) {
	//	ClearEvent();
	//	tree_->Fill();
			continue;
		}
		A0+=1;
		if (fabs(Bmass-Mumumass-1.593) <0.070) {
	//	ClearEvent();
	//	tree_->Fill();
			continue;
		}
		A1+=1;
		
		tree_->Fill();
		a++;

		}
   }
	tree_->AutoSave();
	cout<<"q0 = "<<q0<<endl;
	cout<<"q1 = "<<q1<<endl;
	cout<<"q2 = "<<q2<<endl;
	cout<<"q3 = "<<q3<<endl;
	cout<<"q4 = "<<q4<<endl;
	cout<<"r1 = "<<r1<<endl;
	cout<<"r2 = "<<r2<<endl;
	cout<<"p1 = "<<p1<<endl;
	cout<<"p2 = "<<p2<<endl;
	cout<<"p3 = "<<p3<<endl;
	cout<<"p4 = "<<p4<<endl;
	cout<<"p5 = "<<p5<<endl;
	cout<<"A0 = "<<A0<<endl;
	cout<<"A1 = "<<A1<<endl;
	cout<<"A2 = "<<A2<<endl;
	cout<<"A3 = "<<A3<<endl;
	cout<<"A1+A3 = "<<A1+A3<<endl;
	cout<<"A2+A3 = "<<A2+A3<<endl;
	cout<<endl;
	cout<<"a   = "<<a<<endl;
	cout<<"gen = "<<gen<<endl;



	t_now_.Set();
	printf(" \n \n ---------- End Job ---------- \n" ) ;
	t_now_.Print();
	printf(" duration: %i sec \n \n",
	t_now_.Convert() - t_begin_.Convert() );

}
