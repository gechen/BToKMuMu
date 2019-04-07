#include "cut_optimize.h"
#include "TFile.h"
#include <cmath>

int main() {
//	gROOT->ProcessLine(".L ./cut_optimize.h++");
//  TFile *fsigbg = new TFile("../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT.root");
//  TFile *fsig = new TFile("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT_NoGEN.root");
//  TFile *fsigbg = new TFile("../RootFiles/OPT/Data_2012_8TeV_v4_OPTwVTXLoose.root");
//  TFile *fsig = new TFile("../RootFiles/OPT/MC_Signal_8TeV_v4_OPTwVTXLoose_NoGEN.root");
//  TFile *fsigbg = new TFile("../RootFiles/OPT/Data_2012_8TeV_v4_OPT3.root");
//  TFile *fsig = new TFile("../RootFiles/OPT/MC_Signal_8TeV_v4_OPT3_NoGEN.root");
//  TFile *fsigbg = new TFile("../RootFiles/OPT/Data_2012_bfAnti_cut1.root");
//  TFile *fsig = new TFile("../RootFiles/OPT/BToKMuMu_SignalMC_8TeV_v4_bfAnti_cut1_NoGEN.root");
  TFile *fsigbg = new TFile("../RootFiles/OPT/Data_2012_bfAnti_cutJpsi.root");
  TFile *fsig = new TFile("../RootFiles/OPT/BToJpsiK_MC_8TeV_v4_bfAnti_cutJpsi.root");
//  TFile *fsigbg = new TFile("../RootFiles/OPT/Data_2012_bfAnti_cutPsip.root");
//  TFile *fsig = new TFile("../RootFiles/OPT/BToPsi2SK_MC_8TeV_v4_bfAnti_cutPsip.root");
	// LHCb
//	TFile *fsigbg = new TFile("../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root");
//	TFile *fsig = new TFile("../RootFiles/Files/MC_JpsiK_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root");
//	TFile *fsigbg = new TFile("../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root");
//	TFile *fsig = new TFile("../RootFiles/Files/MC_Psi2SK_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root");
	
	TTree *t_sig = (TTree*)fsig->Get("tree");
	TTree *t_sigbg = (TTree*)fsigbg->Get("tree");
	
	CutOptimize m_cutopt;
	m_cutopt.set_debug(1);
	
	m_cutopt.set_tree_sig(t_sig, 0.006213575);
	m_cutopt.set_tree_sigbg(t_sigbg);
//	add all possible cut options
//	m_cutopt.add_opt(true, "(Mumumass < 3.0969) && abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.14
//	m_cutopt.add_opt(true, "(Mumumass < 3.0969) && abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.01
//
//	m_cutopt.add_opt(true, "(Mumumass < 3.6861) && (Mumumass > 3.0969) && abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.07
//	m_cutopt.add_opt(true, "(Mumumass < 3.6861) && (Mumumass > 3.0969) && abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.07
//
//	m_cutopt.add_opt(true, "(Mumumass > 3.6861) && abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.01
//	m_cutopt.add_opt(true, "(Mumumass > 3.6861) && abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.07

//	Anti-radiation  
//	m_cutopt.add_opt(true, "(Bmass < 5.23) &&  abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.138
//	m_cutopt.add_opt(true, "(Bmass < 5.23) &&  abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005);  // 0.056
//
//	m_cutopt.add_opt(true, "(Bmass > 5.33) &&  abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.09
//	m_cutopt.add_opt(true, "(Bmass > 5.33) &&  abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.06
//
//	m_cutopt.add_opt(true, "(Bmass < 5.33) && (Bmass > 5.23) && abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.09
//	m_cutopt.add_opt(true, "(Bmass < 5.33) && (Bmass > 5.23) && abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.06

// Anti3
//	m_cutopt.add_opt(true, " abs(Bmass - Mumumass - 2.182) >", 0.01, 0.005);  // 0.09
//	m_cutopt.add_opt(true, " abs(Bmass - Mumumass - 1.593) >", 0.01, 0.005); // 0.06
// anti3 control
	m_cutopt.add_opt(true, " abs(Bmass - Mumumass - 2.182) <", 0.01, 0.005);  // 0.09
//	m_cutopt.add_opt(true, " abs(Bmass - Mumumass - 1.593) <", 0.01, 0.005); // 0.06
//	LHCb
//	m_cutopt.add_opt(true, "(Bmass < 5.23) && (Mumumass < 3.15) && (3.097 - Mumumass) <", 0.01, 0.005);  // 0.13
//	m_cutopt.add_opt(true, "(Bmass < 5.23) && (Mumumass < 3.75) && (3.686 - Mumumass) <", 0.01, 0.005); // 0.16
  
// Anti4
//	m_cutopt.add_opt(true, " - abs(Bmass - Mumumass - 2.182) > ", -0.005, 0.001);  // 0.09
//	m_cutopt.add_opt(true, " - abs(Bmass - Mumumass - 1.593) > ", -0.005, 0.001); // 0.06
  
  
  //m_cutopt.add_opt(true, "Bvtxcl >", 0.08, 0.01);
  //m_cutopt.add_opt(true, "Bcosalphabs2D >", 0.9995, 0.0001); 
  //m_cutopt.add_opt(true, "Blxysig >", 8.1, 0.1); 
  //m_cutopt.add_opt(true, "Blxysig >", 6.9, 0.1); 
  //m_cutopt.add_opt(true, "Bctau >", 0.4, 0.1); 
  //m_cutopt.add_opt(true, "abs", "Bmass", "<", 5.279, 0.062, 0.002);  
  //m_cutopt.add_opt(true, "Bpt >", 7.0, 0.2); 
  //m_cutopt.add_opt(true, "Tk_4vec.Pt() >", 1.6, 0.1); 
  //m_cutopt.add_opt(true, "Trkdcasigbs >", 2.0, 0.1); 


  cout << "after optimization" << endl;
  m_cutopt.print();

  m_cutopt.loop();
  cout << "============================" << endl;
  cout << "after anti-radiation" << endl;
  m_cutopt.print();
  return 0;
}
