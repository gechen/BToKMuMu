#include "cut_optimize.h"
#include "TFile.h"
#include <cmath>

int main() {
  //gROOT->ProcessLine(".L ./cut_optimize.h++");
  TFile *fsigbg = new TFile("../RootFiles/bfOPT/Data_2012_bfOPT_cut1.root");
  TFile *fsig = new TFile("../RootFiles/bfOPT/BToKMuMu_SignalMC_8TeV_v4_bfOPT_cut1.root");
  //TFile *fsigbg = new TFile("../RootFiles/bfOPT/Data_2012_8TeV_v4_OPTonlyVTXLoose.root");
  //TFile *fsig = new TFile("../RootFiles/bfOPT/MC_Signal_8TeV_v4_OPTonlyVTXLoose_NoGEN.root");
//  TFile *fsigbg = new TFile("../RootFiles/bfOPT/Data_2012_bfOPT_cutJpsi.root");
//  TFile *fsig = new TFile("../RootFiles/bfOPT/BToJpsiK_MC_8TeV_v4_bfOPT_cutJpsi.root");
  //TFile *fsigbg = new TFile("../RootFiles/bfOPT/Data_2012_bfOPT_cutPsip.root");
  //TFile *fsig = new TFile("../RootFiles/bfOPT/BToPsi2SK_MC_8TeV_v4_bfOPT_cutPsip.root");

  TTree *t_sig = (TTree*)fsig->Get("tree");
  TTree *t_sigbg = (TTree*)fsigbg->Get("tree");

  CutOptimize m_cutopt;
  m_cutopt.set_debug(1);

  m_cutopt.set_tree_sig(t_sig, 0.006213575);
  m_cutopt.set_tree_sigbg(t_sigbg);
  // add all possible cut options
  m_cutopt.add_opt(true, "Bvtxcl >", 0.06, 0.01);
  m_cutopt.add_opt(true, "Blxysig >", 4.0, 0.1); 
  m_cutopt.add_opt(true, "Bcosalphabs2D >", 0.9976, 0.0001); 
  m_cutopt.add_opt(true, "Trkdcasigbs >", 2.0, 0.1); 
  m_cutopt.add_opt(true, "Tk_4vec.Pt() >", 0.8, 0.1); 
/* // HLT   
  m_cutopt.add_opt(true, "Bvtxcl >", 0.08, 0.01);
  m_cutopt.add_opt(true, "Bcosalphabs2D >", 0.9997, 0.0001); 
  m_cutopt.add_opt(true, "Trkdcasigbs >", 2.7, 0.1); 
  m_cutopt.add_opt(true, "Blxysig >", 7.0, 0.1); 
  m_cutopt.add_opt(true, "Tk_4vec.Pt() >", 1.3, 0.1); 
*/  
/*
  m_cutopt.add_opt(true, "Bvtxcl >", 0.08, 0.01);
  m_cutopt.add_opt(true, "Bcosalphabs2D >", 0.9995, 0.0001); 
  m_cutopt.add_opt(true, "Blxysig >", 6.5, 0.1); 
  //m_cutopt.add_opt(true, "Blxysig >", 6.9, 0.1); 
  //m_cutopt.add_opt(true, "Bctau >", 0.4, 0.1); 
  //m_cutopt.add_opt(true, "abs", "Bmass", "<", 5.279, 0.062, 0.002);  
  //m_cutopt.add_opt(true, "Bpt >", 7.0, 0.2); 
  m_cutopt.add_opt(true, "Tk_4vec.Pt() >", 1.1, 0.1); 
  m_cutopt.add_opt(true, "Trkdcasigbs >", 2.5, 0.1); 
*/

  cout << "before optimization" << endl;
  m_cutopt.print();

  m_cutopt.loop();
  cout << "============================" << endl;
  cout << "after optimization" << endl;
  m_cutopt.print();
  return 0;
}
