//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  8 05:05:07 2014 by ROOT version 5.34/14
// from TTree tree/tree
// found on file: BToKMuMu_SignalMC_8TeV_v3-4_cut0.root
//////////////////////////////////////////////////////////

#ifndef MCB2KMuMu_h
#define MCB2KMuMu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class MCB2KMuMu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Triggers;
	TLorentzVector  *Mum_4vec;
   TLorentzVector  *Mup_4vec;
   Double_t        DimuonPt;
   Double_t        Mumumass;
   Double_t        Mumumasserr;
   TLorentzVector  *Tk_4vec;
   Int_t           TkChg;
   Double_t        Trkdcasigbs;
   TLorentzVector  *B_4vec;
   Double_t        Bmass;
   Int_t           Bchg;
   Double_t        Bvtxcl;
   Double_t        Blxysig;
   Double_t        Bcosalphabs;
   Double_t        Bcosalphabs2D;
   Double_t        Bctau;
   Double_t        Q2;
   Double_t        CosThetaL;
   
	Int_t           genBChg;
   Double_t        genBPt;
   Double_t        genBEta;
   Double_t        genBPhi;
   Double_t        genDimuonPt;
   Double_t        genMupPt;
   Double_t        genMupEta;
   Double_t        genMumPt;
   Double_t        genMumEta;
   Int_t           genTkChg;
   Double_t        genTkPt;
   Double_t        genTkEta;
   Double_t        genQ2;
   Double_t        genCosThetaL;

   // List of branches
   TBranch        *b_Triggers;   //!
   TBranch        *b_Mum_4vec;   //!
   TBranch        *b_Mup_4vec;   //!
   TBranch        *b_DimuonPt;   //!
   TBranch        *b_Mumumass;   //!
   TBranch        *b_Mumumasserr;   //!
   TBranch        *b_Tk_4vec;   //!
   TBranch        *b_TkChg;   //!
   TBranch        *b_Trkdcasigbs;   //!
   TBranch        *b_B_4vec;   //!
   TBranch        *b_Bmass;   //!
   TBranch        *b_Bchg;   //!
   TBranch        *b_Bvtxcl;   //!
   TBranch        *b_Blxysig;   //!
   TBranch        *b_Bcosalphabs;   //!
   TBranch        *b_Bcosalphabs2D;   //!
   TBranch        *b_Bctau;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_CosThetaL;   //!
   
	TBranch        *b_genBChg;   //!
   TBranch        *b_genBPt;   //!
   TBranch        *b_genBEta;   //!
   TBranch        *b_genBPhi;   //!
   TBranch        *b_genDimuonPt;   //!
   TBranch        *b_genMupPt;   //!
   TBranch        *b_genMupEta;   //!
   TBranch        *b_genMumPt;   //!
   TBranch        *b_genMumEta;   //!
   TBranch        *b_genTkChg;   //!
   TBranch        *b_genTkPt;   //!
   TBranch        *b_genTkEta;   //!
   TBranch        *b_genQ2;   //!
   TBranch        *b_genCosThetaL;   //!

   MCB2KMuMu(TTree *tree=0);
   virtual ~MCB2KMuMu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
	void SaveGen();
	void ClearEvent();
};

#endif

#ifdef MCB2KMuMu_cxx
MCB2KMuMu::MCB2KMuMu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
	   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../RootFiles/bfOPT/BToKMuMu_SignalMC_8TeV_v4_bfOPT_cut1.root");
	   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./Data_2012_bfOPT_cut1.root");
	   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../RootFiles/OPT/Data_2012_bfAnti_cut1.root");
	   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../RootFiles/OPT/Data_2012_bfAnti_cutJpsi.root");
	  // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../RootFiles/OPT/Data_2012_bfAnti_cutPsip.root");
      if (!f || !f->IsOpen()) {
	      f = new TFile("../../RootFiles/bfOPT/BToKMuMu_SignalMC_8TeV_v4_bfOPT_cut1.root");
	      //f = new TFile("./Data_2012_bfOPT_cut1.root");
	      //f = new TFile("../../RootFiles/OPT/Data_2012_bfAnti_cut1.root");
	      //f = new TFile("../../RootFiles/OPT/Data_2012_bfAnti_cutJpsi.root");
	      //f = new TFile("../../RootFiles/OPT/Data_2012_bfAnti_cutPsip.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

MCB2KMuMu::~MCB2KMuMu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MCB2KMuMu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MCB2KMuMu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MCB2KMuMu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Mum_4vec = 0;
   Mup_4vec = 0;
   Tk_4vec = 0;
   B_4vec = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
   fChain->SetBranchAddress("Mum_4vec", &Mum_4vec, &b_Mum_4vec);
   fChain->SetBranchAddress("Mup_4vec", &Mup_4vec, &b_Mup_4vec);
   fChain->SetBranchAddress("DimuonPt", &DimuonPt, &b_DimuonPt);
   fChain->SetBranchAddress("Mumumass", &Mumumass, &b_Mumumass);
   fChain->SetBranchAddress("Mumumasserr", &Mumumasserr, &b_Mumumasserr);
   fChain->SetBranchAddress("Tk_4vec", &Tk_4vec, &b_Tk_4vec);
   fChain->SetBranchAddress("TkChg", &TkChg, &b_TkChg);
   fChain->SetBranchAddress("Trkdcasigbs", &Trkdcasigbs, &b_Trkdcasigbs);
   fChain->SetBranchAddress("B_4vec", &B_4vec, &b_B_4vec);
   fChain->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
   fChain->SetBranchAddress("Bchg", &Bchg, &b_Bchg);
   fChain->SetBranchAddress("Bvtxcl", &Bvtxcl, &b_Bvtxcl);
   fChain->SetBranchAddress("Blxysig", &Blxysig, &b_Blxysig);
   fChain->SetBranchAddress("Bcosalphabs", &Bcosalphabs, &b_Bcosalphabs);
   fChain->SetBranchAddress("Bcosalphabs2D", &Bcosalphabs2D, &b_Bcosalphabs2D);
   fChain->SetBranchAddress("Bctau", &Bctau, &b_Bctau);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("CosThetaL", &CosThetaL, &b_CosThetaL);
   
   fChain->SetBranchAddress("genBChg", &genBChg, &b_genBChg);
   fChain->SetBranchAddress("genBPt", &genBPt, &b_genBPt);
   fChain->SetBranchAddress("genBEta", &genBEta, &b_genBEta);
   fChain->SetBranchAddress("genBPhi", &genBPhi, &b_genBPhi);
   fChain->SetBranchAddress("genDimuonPt", &genDimuonPt, &b_genDimuonPt);
   fChain->SetBranchAddress("genMupPt", &genMupPt, &b_genMupPt);
   fChain->SetBranchAddress("genMupEta", &genMupEta, &b_genMupEta);
   fChain->SetBranchAddress("genMumPt", &genMumPt, &b_genMumPt);
   fChain->SetBranchAddress("genMumEta", &genMumEta, &b_genMumEta);
   fChain->SetBranchAddress("genTkChg", &genTkChg, &b_genTkChg);
   fChain->SetBranchAddress("genTkPt", &genTkPt, &b_genTkPt);
   fChain->SetBranchAddress("genTkEta", &genTkEta, &b_genTkEta);
   fChain->SetBranchAddress("genQ2", &genQ2, &b_genQ2);
   fChain->SetBranchAddress("genCosThetaL", &genCosThetaL, &b_genCosThetaL);
  
	Notify();
}

Bool_t MCB2KMuMu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MCB2KMuMu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MCB2KMuMu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MCB2KMuMu_cxx
