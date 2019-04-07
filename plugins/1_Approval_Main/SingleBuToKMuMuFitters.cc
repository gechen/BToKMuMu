// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
//////////////////////////////////////////////////////////////////////////
#include "collections.cc"
#include "bmass.cc"
#include "angular_gen.cc"
#include "accXrecoEff.cc"
#include "angular_reco.cc"
#include "angular_1a1b2a.cc"
#include "prior.cc"
#include "angular2D.cc"
#include "systemUncertainties.cc"
#include "FCToyFitJobs_Profiled.cc"
int main(int argc, char** argv) {
//	Tags
	is7TeVCheck = false;   
//	Help message
	if (argc <= 2) {
	    cout<<"Usage       : ./fit Function infile binID"<<endl;
		cout<<"Functions   :"<<endl;
		cout<<"    0. bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg."<<endl;
		cout<<"    1. angular_gen         Derive F_{H} and A_{FB} from cosThetaL distribution at GEN level."<<endl;
		cout<<"       acceptance          Get acceptance map from unfiltered signal GEN, |Mu pT| > 2.8 GeV, |Mu eta| < 2.3."<<endl;
		cout<<"       recoEff             Get reconstruction efficiency map from signal simulation."<<endl;
		cout<<"    2. accXrecoEff         Get efficiency map from signal simulation."<<endl;
		cout<<"    3. angular_reco        Derive F_{H} and A_{FB} from cosThetaL distribution at RECO level."<<endl;
		cout<<"    4. angular_gen_R       Derive F_{H} and A_{FB} from gencosThetaL distribution at GEN level with official MC."<<endl;
		cout<<"    5. angular2D_1a_Sm     Leading step1 to angular2D, determine signal shape from simulation."<<endl;
		cout<<"    6. angular2D_1b_YpPm   Leading step2 to angular2D, determine mass spectrum of peaking bkg from simulation."<<endl;
		cout<<"    7. angular2D_2a_PkPl   Leading step3 to angular2D, determine angular dist. of peaking bkg from simulation."<<endl;
		cout<<"    8. angular2D_prior     Leading step4 to angular2D, fit to data sideband to get initial values of combinatorial bkg."<<endl;
		cout<<"    9. angular2D           Derive F_{H} and A_{FB} by fitting to mass and angular distribution."<<endl;
		cout<<"       For data fitting test:        ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
		cout<<"       For initial values scanning:          ./fit <function> <input.root> <iBin>  "<<endl;
		cout<<"   10. PlotFCN             Plot FCN distribution for each q^{2} bin, find the samllest FCN value"<<endl;
		cout<<"                              and it's initial value, save the fitting results."<<endl;
		cout<<"       For refit with final initial values:  ./fit <function> <input.root> <iBin> refit  "<<endl;
		cout<<"   11. PlotAfbFh           Plot Afb-Fh 2D distribution for each q^{2} bin, and the final results"<<endl;
		cout<<"Remark      :"<<endl;
		cout<<"    1. Outputs will be stored in ./plots, please keep the directory."<<endl;
		cout<<"    2. Fitted parameters will be stored in ./fitParameters/*.txt, please keep the directory."<<endl;
		cout<<"    3. Scaned fitted parameters will be stored in ./OutputValues/bin*/*.txt, please keep these directories."<<endl;
	//	cout<<"    3. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!"<<endl;
		return 0;
	}
//	main
  	if (argc != 4 && argc != 6 && argc != 7 && argc != 8 && argc != 5 && argc != 17 && argc != 14 && argc != 15){
		cout<<"./fit func infile binID"<<endl;
		for (int i = 0; i < 11; i++) {
		//	if (i == 3 || i == 5) continue;
			printf("    Bin %d : %s\n",i,Q2range[i]);
		}
		return 0;
	}
    int nToy=500;
	TString func    = argv[1];
	TString infile  = argv[2];
	int iBin        = atoi(argv[3]);
    if (iBin == 0) nToy = 1000;
	
	if (func == "MultivariateGaussianTest") {
		const char outfile[]="accXrecoEff";
		MultivariateGaussianTest(iBin, outfile);
	} else if (func == "DataMC_CosThetaL_Jpsi") {
        DataMC_CosThetaL_Jpsi(); 
	} else if (func == "DataMC_CosThetaL_Psi2S") {
        DataMC_CosThetaL_Psi2S(); 
	} else if (func == "PlotBlxsig_all") {
		ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        PlotBlxsig_all(iBin);  //// dimuon mass plots after all cuts
	} else if (func == "PlotBlxsig_CosThetaL") {
		ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        PlotBlxsig_CosThetaL(iBin);  //// dimuon mass plots after all cuts
	} else if (func == "PlotDimuon_all") {
		ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        PlotDimuon_all(iBin);  //// dimuon mass plots after all cuts
	} else if (func == "PlotDimuon") {
		PlotDimuon();  //// dimuon mass plots w/o reasonace cut and anti cuts.
	} else if (func == "singleMuonEff") {
		if (iBin < 0 || iBin > 9) return 0;
		const char outfile[]="singleMuonEff";
		singleMuonEff(iBin); 
	} else if (func == "singleMuondRDimuonPt") {
		if (iBin < 0 || iBin > 9) return 0;
		const char outfile[]="singleMuondRDimuonPt";
		singleMuondRDimuonPt(iBin); 
	} else if (func == "OPTEff") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="OPTEff";
		OPTEff(iBin, outfile); 
/////////////////////////////////////////////////////////////////////////////////////////////
	} else if (func == "bmass") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass";
		bmass(iBin, outfile); 
	} else if (func == "bmass_JpsiKGEN") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass_JpsiKGEN";
		bmass_JpsiKGEN(iBin, outfile); 
	} else if (func == "bmass_JpsiK") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass_JpsiK";
		bmass_JpsiK(iBin, outfile); 
	} else if (func == "bmass_Psi2SK") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass_Psi2SK";
		bmass_Psi2SK(iBin, outfile); 
/////////////////////////////////////////////////////////////////////////////////////////////
	} else if (func == "accXrecoEff") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEff(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else if (iBin == 998) {
			Accptance(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "recoEffQ2dR") {
		if (iBin >= 0 && iBin <= 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			recoEffQ2dR(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "recoEffdRDimuonPt") {
		if (iBin >= 0 && iBin <= 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			recoEffdRDimuonPt(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "recoEffDimuonPt") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			recoEffDimuonPt(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "accXrecoEffFIT") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEffFIT(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate gen efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "accXrecoEffSYS") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEffSYS(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate gen efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "accXrecoEffGEN") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEffGEN(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate gen efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "accXrecoEffDimuon") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			accXrecoEffDimuon(iBin);
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate gen efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "accXrecoEffLimited") {
		if (iBin >= 0 && iBin < 11) {
			int Index;
			double EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10; 
			EffP1 = atof(argv[4]);
			EffP2 = atof(argv[5]);
			EffP3 = atof(argv[6]);
			EffP4 = atof(argv[7]);
			EffP5 = atof(argv[8]);
			EffP6 = atof(argv[9]);
			EffP7 = atof(argv[10]);
			EffP8 = atof(argv[11]);
			EffP9 = atof(argv[12]);
			EffP10 = atof(argv[13]);
			Index = atoi(argv[14]);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			cout<<"Index "<<Index<<" = 1: "<<EffP1<<"; 2: "<<EffP2<<"; 3: "<<EffP3<<"; 4: "<<EffP4<<"; 5: "<<EffP5<<"; 6: "<<EffP6<<"; 7: "<<EffP7<<"; 8: "<<EffP8<<"; 9: "<<EffP9<<"; 10: "<<EffP10<<endl;
			accXrecoEffLimited(iBin, EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10, Index );
		//	accXrecoEffLimited(iBin, EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10 );
		}else if (iBin == 999) {
			createAccptanceHist(); // Set to true if no given ./RootFiles/acceptance_8TeV.root
		}else { 
			cout<<"Calculate gen efficiency, iBin counts from 0 to 10; or 999 for acceptance!"<<endl;
			return 0; 
		}
	} else if (func == "DimuonResolution") {
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<endl;
			DimuonResolution(iBin);
		}else { 
			return 0; 
		}
/////////////////////////////////////////////////////////////////////////////////////////////
	} else if (func == "angular2D_1a_Sm" || func == "angular2D_prior" || func == "angular2D_priorSYS" || func == "angular2D_prior_Check"){
		void (*fx)(int, const char*, bool);
		if ( func == "angular2D_1a_Sm" ){
			fx = angular2D_1a_Sm;
		}else if ( func == "angular2D_priorSYS" ){
			fx = angular2D_priorSYS;
		}else if ( func == "angular2D_prior" ){
			fx = angular2D_prior;
		}else if ( func == "angular2D_prior_Check" ){
			fx = angular2D_prior_Check;
		}	
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		fx(iBin,func,true);// By default overwrite exist parameters.
	} else if ( func == "angular2D_1b_YpPm_Jpsi" || func == "angular2D_1b_YpPm_Psi" || func == "angular2D_1b_YpPm_JX" || func == "angular2D_2a_PkPl_Jpsi" || func == "angular2D_2a_PkPl_Psi" || func == "angular2D_2a_PkPl_JX" ) {
		void (*fx)(int, const char*, bool);
		if (func == "angular2D_1b_YpPm_Jpsi" || func == "angular2D_1b_YpPm_Psi"|| func == "angular2D_1b_YpPm_JX" ) {
			fx = angular2D_1b_YpPm;
		}else if (func == "angular2D_2a_PkPl_Jpsi") {
			fx = angular2D_2a_PkPl_Jpsi;
		}else if (func == "angular2D_2a_PkPl_Psi") {
			fx = angular2D_2a_PkPl_Psi;
		}else if (func == "angular2D_2a_PkPl_JX") {
			fx = angular2D_2a_PkPl_JX;
		}
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		fx(iBin,func,true);// By default overwrite exist parameters.
	} else if (func == "angular_gen" || func == "angular2D" || func == "angular_gen_R" || func == "angular_reco" || func == "angular_reco_R" || func == "angular_reco_D" || func == "angular2D_Limited"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			
			float Iafb, Ifh, iafb, ifh;
			int Index = -999;
			if (argc == 7) {   //////////////////////////////   scanning   ///////////////////////////////////////////
				Index  = atoi(argv[4]);
				Iafb   = atof(argv[5]);
				Ifh    = atof(argv[6]);
				cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
				if (func == "angular2D") {
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_gen") {
					vbin = angular_gen_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_gen_R") {
					vbin = angular_gen_R_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_reco_R") {
					vbin = angular_reco_R_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_reco_D") {
					vbin = angular_reco_D_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_reco") {
					vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular2D_Limited") {
					double EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10;
					for (int idex = 1; idex <= 200; idex++) {
						EffP1 = readEffP(iBin, TString::Format("EffP_%d",idex), 0);
						EffP2 = readEffP(iBin, TString::Format("EffP_%d",idex), 1);
						EffP3 = readEffP(iBin, TString::Format("EffP_%d",idex), 2);
						EffP4 = readEffP(iBin, TString::Format("EffP_%d",idex), 3);
						EffP5 = readEffP(iBin, TString::Format("EffP_%d",idex), 4);
						EffP6 = readEffP(iBin, TString::Format("EffP_%d",idex), 5);
						EffP7 = readEffP(iBin, TString::Format("EffP_%d",idex), 6);
						EffP8 = readEffP(iBin, TString::Format("EffP_%d",idex), 7);
						EffP9 = readEffP(iBin, TString::Format("EffP_%d",idex), 8);
						EffP10 = readEffP(iBin, TString::Format("EffP_%d",idex),9);
						cout<<"1: "<<EffP1<<"; 2: "<<EffP2<<"; 3: "<<EffP3<<"; 4: "<<EffP4<<"; 5: "<<EffP5<<"; 6: "<<EffP6<<"; 7: "<<EffP7<<"; 8: "<<EffP8<<"; 9: "<<EffP9<<"; 10: "<<EffP10<<endl;
						vbin = angular2D_Limited_bin(iBin, Iafb, Ifh, Index, EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10, idex );
					}
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				} else {
				//	cout<<"A test of data fitting with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
					return 0;
				}
			} else if (argc == 4) {				
			//	TF1 *f1 = new TF1("f1","gaus",-1,1);
				TF1 *f1 = new TF1("f1","x+1",-1000,1000);
				TF1 *f2 = new TF1("f2","gaus",-10,10);
				f2->SetParameters(1,0,0.1);
				int NIndex = 1000;
				if (func == "angular_gen") NIndex = 50;
				else if (func == "angular_gen_R") NIndex = 200;
				else if (func == "angular_reco") NIndex = 200;
				else if (func == "angular_reco_R") NIndex = 200;
				else if (func == "angular_reco_D") NIndex = 200;
				else if (func == "angular2D_Limited") NIndex = 200;
				else if (func == "angular2D") {
					if (iBin == 0 || iBin == 6) NIndex = 200;
					else if (iBin == 8 || iBin == 10) NIndex = 600;
					else if (iBin == 7) NIndex = 1000;
					else NIndex = 400;
				}
				for (int Index = 0; Index < NIndex; Index+=1) {    
					Ifh = f1->GetRandom();
					Iafb = f2->GetRandom();
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
					if (func == "angular2D") {
						vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_gen") {
						vbin = angular_gen_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_gen_R") {
						vbin = angular_gen_R_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco_R") {
						vbin = angular_reco_R_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco_D") {
						vbin = angular_reco_D_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco") {
						vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular2D_Limited") {
						vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
					}
				}
			} else if (argc == 5) {  //////////////////////////////   refit ////////////////////////////////////
				TString type  = argv[4];
				if ( type == "refit") {
					Index = -1;
					if (func == "angular2D") {
						const char outfile[] = "angular2D";
					//	Iafb = readParam(iBin, "Iafb", 0);
					//	Ifh  = readParam(iBin, "Ifh", 0);
					//	cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
					
					// Cconstrained!!!	
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);   // Fitted results!
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
					//	iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 1);  // Selected fitted values!
					//	ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 1);
//	bin 0,1,6,7,8,9
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);  // Selected initial values!
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
					// Unconstrained!!!	
//	bin 2,4,10
////						Iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);  // Selected initial values!  Unconstrained
////						Ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
					//	Iafb = tan( iafb * TMath::Pi() / tan( (ifh - 1.5) * TMath::Pi() / 3.) ) ; // Selected initial values!  Unconstrained
					//	Ifh  = tan( (ifh - 1.5) * TMath::Pi() / 3. );
					//	Iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 1);  // Selected fitted values!   Unconstrained
					//	Ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 1);
						vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular2D_Limited") {
						const char outfile[] = "angular2D_Limited";
						double EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10; 
    					for (int idex = 1; idex <= 200; idex++) {
    						EffP1 = readEffP(iBin, TString::Format("EffP_%d",idex), 0);
    						EffP2 = readEffP(iBin, TString::Format("EffP_%d",idex), 1);
    						EffP3 = readEffP(iBin, TString::Format("EffP_%d",idex), 2);
    						EffP4 = readEffP(iBin, TString::Format("EffP_%d",idex), 3);
    						EffP5 = readEffP(iBin, TString::Format("EffP_%d",idex), 4);
    						EffP6 = readEffP(iBin, TString::Format("EffP_%d",idex), 5);
    						EffP7 = readEffP(iBin, TString::Format("EffP_%d",idex), 6);
    						EffP8 = readEffP(iBin, TString::Format("EffP_%d",idex), 7);
    						EffP9 = readEffP(iBin, TString::Format("EffP_%d",idex), 8);
    						EffP10 = readEffP(iBin, TString::Format("EffP_%d",idex),9);
    					//	Iafb = readParam(iBin, "Iafb", 0);
    					//	Ifh  = readParam(iBin, "Ifh", 0);
    					//	cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
    					
    					// Cconstrained!!!	
    					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);   // Fitted results!
    					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
    					//	iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 1);  // Selected fitted values!
    					//	ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 1);
    						iafb = readEffP(iBin, TString::Format("Iafb_%s_EffP_%d",outfile,idex), 0);  // Selected initial values!
    						ifh  = readEffP(iBin, TString::Format("Ifh_%s_EffP_%d",outfile,idex), 0);
    						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
    						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
    					// Unconstrained!!!	
    					//	Iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);  // Selected initial values!  Unconstrained
    					//	Ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
    					//	Iafb = tan( iafb * TMath::Pi() / tan( (ifh - 1.5) * TMath::Pi() / 3.) ) ; // Selected initial values!  Unconstrained
    					//	Ifh  = tan( (ifh - 1.5) * TMath::Pi() / 3. );
    					//	Iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 1);  // Selected fitted values!   Unconstrained
    					//	Ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 1);
    					//	vbin = angular2D_Limited_bin(iBin, Iafb, Ifh, Index);
    					//	vbin = angular2D_Limited_bin(iBin, Iafb, Ifh, Index, EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10 );
    						cout<<"1: "<<EffP1<<"; 2: "<<EffP2<<"; 3: "<<EffP3<<"; 4: "<<EffP4<<"; 5: "<<EffP5<<"; 6: "<<EffP6<<"; 7: "<<EffP7<<"; 8: "<<EffP8<<"; 9: "<<EffP9<<"; 10: "<<EffP10<<endl;
    						vbin = angular2D_Limited_bin(iBin, Iafb, Ifh, Index, EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10, idex );
    					}
					} else if (func == "angular_gen") {
						const char outfile[] = "angular_gen";
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular_gen_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_gen_R") {
						const char outfile[] = "angular_gen_R";
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular_gen_R_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco_R") {
						const char outfile[] = "angular_reco_R";
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
					//	Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
					//	Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						cout<<Iafb<<" "<<Ifh<<endl;
						vbin = angular_reco_R_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco_D") {
						const char outfile[] = "angular_reco_D";
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
					//	Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
					//	Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular_reco_D_bin(iBin, Iafb, Ifh, Index);
					} else if (func == "angular_reco") {
						const char outfile[] = "angular_reco";
					//	Iafb = readParam(iBin, TString::Format("afb_%s",outfile), 0);
					//	Ifh  = readParam(iBin, TString::Format("fh_%s",outfile), 0);
					//	Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
					//	Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
					}
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb_"<<func<<" = "<<Iafb<<" >>>>>> Ifh_"<<func<<" = "<<Ifh<<endl;
				} else {
					cout<<"Refit with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> refit "<<endl;
					return 0;
				}
			}
		}else if (iBin == 997) {
			if (func == "angular2D") {
				const char outfile[] = "angular2D";
				angular2D_LHCb(outfile);
            }
		}else if (iBin == 999) {
			if (func == "angular2D") {
				const char outfile[] = "angular2D";
				angular2D(outfile);
			} else if (func == "angular_gen") {
				const char outfile[] = "angular_gen";
				angular_gen(outfile);
			} else if (func == "angular_gen_R") {
				const char outfile[] = "angular_gen_R";
				angular_gen_R(outfile);
			} else if (func == "angular_reco_R") {
				const char outfile[] = "angular_reco_R";
				angular_reco_R(outfile);
			} else if (func == "angular_reco_D") {
				const char outfile[] = "angular_reco_D";
				angular_reco_D(outfile);
			} else if (func == "angular_reco") {
				const char outfile[] = "angular_reco";
				angular_reco(outfile);
			}
		}else { 
			cout<<"Refit data, iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
	} else if (func == "PlotFCN_Limited") {
		TString label = argv[4];
		if (label == "angular2D_Limited") {
			const char outfile[] = "angular2D_Limited";
			PlotFCN_Limited(iBin, outfile);
		}
	} else if (func == "PlotFCN") {
		TString label = argv[4];
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_gen") {
			const char outfile[] = "angular_gen";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_gen_R") {
			const char outfile[] = "angular_gen_R";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_reco_R") {
			const char outfile[] = "angular_reco_R";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_reco_D") {
			const char outfile[] = "angular_reco_D";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			PlotFCN(iBin, outfile);
		}
	//	const char outfile[]="FCN";
	} else if (func == "Plot2DFCN") {
		TString label = argv[4];
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			Plot2DFCN(iBin, outfile);
		} else if (label == "angular_gen") {
			const char outfile[] = "angular_gen";
			Plot2DFCN(iBin, outfile);
		} else if (label == "angular_gen_R") {
			const char outfile[] = "angular_gen_R";
			Plot2DFCN(iBin, outfile);
		} else if (label == "angular_reco_R") {
			const char outfile[] = "angular_reco_R";
			Plot2DFCN(iBin, outfile);
		} else if (label == "angular_reco_D") {
			const char outfile[] = "angular_reco_D";
			Plot2DFCN(iBin, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			Plot2DFCN(iBin, outfile);
		}
	//	const char outfile[]="FCN";
	} else if (func == "PlotAfbFh_i") {	
		TString label = argv[4];
		int Index;
		if (iBin >= 0 && iBin < 11) {
			Index = -999;
		}else if (iBin == 999) {
			Index = -1;
		}else { 
			cout<<"Plot Afb and Fh 2D distribution: for each q2 iBin counts from 0 to 10; or 999 to plot the refitted results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			PlotAfbFh_i(iBin, Index, outfile);
		} else if (label == "angular_gen") {
			const char outfile[] = "angular_gen";
			PlotAfbFh_i(iBin, Index, outfile);
		} else if (label == "angular_gen_R") {
			const char outfile[] = "angular_gen_R";
			PlotAfbFh_i(iBin, Index, outfile);
		} else if (label == "angular_reco_R") {
			const char outfile[] = "angular_reco_R";
			PlotAfbFh_i(iBin, Index, outfile);
		} else if (label == "angular_reco_D") {
			const char outfile[] = "angular_reco_D";
			PlotAfbFh_i(iBin, Index, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			PlotAfbFh_i(iBin, Index, outfile);
		}
	} else if (func == "PlotAfbFh_f") {	
		TString label = argv[4];
		int Index;
		if (iBin >= 0 && iBin < 11) {
			Index = -999;
		}else if (iBin == 999) {
			Index = -1;
		}else if (iBin == 998) {
			Index = -2;
		}else { 
			cout<<"Plot Afb and Fh 2D distribution: for each q2 iBin counts from 0 to 10; or 999 to plot the refitted results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_gen") {
			const char outfile[] = "angular_gen";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_gen_R") {
			const char outfile[] = "angular_gen_R";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_reco_R") {
			const char outfile[] = "angular_reco_R";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_reco_D") {
			const char outfile[] = "angular_reco_D";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			PlotAfbFh_f(iBin, Index, outfile);
		}
	} else if (func == "angular2D_Limited_FCN") {	
		TString label = argv[4];
		int Index;
		if (iBin >= 0 && iBin < 11) {
			Index = -999;
		}else if (iBin == 999) {
			Index = -1;
		}else if (iBin == 998) {
			Index = -2;
		}else { 
			cout<<"Plot Afb and Fh 2D distribution: for each q2 iBin counts from 0 to 10; or 999 to plot the refitted results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
		const char outfile[] = "angular2D_Limited";
		angular2D_Limited_FCN(iBin, Index, outfile);
    } else if (func == "angular2D_Profiled_Afb" || func == "angular2D_Profiled_Fh" || func == "angular2D_NLL"){
  		if (iBin >= 0 && iBin < 11) {
      		ch->Add(infile.Data());
      		if (ch == NULL) gSystem->Exit(0);
      		std::vector<double> vbin;
      			
      		double Iafb, Ifh, iafb, ifh;
  			int Index = -999;
  			if (argc == 7 || argc == 8) {   //////////////////////////////   scanning   ///////////////////////////////////////////
    			Index  = atoi(argv[4]);
    			Iafb   = atof(argv[5]);
    			Ifh    = atof(argv[6]);
                int idex = -1;
				cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
				if (func == "angular2D_NLL") {
				  vbin = angular2D_NLL_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular2D_Profiled_Afb") {
					idex = atof(argv[7]);
				  vbin = angular2D_Profiled_Afb_bin(iBin, Iafb, Ifh, Index, idex);
				} else if (func == "angular2D_Profiled_Fh") {
					idex = atof(argv[7]);
				  vbin = angular2D_Profiled_Fh_bin(iBin, Iafb, Ifh, Index, idex);
				} else {
				//	cout<<"A test of data fitting with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
					return 0;
				}
  			} else if (argc == 5 || argc == 6) {  //////////////////////////////   refit ////////////////////////////////////
    			TString type  = argv[4];
    			if ( type == "refit") {
      				Index = -1;
    				if (func == "angular2D_NLL") {
						const char outfile[] = "angular2D_NLL";
						iafb = readEffP(iBin, TString::Format("Iafb_%s_Pro",outfile), 0);  // Selected initial values!
						ifh  = readEffP(iBin, TString::Format("Ifh_%s_Fix",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular2D_NLL_bin(iBin, Iafb, Ifh, Index);
                        //vbin = angular2D_bin(iBin, iafb, ifh, Index);   // unconstrained fitting
  					} else if (func == "angular2D_Profiled_Afb") {
                        int idex = -1;
        				idex = atof(argv[5]);
						const char outfile[] = "angular2D_Profiled_Afb";
						iafb = readEffP(iBin, TString::Format("Iafb_%s_Pro_%d",outfile,idex), 0);  // Selected initial values!
						ifh  = readEffP(iBin, TString::Format("Ifh_%s_Fix_%d",outfile,idex), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular2D_Profiled_Afb_bin(iBin, Iafb, Ifh, Index, idex);
                        //vbin = angular2D_bin(iBin, iafb, ifh, Index);   // unconstrained fitting
      				} else if (func == "angular2D_Profiled_Fh") {
                        int idex = -1;
        				idex = atof(argv[5]);
        				const char outfile[] = "angular2D_Profiled_Fh";
        				iafb = readEffP(iBin, TString::Format("Iafb_%s_Fix_%d",outfile,idex), 0);  // Selected initial values!
						ifh  = readEffP(iBin, TString::Format("Ifh_%s_Pro_%d",outfile,idex), 0);
        				Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular2D_Profiled_Fh_bin(iBin, Iafb, Ifh, Index, idex);
                        //vbin = angular2D_bin(iBin, iafb, ifh, Index);   // unconstrained fitting
      				} else {
      				    return 0;
    				}
    				cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb_"<<func<<" = "<<Iafb<<" >>>>>> Ifh_"<<func<<" = "<<Ifh<<endl;
  				} else {
  					cout<<"Refit with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> refit "<<endl;
  					return 0;
  				}
  			}
		} else { 
  			cout<<"Refit data, iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
  			cout<<" Please check the Usage : ./fit"<<endl;
  			return 0; 
		}
/////////////////////////////////////////////////////////////////////////////////////////////
  	} else if (func == "PlotFCN_Profiled") {
    		TString label = argv[4];
    		if (label == "angular2D_Profiled_Afb") {
      			const char outfile[] = "angular2D_Profiled_Afb";
      			PlotFCN_Profiled_Afb(iBin, outfile);
    		} else if (label == "angular2D_Profiled_Fh") {
      			const char outfile[] = "angular2D_Profiled_Fh";
      			PlotFCN_Profiled_Fh(iBin, outfile);
        }
/////////////////////////////////////////////////////////////////////////////////////////////
  	} else if (func == "angular2D_Toy") {
  	   	ch->Add(infile.Data());
        summarypath=GetDirectory(argv[2]);
        cout<<summarypath<<endl;
  		  if (ch == NULL) gSystem->Exit(0);
  		  angular2D_Toy_bin(iBin); 
  	} else if (func == "angular2D_Toy_unCons") {
  		  ch->Add(infile.Data());
        summarypath=GetDirectory(argv[2]);
        cout<<summarypath<<endl;
  		  if (ch == NULL) gSystem->Exit(0);
  		  angular2D_Toy_unCons_bin(iBin); 
/////////////////////////////////////////////////////////////////////////////////////////////
    } else if (func == "createToys"){
        scaleFactor = 100;
        double expNCombBkg[11] = { 518.5 , 1125.1 , 1586.7 , 0 , 755.8 , 0 , 248.2 , 242.8 , 286.1 , 2355.5 , 4838.5 };
        double expNSig[11]     = { 173.5 ,  330.9 ,  783.3 , 0 , 360.2 , 0 , 214.8 , 262.1 , 225.9 ,  779.5 , 2275.4 };
        for (int iBin = 0; iBin < 11; iBin++) {
              genToySignal(iBin,expNSig[iBin]*scaleFactor);// (iBin, nEvts)
              genToyCombBkg(iBin,expNCombBkg[iBin]*scaleFactor,"combBkgToy");
        }
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        splitMCSamples();
        //rndPickMCSamples(scaleFactor);// (nSets,lumis)
    } else if (func == "ToyVali"){
        // Get parameters for the q2 Bin
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        if (!wspace) return 1;
        double  nsig = wspace->var("nsig")->getVal();
        double  nbkgComb = wspace->var("nbkgComb")->getVal();
        owspacepath=TString::Format("./NLL/bin%d",iBin);
        //genToySignal(iBin, nsig*100, 0.10,0.40); // v0
        //genToySignal(iBin, nsig*100, 0.01,0.05); // v1
        genToySignal(iBin, nsig*100, -0.30,0.80); // v2
        genToyCombBkg(iBin,nbkgComb*100);
/////////////////////////////////////////////////////////////////////////////////////////////
    } else if (func == "createFCToys"){
        createFCToys(iBin,nToy);
    } else if (func == "FCScan"){
        owspacepath=TString::Format("./ToyMC");
        int whichBin=-1;
        whichBin=iBin;
        harvestFCFitResults(whichBin,nToy);
        cout<<"----------------"<<endl;
//        getFCInterval(whichBin,nToy);
        getFCInterval2(whichBin,nToy);
/////////////////////////////////////////////////////////////////////////////////////////////
////    } else if (func == "createNLLToys"){
////        createNLLToys(iBin,nToy);
////        NLLwFix(iBin,nToy); 
////    } else if (func == "NLLScan"){
////        owspacepath=TString::Format("./NLLMC");
////        int whichBin=-1;
////        whichBin=iBin;
////        harvestNLLFitResults(whichBin,nToy);
////        cout<<"----------------"<<endl;
////        getNLLResults(whichBin,nToy);
/////////////////////////////////////////////////////////////////////////////////////////////
	} else if (func == "Average") {
		if (iBin != 999) return 0;
		const char outfile[] = "Average";
		Average(outfile);
	} else if (func == "BranchFraction_new") {
		if (iBin != 999) return 0;
		const char outfile[] = "BranchFraction_new";
		BranchFraction_new(iBin, outfile);
	} else if (func == "BranchFraction") {
		if (iBin != 999) return 0;
		const char outfile[] = "BranchFraction";
		BranchFraction(iBin, outfile);
	} else if (func == "read_pdf") {
		read_pdf();
/////////////////////////////////////////////////////////////////////////////////////////////
    } else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
