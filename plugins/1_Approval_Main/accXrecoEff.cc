// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------

void OPTEff(int iBin, const char outfile[] = "OPTEff")
{//{{{
	setTDRStyle();
	printf("Evaluate OPT efficiency for bin#%d\n",iBin);
	TLorentzVector *Tk_4vec;
    Tk_4vec = 0;
	double Bmass = 0;
	double Bvtxcl = 0;
	double Bcosalphabs2D = 0;
	double Blxysig = 0;
	double Trkdcasigbs = 0;
	double Q2 = 0;
	double CosThetaL = 0;
//	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Bvtxcl"        , 1);
	ch->SetBranchStatus("Blxysig"       , 1);
	ch->SetBranchStatus("Bcosalphabs2D" , 1);
	ch->SetBranchStatus("Trkdcasigbs"   , 1);
	ch->SetBranchStatus("Tk_4vec"       , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
    ch->SetBranchAddress("Bmass"        , &Bmass);
	ch->SetBranchAddress("Bvtxcl"           , &Bvtxcl);
	ch->SetBranchAddress("Blxysig"           , &Blxysig);
	ch->SetBranchAddress("Bcosalphabs2D"           , &Bcosalphabs2D);
	ch->SetBranchAddress("Trkdcasigbs"           , &Trkdcasigbs);
	ch->SetBranchAddress("Tk_4vec"           , &Tk_4vec);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	
//	Fill histograms
    TH1D *h1d_1_Bvtxcl = new TH1D("h1d_1_Bvtxcl","Bvtxcl",20,0.0,0.2);
    TH1D *h1d_1_Blxysig = new TH1D("h1d_1_Blxysig","Blxysig",50,3.0,18.);
    TH1D *h1d_1_Trkpt = new TH1D("h1d_1_Trkpt","Trkpt",20,0.6,2.6);
    TH1D *h1d_1_Trkdcabs = new TH1D("h1d_1_Trkdcabs","Trkdcasigbs",35,1.0,4.5);
    TH1D *h1d_1_Bcos2D = new TH1D("h1d_1_Bcos2D","Bcosalphabs2D",10,0.9990,1.0);
    TH1D *h1d_2_Bvtxcl = new TH1D("h1d_2_Bvtxcl","Bvtxcl",20,0.0,0.2);
    TH1D *h1d_2_Blxysig = new TH1D("h1d_2_Blxysig","Blxysig",50,3.0,18.);
    TH1D *h1d_2_Trkpt = new TH1D("h1d_2_Trkpt","Trkpt",20,0.6,2.6);
    TH1D *h1d_2_Trkdcabs = new TH1D("h1d_2_Trkdcabs","Trkdcasigbs",35,1.0,4.5);
    TH1D *h1d_2_Bcos2D = new TH1D("h1d_2_Bcos2D","Bcosalphabs2D",10,0.9990,1.0);
   
    TH2D *h2d_2_Bvtxcl = new TH2D("h2d_2_Bvtxcl","Bmass vs Bvtxcl",20,5.10,5.60,20,0.0,0.2);
    TH2D *h2d_2_Blxysig = new TH2D("h2d_2_Blxysig","Bmass vs Blxysig",20,5.10,5.60,50,3.0,18.);
    TH2D *h2d_2_Trkpt = new TH2D("h2d_2_Trkpt","Bmass vs Trkpt",20,5.10,5.60,20,0.6,2.6);
    TH2D *h2d_2_Trkdcabs = new TH2D("h2d_2_Trkdcabs","Bmass vs Trkdcasigbs",20,5.10,5.60,35,1.0,4.5);
    TH2D *h2d_2_Bcos2D = new TH2D("h2d_2_Bcos2D","Bmass vs Bcosalphabs2D",20,5.10,5.60,10,0.9990,1.0);
	Long64_t nentries1 = ch->GetEntriesFast();
    int N_Bvtxcl[20], N_Blxysig[20], N_Trkpt[20], N_Trkdcabs[20], N_Bcos2D[20];
    int Nbkg;
    double Eff_Bvtxcl[20], Eff_Blxysig[20], Eff_Trkpt[20], Eff_Trkdcabs[20], Eff_Bcos2D[20];
    double Eff_Bvtxcl_err[20], Eff_Blxysig_err[20], Eff_Trkpt_err[20], Eff_Trkdcabs_err[20], Eff_Bcos2D_err[20];
    double x_Bvtxcl[20], x_Blxysig[20], x_Trkpt[20], x_Trkdcabs[20], x_Bcos2D[20];
    double x_err[20];
    for (int i=0; i<20; i++) {
        x_err[i]=0.;
        N_Bvtxcl[i]=0;
        N_Blxysig[i]=0;
        N_Trkpt[i]=0;
        N_Trkdcabs[i]=0;
        N_Bcos2D[i]=0;
        Eff_Bvtxcl[i]=0.;
        Eff_Blxysig[i]=0.;
        Eff_Trkpt[i]=0.;
        Eff_Trkdcabs[i]=0.;
        Eff_Bcos2D[i]=0.;
        x_Bvtxcl[i]=0.01 + 0.01*i;
        x_Blxysig[i]=3.0 + 0.75*i;
        x_Trkpt[i]=0.6 + 0.1*i;
        x_Trkdcabs[i]=1.0 + 0.2*i;
        x_Bcos2D[i]=0.9980 + 0.0001*i;
    }
   
    for (Long64_t entry = 0; entry < ch->GetEntries(); entry++) {
    //for (Long64_t entry = 0; entry < nentries1; entry++) {
        ch->GetEntry(entry);
		if (Bmass != -999 ) { // HLT 2016-04-19
            Nbkg++;
            h1d_1_Bvtxcl->Fill(Bvtxcl);
            h1d_1_Blxysig->Fill(Blxysig);
            h1d_1_Trkpt->Fill(Tk_4vec->Pt());
            h1d_1_Trkdcabs->Fill(Trkdcasigbs);
            h1d_1_Bcos2D->Fill(Bcosalphabs2D);
            if (Bcosalphabs2D > 0.9997 && Tk_4vec->Pt() > 1.3 && Trkdcasigbs > 3.3 && Blxysig > 10.6) {
                h2d_2_Bvtxcl->Fill(Bmass,Bvtxcl);
                h1d_2_Bvtxcl->Fill(Bvtxcl);
                for (int j=0; j<20; j++) {
                if (Bvtxcl > x_Bvtxcl[j]) N_Bvtxcl[j]++;
                }
            }  
            if (Bcosalphabs2D > 0.9997 && Tk_4vec->Pt() > 1.3 && Trkdcasigbs > 3.3 && Bvtxcl > 0.12) {
                h2d_2_Blxysig->Fill(Bmass,Blxysig);            
                h1d_2_Blxysig->Fill(Blxysig);
                for (int j=0; j<20; j++) {
                    if (Blxysig > x_Blxysig[j]) N_Blxysig[j]++;
                }
            }
            if (Bcosalphabs2D > 0.9997 && Tk_4vec->Pt() > 1.3 && Bvtxcl > 0.12 && Blxysig > 10.6) {
                h2d_2_Trkdcabs->Fill(Bmass,Trkdcasigbs);
                h1d_2_Trkdcabs->Fill(Trkdcasigbs);
                for (int j=0; j<20; j++) {
                    if (Trkdcasigbs > x_Trkdcabs[j]) N_Trkdcabs[j]++;
                }
            }
            if (Bcosalphabs2D > 0.9997 && Trkdcasigbs > 3.3 && Bvtxcl > 0.12 && Blxysig > 10.6) {
                h2d_2_Trkpt->Fill(Bmass,Tk_4vec->Pt());
                h1d_2_Trkpt->Fill(Tk_4vec->Pt());
                for (int j=0; j<20; j++) {
                    if (Tk_4vec->Pt() > x_Trkpt[j]) N_Trkpt[j]++;
                }
            }
            if (Tk_4vec->Pt() > 1.3 && Trkdcasigbs > 3.3 && Bvtxcl > 0.12 && Blxysig > 10.6) {
                h2d_2_Bcos2D->Fill(Bmass,Bcosalphabs2D);
                h1d_2_Bcos2D->Fill(Bcosalphabs2D);
                for (int j=0; j<20; j++) {
                    if (Bcosalphabs2D > x_Bcos2D[j]) N_Bcos2D[j]++;
                }
            }
        }
    }
    for (int i=0; i<20; i++) {
        Eff_Bvtxcl[i]=N_Bvtxcl[i]*1.0 / Nbkg;
        Eff_Blxysig[i]=N_Blxysig[i]*1.0 / Nbkg;
        Eff_Trkdcabs[i]=N_Trkdcabs[i]*1.0 / Nbkg;
        Eff_Trkpt[i]=N_Trkpt[i]*1.0 / Nbkg;
        Eff_Bcos2D[i]=N_Bcos2D[i]*1.0 / Nbkg;
		
        Eff_Bvtxcl_err[i]=sqrt(Eff_Bvtxcl[i]*(1-Eff_Bvtxcl[i])/N_Bvtxcl[i]);
        Eff_Blxysig_err[i]=sqrt(Eff_Blxysig[i]*(1-Eff_Blxysig[i])/N_Blxysig[i]);
        Eff_Trkdcabs_err[i]=sqrt(Eff_Trkdcabs[i]*(1-Eff_Trkdcabs[i])/N_Trkdcabs[i]);
        Eff_Trkpt_err[i]=sqrt(Eff_Trkpt[i]*(1-Eff_Trkpt[i])/N_Trkpt[i]);
        Eff_Bcos2D_err[i]=sqrt(Eff_Bcos2D[i]*(1-Eff_Bcos2D[i])/N_Bcos2D[i]);
    }
    TGraphErrors *h1d_3_Bvtxcl = new TGraphErrors(20,x_Bvtxcl,Eff_Bvtxcl,x_err,Eff_Bvtxcl_err);
    TGraphErrors *h1d_3_Blxysig = new TGraphErrors(20,x_Blxysig,Eff_Blxysig,x_err,Eff_Blxysig_err);
    TGraphErrors *h1d_3_Trkpt = new TGraphErrors(20,x_Trkpt,Eff_Trkpt,x_err,Eff_Trkpt_err);
    TGraphErrors *h1d_3_Trkdcabs = new TGraphErrors(20,x_Trkdcabs,Eff_Trkdcabs,x_err,Eff_Trkdcabs_err);
    TGraphErrors *h1d_3_Bcos2D = new TGraphErrors(20,x_Bcos2D,Eff_Bcos2D,x_err,Eff_Bcos2D_err);
   
    h1d_3_Bvtxcl->GetYaxis()->SetRangeUser(0.0,0.35);  // 0.32
    h1d_3_Blxysig->GetYaxis()->SetRangeUser(0.0,0.35);
    h1d_3_Trkpt->GetYaxis()->SetRangeUser(0.0,0.35);
    h1d_3_Trkdcabs->GetYaxis()->SetRangeUser(0.0,0.35);
    h1d_3_Bcos2D->GetYaxis()->SetRangeUser(0.0,0.35);
   
    h1d_3_Bvtxcl->GetXaxis()->SetTitle("IF Bvtxcl > ");
    h1d_3_Bvtxcl->GetYaxis()->SetTitle("Efficiency");
    h1d_3_Blxysig->GetXaxis()->SetTitle("IF Blxysig >");
    h1d_3_Blxysig->GetYaxis()->SetTitle("Efficiency");
    h1d_3_Trkpt->GetXaxis()->SetTitle("IF Trkpt >");
    h1d_3_Trkpt->GetYaxis()->SetTitle("Efficiency");
    h1d_3_Trkdcabs->GetXaxis()->SetTitle("IF Trkdcasigbs >");
    h1d_3_Trkdcabs->GetYaxis()->SetTitle("Efficiency");
    h1d_3_Bcos2D->GetXaxis()->SetTitle("IF Bcosalphabs2D >");
    h1d_3_Bcos2D->GetYaxis()->SetTitle("Efficiency");

//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
/*
    h2d_2_Bvtxcl->Draw("COLZ");
	canvas.Update();
	canvas.Print("./plots/OPTEff_Bvtxcl.pdf");
    h2d_2_Blxysig->Draw("COLZ");
	canvas.Update();
	canvas.Print("./plots/OPTEff_Blxysig.pdf");
    h2d_2_Trkpt->Draw("COLZ");
	canvas.Update();
	canvas.Print("./plots/OPTEff_Trkpt.pdf");
    h2d_2_Trkdcabs->Draw("COLZ");
	canvas.Update();
	canvas.Print("./plots/OPTEff_Trkdcabs.pdf");
    h2d_2_Bcos2D->Draw("COLZ");
	canvas.Update();
	canvas.Print("./plots/OPTEff_Bcos2D.pdf");
*/
    h1d_3_Bvtxcl->Draw("APE");
	canvas.Update();
	canvas.Print("./plots/OPT3_Bvtxcl.pdf");
    h1d_3_Blxysig->Draw("APE");
	canvas.Update();
	canvas.Print("./plots/OPT3_Blxysig.pdf");
    h1d_3_Trkpt->Draw("APE");
	canvas.Update();
	canvas.Print("./plots/OPT3_Trkpt.pdf");
    h1d_3_Trkdcabs->Draw("APE");
	canvas.Update();
	canvas.Print("./plots/OPT3_Trkdcabs.pdf");
    h1d_3_Bcos2D->Draw("APE");
	canvas.Update();
	canvas.Print("./plots/OPT3_Bcos2D.pdf");
/*
    h1d_1_Bvtxcl->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT1_Bvtxcl.pdf");
    h1d_1_Blxysig->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT1_Blxysig.pdf");
    h1d_1_Trkpt->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT1_Trkpt.pdf");
    h1d_1_Trkdcabs->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT1_Trkdcabs.pdf");
    h1d_1_Bcos2D->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT1_Bcos2D.pdf");
	
    h1d_2_Bvtxcl->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT2_Bvtxcl.pdf");
    h1d_2_Blxysig->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT2_Blxysig.pdf");
    h1d_2_Trkpt->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT2_Trkpt.pdf");
    h1d_2_Trkdcabs->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT2_Trkdcabs.pdf");
    h1d_2_Bcos2D->Draw("PE");
	canvas.Update();
	canvas.Print("./plots/OPT2_Bcos2D.pdf");
*/   
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
void createAccptanceHist() // create acceptance histogram from UNFILTERED GEN.
{//{{{
	setTDRStyle();
//	double ngen[11]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	double nreco[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double accUpperBound = 0.050;
	double gQ2 = 0;
	double gCosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
    int    count_G[11];
    int    count_g[11];
    for (int i = 0; i < 11; i++) { count_G[i] = 0; count_g[i] = 0;}
	TChain *treein=new TChain("tree");
	treein->Add("../RootFiles/MC_GENOnly/BToKMuMu_GENOnly_8TeV_genonly_v5+7.root");
	if (treein == NULL) gSystem->Exit(0);
	treein->SetBranchStatus("*",0);
	treein->SetBranchStatus("genQ2"         , 1);
	treein->SetBranchStatus("genCosThetaL"  , 1);
	treein->SetBranchStatus("genMu*"        , 1);
	treein->SetBranchAddress("genQ2"        , &gQ2);
	treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	treein->SetBranchAddress("genMupPt"     , &gmuppt);
	treein->SetBranchAddress("genMupEta"    , &gmupeta);
	treein->SetBranchAddress("genMumPt"     , &gmumpt);
	treein->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Create histograms
	TFile *fout = new TFile("./RootFiles/acceptance_8TeV.root","RECREATE");
//	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
//	TH1F *h1_ngen[11];
//	TH1F *h1_nacc[11];
//	TH1F *h1_acc[11];
	TH1F *h1_ngen_fine[11];
	TH1F *h1_nacc_fine[11];
	TH1F *h1_acc_fine[11];
	for(int iBin = 0; iBin < 11; iBin++){
	//	if (iBin == 3 || iBin == 5) continue;    ///////////////////////////////////////////////////////////// 12-23
//		h1_ngen[iBin] = new TH1F(TString::Format("h1_ngen_bin%d",iBin),"h1_ngen",6,thetaLBins);
//		h1_nacc[iBin] = new TH1F(TString::Format("h1_nacc_bin%d",iBin) ,"h1_nacc" ,6,thetaLBins); 
//		h1_acc [iBin] = new TH1F(TString::Format("h1_acc_bin%d",iBin),"",6,thetaLBins);
		h1_ngen_fine[iBin] = new TH1F(TString::Format("h1_ngen_fine_bin%d",iBin),"h1_ngen",20,-1,1);
		h1_nacc_fine[iBin] = new TH1F(TString::Format("h1_nacc_fine_bin%d",iBin) ,"h1_nacc" ,20,-1,1); 
		h1_acc_fine[iBin]  = new TH1F(TString::Format("h1_acc_fine_bin%d",iBin),"",20,-1,1);
//		h1_ngen[iBin]->SetTitleOffset(1.3,"XY");
//		h1_ngen[iBin]->SetXTitle("genCosThetaL");
//		h1_ngen[iBin]->SetYTitle("Generated events");
//		h1_nacc[iBin]->SetTitleOffset(1.3,"XY");
//		h1_nacc[iBin]->SetXTitle("genCosThetaL");
//		h1_nacc[iBin]->SetYTitle("Events in acceptance");
//		h1_acc [iBin]->SetStats(0);
//		h1_acc [iBin]->SetMinimum(0.);
//		h1_acc [iBin]->SetMaximum(accUpperBound);
//		h1_acc [iBin]->SetTitleOffset(1.3,"XY");
//		h1_acc [iBin]->SetXTitle("genCosThetaL");
//		h1_acc [iBin]->SetYTitle("Acceptance");
		
	//	h1_ngen_fine[iBin]->SetTitleOffset(1.3,"XY");
		h1_ngen_fine[iBin]->SetXTitle("cos#theta_{l}^{gen}");
		h1_ngen_fine[iBin]->SetYTitle("Generated events");
	//	h1_nacc_fine[iBin]->SetTitleOffset(1.3,"XY");
		h1_nacc_fine[iBin]->SetXTitle("cos#theta_{l}^{gen}");
		h1_nacc_fine[iBin]->SetYTitle("Events in acceptance");
		h1_acc_fine [iBin]->SetStats(0);
		h1_acc_fine [iBin]->SetMinimum(0.);
		h1_acc_fine [iBin]->SetMaximum(accUpperBound);
	//	h1_acc_fine [iBin]->SetTitleOffset(1.3,"XY");
		h1_acc_fine [iBin]->SetXTitle("cos#theta_{l}^{gen}");
		h1_acc_fine [iBin]->SetYTitle("Acceptance");
	}

//	Fill histograms
	// Read data
	for (int entry = 0; entry < treein->GetEntries(); entry++) {
		treein->GetEntry(entry);
		for(int iBin = 0; iBin < 11; iBin++){
		//	if (iBin == 3 || iBin == 5) continue;    //////////////////////////////////// 12-23
			if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
            count_g[iBin] = count_g[iBin] + 1;
//			h1_ngen[iBin]->Fill(gCosThetaL);
			h1_ngen_fine[iBin]->Fill(gCosThetaL);
//			ngen[iBin] = ngen[iBin] + 1;   
			if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
		//	if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 2.8 && gmuppt > 2.8 ){
      	        count_G[iBin] = count_G[iBin] + 1;
//				h1_nacc[iBin]->Fill(gCosThetaL);
				h1_nacc_fine[iBin]->Fill(gCosThetaL);
//				nreco[iBin] = nreco[iBin] + 1;   
			}
		}
	}
	for(int iBin = 0; iBin < 11; iBin++){
        cout<<"Acc_"<<iBin<<" = "<<count_G[iBin]<<" / "<<count_g[iBin]<<" = "<<1.0*count_G[iBin]/count_g[iBin]<<endl<<endl;
        const double aaa = h1_nacc_fine[iBin]->GetEntries();
        const double bbb = h1_ngen_fine[iBin]->GetEntries();
        cout<<"New_"<<iBin<<" = "<<aaa<<" / "<<bbb<<" = "<<aaa/bbb<<endl<<endl;
//	cout<<"Acceptance["<<iBin<<"] = "<<nreco[iBin]<<" / "<<ngen[iBin]<<" = "<<nreco[iBin]/ngen[iBin]<<endl;
//	if (iBin == 3 || iBin == 5) continue;  /////////////////////////////////////////////////////////////// 12-23
//	Calculate acceptance
//	h1_acc[iBin]->SetAxisRange(0.,1.,"Y");
//	for (int i = 1; i <= 6; i++) {
//	Fill acceptance
//	if (h1_ngen[iBin]->GetBinContent(i) == 0) {
//		printf("WARNING: Acceptance(%d)=%f/%f\n",i,h1_nacc[iBin]->GetBinContent(i),h1_ngen[iBin]->GetBinContent(i));
//		h1_acc[iBin]->SetBinContent(i,0.);
//		h1_acc[iBin]->SetBinError(i,1.);
//	}else{
//		h1_acc[iBin]->SetBinContent(i,h1_nacc[iBin]->GetBinContent(i)/h1_ngen[iBin]->GetBinContent(i));
//		if (h1_nacc[iBin]->GetBinContent(i) != 0){
//		    h1_acc[iBin]->SetBinError(i,sqrt(h1_acc[iBin]->GetBinContent(i)*(1.-h1_acc[iBin]->GetBinContent(i))/h1_ngen[iBin]->GetBinContent(i)));
//		}else{
//			h1_acc[iBin]->SetBinError(i,0.);
//		}
//	}
//}
//printf("INFO: h1_acc_bin%d built.\n",iBin);
//		
//	h1_acc_fine[iBin]->SetAxisRange(0.,1.,"Y");
    int ggg = 0, GGG = 0;
	for (int i = 1; i <= 20; i++) {//L
    //	Fill acceptance
		if (h1_nacc_fine[iBin]->GetBinContent(i) == 0 || h1_ngen_fine[iBin]->GetBinContent(i) == 0) {
			h1_acc_fine[iBin]->SetBinContent(i,0.);
			h1_acc_fine[iBin]->SetBinError(i,1.);
		}else{
			h1_acc_fine[iBin]->SetBinContent(i,h1_nacc_fine[iBin]->GetBinContent(i)/h1_ngen_fine[iBin]->GetBinContent(i));
			if (h1_nacc_fine[iBin]->GetBinContent(i) != 0){
				h1_acc_fine[iBin]->SetBinError(i,sqrt(h1_acc_fine[iBin]->GetBinContent(i)*(1.-h1_acc_fine[iBin]->GetBinContent(i))/h1_ngen_fine[iBin]->GetBinContent(i)));
			}else{
				h1_acc_fine[iBin]->SetBinError(i,1.);
			}
		}
        ggg = ggg + h1_nacc_fine[iBin]->GetBinContent(i);
        GGG = GGG + h1_ngen_fine[iBin]->GetBinContent(i);
//		if ( i == 20 && iBin == 0) { cout<<h1_acc_fine[iBin]->GetBinContent(i)<<endl<<h1_nacc_fine[iBin]->GetBinContent(i)<<endl<<h1_ngen_fine[iBin]->GetBinContent(i)<<endl; }
	}
        cout<<"GGG_"<<iBin<<" = "<<ggg<<" / "<<GGG<<" = "<<1*ggg/GGG<<endl<<endl;
	printf("INFO: h1_acc_fine_bin%d built.\n",iBin);
//	Draw FitResult
	TString f1_model_format_0 = "[0] + [1]*exp(-0.5* ( ((x-[2])/[3])**2)) ";
	const int nPar = 4;
	TF1 *f1_model = new TF1 ("f1_model", f1_model_format_0, -1., 1.);
	if (iBin==0 || iBin==1 || iBin==9) f1_model->FixParameter(0,0.);
//	f1_model->SetParameter(1,0.01);
//	f1_model->SetParameter(2,0.);
	f1_model->SetParameter(3,10);
		
	TCanvas canvas("canvas");
	h1_acc_fine[iBin]->SetMinimum(0.);
	h1_acc_fine[iBin]->SetXTitle("cos#it{#theta_{l}}^{gen}");
	h1_acc_fine[iBin]->SetYTitle("Acceptance / 0.1");
	h1_acc_fine[iBin]->SetStats(0);
	h1_acc_fine[iBin]->SetTitleOffset(1.15,"Y");
	//	h1_acc_fine[iBin]->SetMaximum(effUpperBound2[iBin]);
	//	if (iBin == 0) h1_acc_fine[iBin]->SetMaximum(0.01);
	h1_acc_fine[iBin]->Draw("PE1");
	//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	//	h1_acc_fine[iBin]->Draw();
	h1_acc_fine[iBin]->Fit(f1_model,"R"); //// 09-09
   
	Double_t matrix[4][4];
	gMinuit->mnemat(&matrix[0][0],4);
	cout<<"Accepatance ERROR MATRIX: "<<endl;
	cout<<"                 1             2            3            4       "<<endl;
	for (int i = 0; i<4; i++) {
		cout<<"AccErr"<<i+1<<"  ";
		for (int j =0; j<4; j++) {
			cout<<matrix[i][j]<<" ";
		}
		cout<<endl;
	}
	
	f1_model->SetTitle("");
//	f1_model->SetMaximum(effUpperBound[iBin]); //03-11
	f1_model->SetLineWidth(2);
//	f1_model->SetRange(-0.89,0.79);
//	if (iBin == 1) f1_model->SetRange(-0.89,0.89);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
//	Draw compare
	TLatex *latex = new TLatex();
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t1->DrawLatex(.38,.90,TString::Format("Preliminary"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	t1->DrawLatex(.22,.90,TString::Format("Simulation"));
	t2->DrawLatex(.18,.82,TString::Format("#it{q^{2}}:"));
	t2->DrawLatex(.22,.82,Q2String[iBin]);
    if (iBin == 10) {
        t2->DrawLatex(.22,.77, " 14.18-22.00) GeV^{2}");
    }
	t2->DrawLatex(.84,.90,TString::Format("8TeV"));
	
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	double errPar1[nPar], errPar2[nPar], errPar3[nPar], errPar4[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
		errPar1[iPar]   = matrix[0][iPar];
		errPar2[iPar]   = matrix[1][iPar];
		errPar3[iPar]   = matrix[2][iPar];
		errPar4[iPar]   = matrix[3][iPar];
	}
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"acc",   arrPar,   4);  // f1_model_format_0
	writeParam(iBin,"accErr",arrParErr,4);  // f1_model_format_0
	writeParam(iBin,"AccErr1",       errPar1,  nPar);  
	writeParam(iBin,"AccErr2",       errPar2,  nPar);  
	writeParam(iBin,"AccErr3",       errPar3,  nPar);  
	writeParam(iBin,"AccErr4",       errPar4,  nPar);  
//	return output.c_str();	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_accL_fine_bin%d.pdf",iBin));
	}
//	cout<<" gen = "<<gen<<endl;   ///////////////////////////////////////////////////////////////////////////////// 12-23
//	cout<<"reco = "<<reco<<endl;
	fout->Write();
	fout->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void Accptance() // create acceptance histogram from UNFILTERED GEN.
{//{{{
	setTDRStyle();
	Long64_t ngen[6]  = {0, 0, 0, 0, 0, 0};
	Long64_t nacc[6]  = {0, 0, 0, 0, 0, 0};
	Long64_t nGEN[6]  = {0, 0, 0, 0, 0, 0};
	Long64_t nRECO[6] = {0, 0, 0, 0, 0, 0};
//	double ngen[11]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	double nacc[11]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	double nGEN[11]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	double nRECO[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	TLorentzVector  *Mum_4vec;
    TLorentzVector  *Mup_4vec;
    TBranch        *b_Mum_4vec;   //!
    TBranch        *b_Mup_4vec;   //!
    Mum_4vec = 0;
    Mup_4vec = 0;
    Double_t        Bmass;
    TBranch        *b_Bmass;   //!
	double gBEta  = 0;
	double gQ2    = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	double GQ2    = 0;
	double GBEta  = 0;
	double Gmuppt = 0;
	double Gmupeta= 0;
	double Gmumpt = 0;
	double Gmumeta= 0;
    int    count_G[11];
    int    count_g[11];
    for (int i = 0; i < 11; i++) { count_G[i] = 0; count_g[i] = 0;}
	TChain *treein=new TChain("tree");
//	treein->Add("../RootFiles/MC_GENOnly/BToJpsiK_MC_GENOnly_8TeV_genonly_v5.root");
	treein->Add("../RootFiles/MC_GENOnly/BToKMuMu_GENOnly_8TeV_genonly_v5+7.root");
	if (treein == NULL) gSystem->Exit(0);
	TChain *treeIN=new TChain("tree");
//	treeIN->Add("../RootFiles/Files/MC_JpsiK_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT+Anti3_GEN.root");
	treeIN->Add("../CutASelection/RootFiles/MC_Signal_8TeV_v4_AllCut.root");
	if (treeIN == NULL) gSystem->Exit(0);
	treein->SetBranchStatus("*",0);
	treein->SetBranchStatus("genQ2"         , 1);
	treein->SetBranchStatus("genBEta"       , 1);
	treein->SetBranchStatus("genMu*"        , 1);
	treein->SetBranchAddress("genQ2"        , &gQ2);
	treein->SetBranchAddress("genBEta"      , &gBEta);
	treein->SetBranchAddress("genMupPt"     , &gmuppt);
	treein->SetBranchAddress("genMupEta"    , &gmupeta);
	treein->SetBranchAddress("genMumPt"     , &gmumpt);
	treein->SetBranchAddress("genMumEta"    , &gmumeta);
	treeIN->SetBranchStatus("*",0);
	treeIN->SetBranchStatus("Mum_4vec"      , 1);
//	treeIN->SetBranchStatus("Mup_4vec"      , 1);
//	treeIN->SetBranchStatus("genQ2"         , 1);
	treeIN->SetBranchStatus("genMu*"        , 1);
    treeIN->SetBranchAddress("Mum_4vec", &Mum_4vec, &b_Mum_4vec);
    treeIN->SetBranchAddress("Mup_4vec", &Mup_4vec, &b_Mup_4vec);
    treeIN->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
//	treeIN->SetBranchAddress("Mum_4vec"     , &Mum_4vec);
//	treeIN->SetBranchAddress("Mup_4vec"     , &Mup_4vec);
	treeIN->SetBranchAddress("genQ2"        , &GQ2);
	treeIN->SetBranchAddress("genBEta"      , &GBEta);
	treeIN->SetBranchAddress("genMupPt"     , &Gmuppt);
	treeIN->SetBranchAddress("genMupEta"    , &Gmupeta);
	treeIN->SetBranchAddress("genMumPt"     , &Gmumpt);
	treeIN->SetBranchAddress("genMumEta"    , &Gmumeta);
	// Read data
	for (Long64_t entry = 0; entry < treein->GetEntries(); entry++) {
		treein->GetEntry(entry);
//	   // ---------------------------------------------
		for(int iBin = 0; iBin < 11; iBin++){
		//	if (iBin == 3 || iBin == 5) continue;    //////////////////////////////////// 12-23
			if (gQ2 > Q2rangeup[iBin] || gQ2 < Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
	//	    if (gQ2 > Q2rangeup[3] || gQ2 < Q2rangedn[3]) continue;
	//	    if (gQ2 > Q2rangeup[5] || gQ2 < Q2rangedn[5]) continue;
//		    if (gQ2 > Q2rangeup[10] || gQ2 <= Q2rangedn[10]) continue;
//		    if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
//		    if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
	//	    if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 ) {
             count_g[iBin] = count_g[iBin] + 1;
	//	    if ( fabs(gBEta) < 2.5 ) {
			ngen[5] = ngen[5] + 1;   
	//	    }
	//	    if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 2.8 && gmuppt > 2.8 ){
		    if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
	//	    if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 && fabs(gBEta) < 2.5 ){
      	    count_G[iBin] = count_G[iBin] + 1;
			nacc[5] = nacc[5] + 1;   
            }
	    }   // ---------------------------------------------
    }
	for(int iBin = 0; iBin < 11; iBin++){
        cout<<"Acc_"<<iBin<<" = "<<count_G[iBin]<<" / "<<count_g[iBin]<<" = "<<1.0*count_G[iBin]/count_g[iBin]<<endl<<endl;
    }
/*	cout<<"Acceptance0 = "<<nacc[0]<<" / "<<ngen[0]<<" = "<<1.*nacc[0]/ngen[0]<<endl;
	cout<<"Acceptance1 = "<<nacc[1]<<" / "<<ngen[1]<<" = "<<1.*nacc[1]/ngen[1]<<endl;
	cout<<"Acceptance2 = "<<nacc[2]<<" / "<<ngen[2]<<" = "<<1.*nacc[2]/ngen[2]<<endl;
	cout<<"Acceptance3 = "<<nacc[3]<<" / "<<ngen[3]<<" = "<<1.*nacc[3]/ngen[3]<<endl;
	cout<<"Acceptance4 = "<<nacc[4]<<" / "<<ngen[4]<<" = "<<1.*nacc[4]/ngen[4]<<endl;
*/	cout<<"Acceptance5 = "<<nacc[5]<<" / "<<ngen[5]<<" = "<<1.*nacc[5]/ngen[5]<<endl;
	
/*   
    for (Long64_t entry = 0; entry < treeIN->GetEntries(); entry++) {
		treeIN->GetEntry(entry);
		if ( Bmass != -999 ){
			if ( fabs(Mum_4vec->Eta()) < 0.6 ) {
				nRECO[0] = nRECO[0] + 1;   
			} else if ( fabs(Mum_4vec->Eta()) < 1.1  && fabs(Mum_4vec->Eta()) > 0.6 ) {
				nRECO[1] = nRECO[1] + 1;   
			} else if ( fabs(Mum_4vec->Eta()) < 1.45 && fabs(Mum_4vec->Eta()) > 1.1 ) {
				nRECO[2] = nRECO[2] + 1;   
			} else if ( fabs(Mum_4vec->Eta()) < 1.8  && fabs(Mum_4vec->Eta()) > 1.45) {
				nRECO[3] = nRECO[3] + 1;   
			} else if ( fabs(Mum_4vec->Eta()) < 2.4  && fabs(Mum_4vec->Eta()) > 1.8 ) {
				nRECO[4] = nRECO[4] + 1;   
			} 
		//	if ( fabs(Mum_4vec->Eta()) < 2.2 ) {
				nRECO[5] = nRECO[5] + 1;   
		//	}
	    }
		if ( fabs(Gmumeta) > 2.5 || fabs(Gmupeta) > 2.5 || Gmumpt < 3.5 || Gmuppt < 3.5 || fabs(GBEta) > 2.5 ) continue;
	//	if ( fabs(Gmumeta) > 2.4 || fabs(Gmupeta) > 2.4 || Gmumpt < 3.5 || Gmuppt < 3.5 || fabs(GBEta) > 2.4 ) continue;
	//	if ( fabs(Gmumeta) > 2.5 || fabs(Gmupeta) > 2.5 || Gmumpt < 3.5 || Gmuppt < 3.5 ) continue;
	//	if ( fabs(Gmumeta) > 2.2 || fabs(Gmupeta) > 2.2 || Gmumpt < 3.5 || Gmuppt < 3.5 ) continue;
	//	if (GQ2 > Q2rangeup[3] || GQ2 < Q2rangedn[3]) continue;
	//	if (GQ2 > Q2rangeup[5] || GQ2 < Q2rangedn[5]) continue;
//		if (GQ2 > Q2rangeup[10] || GQ2 <= Q2rangedn[10]) continue;
		if (GQ2 < Q2rangeup[3] && GQ2 > Q2rangedn[3]) continue;
		if (GQ2 < Q2rangeup[5] && GQ2 > Q2rangedn[5]) continue;
		if ( fabs(Gmumeta) < 0.6 ) {
			nGEN[0] = nGEN[0] + 1;   
		} else if ( fabs(Gmumeta) < 1.1  && fabs(Gmumeta) > 0.6 ) {
			nGEN[1] = nGEN[1] + 1;   
		} else if ( fabs(Gmumeta) < 1.45 && fabs(Gmumeta) > 1.1 ) {
			nGEN[2] = nGEN[2] + 1;   
		} else if ( fabs(Gmumeta) < 1.8  && fabs(Gmumeta) > 1.45) {		
			nGEN[3] = nGEN[3] + 1;   
		} else if ( fabs(Gmumeta) < 2.4  && fabs(Gmumeta) > 1.8 ) {
			nGEN[4] = nGEN[4] + 1;   
		} 
	//	if ( fabs(Gmumeta) < 2.5 && fabs(Gmupeta) < 2.5 && Gmumpt > 3.5 && Gmuppt > 3.5 ){
			nGEN[5] = nGEN[5] + 1;   
	//	}
	//	if ( fabs(Gmumeta) < 2.5 && fabs(Gmupeta) < 2.5 && Gmumpt > 3.5 && Gmuppt > 3.5 ){
	//	cout<<Bmass<<endl;
	}
	cout<<"RecoEfficiency0 = "<<nRECO[0]<<" / "<<nGEN[0]<<" = "<<1.0*nRECO[0]/nGEN[0]<<endl;
	cout<<"RecoEfficiency1 = "<<nRECO[1]<<" / "<<nGEN[1]<<" = "<<1.0*nRECO[1]/nGEN[1]<<endl;
	cout<<"RecoEfficiency2 = "<<nRECO[2]<<" / "<<nGEN[2]<<" = "<<1.0*nRECO[2]/nGEN[2]<<endl;
	cout<<"RecoEfficiency3 = "<<nRECO[3]<<" / "<<nGEN[3]<<" = "<<1.0*nRECO[3]/nGEN[3]<<endl;
	cout<<"RecoEfficiency4 = "<<nRECO[4]<<" / "<<nGEN[4]<<" = "<<1.0*nRECO[4]/nGEN[4]<<endl;
	cout<<"RecoEfficiency5 = "<<nRECO[5]<<" / "<<nGEN[5]<<" = "<<1.0*nRECO[5]/nGEN[5]<<endl;
*/	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void singleMuonEff(int iBin)
{//
	setTDRStyle();
	TFile f("../RootFiles/singleMuonEff/singleMuonEff_noTracking_L3ptg2_final.root", "read");
    TF1* user_f1 = (TF1*)f.Get(TString::Format("fitTotEff_DATA_pt_etaBin%d",iBin));
//	Draw FitResult
	TString f1_model_format_0 = "[0] + [1]*1.0 /(1.0 + exp( - ([2]*x+[3]) ) ) ";
	const int nPar = 4;
	//TF1 *f1_model = new TF1 ("f1_model", f1_model_format_0, 3.45, 30.45);
	TF1 *f1_model = new TF1 ("f1_model", f1_model_format_0, 3.3, 30.5);
	f1_model->SetParameter(0,0.1);
	f1_model->SetParameter(1,5.0);
	f1_model->SetParameter(2,0.1);
	f1_model->SetParameter(3,10);
   
    double xx[540], yy[540];
    for (int i=0; i<540; i++) {
        xx[i]=i/20.+3.45;
        yy[i]=user_f1->Eval(xx[i]);
    }
    TCanvas canvas("canvas");
	TH1F *frame = new TH1F("frame","",32,0.,32.);
	frame->SetStats(kFALSE);
	frame->SetTitle("");
	frame->Draw();
	frame->SetXTitle("p^{#mu}_{T} [GeV]");
	frame->SetYTitle("Single #mu efficiency / 0.05[GeV]");
	
    TGraph *h1=new TGraph(540,xx, yy);
	
   //h1->Draw("ACP");
    h1->Draw(" SAME P");
	h1->SetMarkerColor(4);
	h1->SetMarkerStyle(24);
	h1->SetMarkerSize(0.8);
    h1->Fit(f1_model,"R");
    h1->SetTitle("");
	h1->SetMaximum(1.05); 
	h1->SetMinimum(0.0); 
	h1->SetLineWidth(1);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
/*
	TPaveText* paveText = new TPaveText( 0.77, 0.15, 0.87, 0.25, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(0);
    paveText->AddText(Form("#eta bin %d ", iBin));
    paveText->Draw();
*/   
    TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);	
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	//t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
    float Etarangedn[10] = {0.0, 0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8};
    float Etarangeup[10] = {0.2, 0.3, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
    t1->DrawLatex(.62,.25,TString::Format(" %.1f < #left|#eta(#mu)#right| < %.1f", Etarangedn[iBin], Etarangeup[iBin]));
	
    double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
	}
	writeParam(10,TString::Format("singleEff_bin%d",iBin),    arrPar,   4);  // f1_model_format_0
	writeParam(10,TString::Format("singleEffErr_bin%d",iBin), arrParErr,4);  // f1_model_format_0
   
	canvas.Update();
	canvas.Print(TString::Format("./plots/singleMuonEff_etaBin%d.pdf",iBin));
}
////////////////////////////////////////////////////////////////////////////////////////////
void singleMuondRDimuonPt(int iBin)
{//{{{
//	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	TLorentzVector B_4vec;
	TLorentzVector *Mup_4vec;
	TLorentzVector *Mum_4vec;
    Double_t        Bmass;
    Double_t        Q2;
    Double_t        MupPt;
    Double_t        MumPt;
    Double_t        MupEta;
    Double_t        MumEta;
    Double_t        DimuonPt;
    Double_t        MumudR;
    Int_t           Triggers;
    Double_t        genQ2;
    Double_t        genDimuonPt;
    Double_t        genMumudR;
    Double_t        genMupPt;
    Double_t        genMumPt;
    Double_t        genMupEta;
    Double_t        genMumEta;
    Double_t        genMupPhi;
    Double_t        genMumPhi;
 	
    TBranch        *b_Mum_4vec;   //!
    TBranch        *b_Mup_4vec;   //!
    TBranch        *b_Bmass;   //!
    TBranch        *b_Q2;   //!
    TBranch        *b_MupPt;   //!
    TBranch        *b_MumPt;   //!
    TBranch        *b_MupEta;   //!
    TBranch        *b_MumEta;   //!
    TBranch        *b_DimuonPt;   //!
    TBranch        *b_MumudR;   //!
    TBranch        *b_Triggers;   //!
    TBranch        *b_genQ2;   //!
    TBranch        *b_genDimuonPt;   //!
    TBranch        *b_genMumudR;   //!
    TBranch        *b_genMupPt;   //!
    TBranch        *b_genMumPt;   //!
    TBranch        *b_genMupEta;   //!
    TBranch        *b_genMumEta;   //!
    TBranch        *b_genMupPhi;   //!
    TBranch        *b_genMumPhi;   //!
    TTree    *tree1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    if (!f1 || !f1->IsOpen()) {
        f1 = new TFile("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    }
    f1->GetObject("tree",tree1);
    Mum_4vec = 0;
    Mup_4vec = 0;
    tree1->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
    tree1->SetBranchAddress("Q2", &Q2, &b_Q2);
    tree1->SetBranchAddress("MupPt", &MupPt, &b_MupPt);
    tree1->SetBranchAddress("MumPt", &MumPt, &b_MumPt);
    tree1->SetBranchAddress("MupEta", &MupEta, &b_MupEta);
    tree1->SetBranchAddress("MumEta", &MumEta, &b_MumEta);
    tree1->SetBranchAddress("DimuonPt", &DimuonPt, &b_DimuonPt);
    tree1->SetBranchAddress("MumudR", &MumudR, &b_MumudR);
    tree1->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
    tree1->SetBranchAddress("genQ2", &genQ2, &b_genQ2);
    tree1->SetBranchAddress("genDimuonPt", &genDimuonPt, &b_genDimuonPt);
    tree1->SetBranchAddress("genMumudR", &genMumudR, &b_genMumudR);
    tree1->SetBranchAddress("Mum_4vec", &Mum_4vec, &b_Mum_4vec);
    tree1->SetBranchAddress("Mup_4vec", &Mup_4vec, &b_Mup_4vec);
    tree1->SetBranchAddress("genMupPt", &genMupPt, &b_genMupPt);
    tree1->SetBranchAddress("genMumPt", &genMumPt, &b_genMumPt);
    tree1->SetBranchAddress("genMupEta", &genMupEta, &b_genMupEta);
    tree1->SetBranchAddress("genMumEta", &genMumEta, &b_genMumEta);
    tree1->SetBranchAddress("genMupPhi", &genMupPhi, &b_genMupPhi);
    tree1->SetBranchAddress("genMumPhi", &genMumPhi, &b_genMumPhi);
	Long64_t nentries1 = tree1->GetEntriesFast();
   
//	float thetaLBins[23]={6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
	int PtBins = 11;// The same value as h_acc
    int nLBins = 26;
    TH2D *h2d_2_fine = new TH2D ("h2d_2_fine","reco Efficiency w/ HLT",PtBins,6.,28.,26,0.01,1.31);
    TH2D *h2d_2_acc = new TH2D ("h2d_2_acc","",PtBins,6.,28.,26,0.01,1.31);
    TH2D *h2d_2_reco = new TH2D ("h2d_2_reco","",PtBins,6.,28.,26,0.01,1.31);
	h2d_2_acc->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_reco->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_fine->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
    TH2D *h2d_3_fine = new TH2D ("h2d_3_fine","reco Efficiency w/o HLT",PtBins,6.,28.,26,0.01,1.31);
    TH2D *h2d_3_acc = new TH2D ("h2d_3_acc","",PtBins,6.,28.,26,0.01,1.31);
    TH2D *h2d_3_reco = new TH2D ("h2d_3_reco","",PtBins,6.,28.,26,0.01,1.31);
	h2d_3_acc->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_reco->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_fine->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");

    double TmPt[20][30], TpPt[20][30], TmEta[20][30], TpEta[20][30];
    double AmPt[20][30], ApPt[20][30], AmEta[20][30], ApEta[20][30];
    int Tm[20][30], Tp[20][30];
    double sigM[20][30], sigP[20][30];
    for (int i=0; i<20; i++) {
        for (int j=0; j<30; j++) {
            TmPt[i][j]=0; TpPt[i][j]=0; TmEta[i][j]=0; TpEta[i][j]=0;
            AmPt[i][j]=0; ApPt[i][j]=0; AmEta[i][j]=0; ApEta[i][j]=0;
            Tm[i][j]=0;   Tp[i][j]=0;   sigM[i][j]=0;  sigP[i][j]=0;
        }
    }

    for (Long64_t entry = 0; entry < nentries1; entry++) {
        tree1->GetEntry(entry);
		if (genQ2 > Q2rangeup[iBin] || genQ2 <= Q2rangedn[iBin]) continue;
		if (iBin == 10) {  ///////////////////////////  2015-04-29
		    if (genQ2 < Q2rangeup[3] && genQ2 > Q2rangedn[3]) continue;
		    if (genQ2 < Q2rangeup[5] && genQ2 > Q2rangedn[5]) continue;
	    }
        //if ( fabs(genMumEta) < 2.5 && fabs(genMupEta) < 2.5 && genMumPt > 3.5 && genMumPt > 3.5 && genDimuonPt < 28.){
        if ( fabs(genMumEta) < 0.2 && fabs(genMupEta) < 0.2 && genMumPt > 3.5 && genMumPt > 3.5 && genDimuonPt < 28.){
            h2d_2_acc->Fill(genDimuonPt,genMumudR);
            h2d_3_acc->Fill(genDimuonPt,genMumudR);
		}
		//if (Bmass != -999 && DimuonPt < 28.) { // w/o HLT 2016-04-19
		if (fabs(MumEta) < 0.2 && fabs(MupEta) < 0.2 && Bmass != -999 && DimuonPt < 28.) { // w/o HLT 2016-04-19
            h2d_3_reco->Fill(DimuonPt,MumudR);
		}
		//if (Bmass != -999 && Triggers == 1 && DimuonPt < 28.) { // w/ HLT 2016-04-19
		if (fabs(MumEta) < 0.2 && fabs(MupEta) < 0.2 && Bmass != -999 && Triggers == 1 && DimuonPt < 28.) { // w/ HLT 2016-04-19
            h2d_2_reco->Fill(DimuonPt,MumudR);
            for (int i=1; i<=14; i++) {
                for (int j=1; j<=26; j++) {
                    if (DimuonPt > (6.0 + i * 2.0) && DimuonPt < (6.0 + (i+1) * 2.0) && MumudR > (0.01 + j * 0.05) && MumudR < (0.01 + (j+1) * 0.05)) {
                        TmPt[i][j]=TmPt[i][j]+fabs(MumPt); TpPt[i][j]=TpPt[i][j]+fabs(MupPt); TmEta[i][j]=TmEta[i][j]+fabs(MumEta); TpEta[i][j]=TpEta[i][j]+fabs(MupEta);
                        Tm[i][j]=Tm[i][j]+1;   Tp[i][j]=Tp[i][j]+1;   
                    }
                }
            }
        }
    }
	TFile f("../RootFiles/singleMuonEff/singleMuonEff_noTracking_L3ptg2_final.root", "read");
    TF1* user_f0 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin0");
    TF1* user_f1 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin1");
    TF1* user_f2 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin2");
    TF1* user_f3 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin3");
    TF1* user_f4 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin4");
    TF1* user_f5 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin5");
    TF1* user_f6 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin6");
    TF1* user_f7 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin7");
    TF1* user_f8 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin8");
    TF1* user_f9 = (TF1*)f.Get("fitTotEff_DATA_pt_etaBin9");
    for (int i=1; i<14; i++) {
        for (int j=1; j<26; j++) {
            AmPt[i][j]=TmPt[i][j]/Tm[i][j]; ApPt[i][j]=TpPt[i][j]/Tp[i][j]; AmEta[i][j]=TmEta[i][j]/Tm[i][j]; ApEta[i][j]=TpEta[i][j]/Tp[i][j];
            if (fabs(AmEta[i][j]) < 0.2)  sigM[i][j]= user_f0->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 0.2 && fabs(AmEta[i][j]) < 0.3) sigM[i][j]= user_f1->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 0.3 && fabs(AmEta[i][j]) < 0.6) sigM[i][j]= user_f2->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 0.6 && fabs(AmEta[i][j]) < 0.8) sigM[i][j]= user_f3->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 0.8 && fabs(AmEta[i][j]) < 1.0) sigM[i][j]= user_f4->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 1.0 && fabs(AmEta[i][j]) < 1.2) sigM[i][j]= user_f5->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 1.2 && fabs(AmEta[i][j]) < 1.4) sigM[i][j]= user_f6->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 1.4 && fabs(AmEta[i][j]) < 1.6) sigM[i][j]= user_f7->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 1.6 && fabs(AmEta[i][j]) < 1.8) sigM[i][j]= user_f8->Eval(AmPt[i][j]);
            if (fabs(AmEta[i][j]) > 1.8 && fabs(AmEta[i][j]) < 2.0) sigM[i][j]= user_f9->Eval(AmPt[i][j]);
            if (fabs(ApEta[i][j]) < 0.2)  sigP[i][j]= user_f0->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 0.2 && fabs(ApEta[i][j]) < 0.3) sigP[i][j]= user_f1->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 0.3 && fabs(ApEta[i][j]) < 0.6) sigP[i][j]= user_f2->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 0.6 && fabs(ApEta[i][j]) < 0.8) sigP[i][j]= user_f3->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 0.8 && fabs(ApEta[i][j]) < 1.0) sigP[i][j]= user_f4->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 1.0 && fabs(ApEta[i][j]) < 1.2) sigP[i][j]= user_f5->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 1.2 && fabs(ApEta[i][j]) < 1.4) sigP[i][j]= user_f6->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 1.4 && fabs(ApEta[i][j]) < 1.6) sigP[i][j]= user_f7->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 1.6 && fabs(ApEta[i][j]) < 1.8) sigP[i][j]= user_f8->Eval(ApPt[i][j]);
            if (fabs(ApEta[i][j]) > 1.8 && fabs(ApEta[i][j]) < 2.0) sigP[i][j]= user_f9->Eval(ApPt[i][j]);
        }
    }
/*   
    cout<<"2, 5 - f6: mPt = "<<user_f6->Eval(4.676)<<" [ "<<user_f6->Eval(4.676-0.6337)<<", "<<user_f6->Eval(4.676+0.6337)<<"]"<<endl;
    cout<<"2, 5 - f6: pPt = "<<user_f6->Eval(4.687)<<" [ "<<user_f6->Eval(4.687-0.6541)<<", "<<user_f6->Eval(4.687+0.6541)<<"]"<<endl;
    cout<<"3, 5 - f4: mPt = "<<user_f4->Eval(5.536)<<" [ "<<user_f4->Eval(5.536-1.0030)<<", "<<user_f4->Eval(5.536+1.0030)<<"]"<<endl;
    cout<<"3, 5 - f4: pPt = "<<user_f4->Eval(5.587)<<" [ "<<user_f4->Eval(5.587-0.9998)<<", "<<user_f4->Eval(5.587+0.9998)<<"]"<<endl;
    cout<<"4, 5 - f4: mPt = "<<user_f4->Eval(6.509)<<" [ "<<user_f4->Eval(6.509-1.4130)<<", "<<user_f4->Eval(6.509+1.4130)<<"]"<<endl;
    cout<<"4, 5 - f4: pPt = "<<user_f4->Eval(6.561)<<" [ "<<user_f4->Eval(6.561-1.4190)<<", "<<user_f4->Eval(6.561+1.4190)<<"]"<<endl;
    cout<<"5, 5 - f4: mPt = "<<user_f4->Eval(7.527)<<" [ "<<user_f4->Eval(7.527-1.8500)<<", "<<user_f4->Eval(7.527+1.8500)<<"]"<<endl;
    cout<<"5, 5 - f4: pPt = "<<user_f4->Eval(7.520)<<" [ "<<user_f4->Eval(7.520-1.8580)<<", "<<user_f4->Eval(7.520+1.8580)<<"]"<<endl;
 
    cout<<"2, 5 - f5: mPt = "<<user_f5->Eval(4.676)<<" [ "<<user_f5->Eval(4.676-0.6337)<<", "<<user_f5->Eval(4.676+0.6337)<<"]"<<endl;
    cout<<"2, 5 - f5: pPt = "<<user_f5->Eval(4.687)<<" [ "<<user_f5->Eval(4.687-0.6541)<<", "<<user_f5->Eval(4.687+0.6541)<<"]"<<endl;
    cout<<"3, 5 - f3: mPt = "<<user_f3->Eval(5.536)<<" [ "<<user_f3->Eval(5.536-1.0030)<<", "<<user_f3->Eval(5.536+1.0030)<<"]"<<endl;
    cout<<"3, 5 - f3: pPt = "<<user_f3->Eval(5.587)<<" [ "<<user_f3->Eval(5.587-0.9998)<<", "<<user_f3->Eval(5.587+0.9998)<<"]"<<endl;
    cout<<"4, 5 - f3: mPt = "<<user_f3->Eval(6.509)<<" [ "<<user_f3->Eval(6.509-1.4130)<<", "<<user_f3->Eval(6.509+1.4130)<<"]"<<endl;
    cout<<"4, 5 - f3: pPt = "<<user_f3->Eval(6.561)<<" [ "<<user_f3->Eval(6.561-1.4190)<<", "<<user_f3->Eval(6.561+1.4190)<<"]"<<endl;
    cout<<"5, 5 - f3: mPt = "<<user_f3->Eval(7.527)<<" [ "<<user_f3->Eval(7.527-1.8500)<<", "<<user_f3->Eval(7.527+1.8500)<<"]"<<endl;
    cout<<"5, 5 - f3: pPt = "<<user_f3->Eval(7.520)<<" [ "<<user_f3->Eval(7.520-1.8580)<<", "<<user_f3->Eval(7.520+1.8580)<<"]"<<endl;
 
    cout<<"2, 5 - f7: mPt = "<<user_f7->Eval(4.676)<<" [ "<<user_f7->Eval(4.676-0.6337)<<", "<<user_f7->Eval(4.676+0.6337)<<"]"<<endl;
    cout<<"2, 5 - f7: pPt = "<<user_f7->Eval(4.687)<<" [ "<<user_f7->Eval(4.687-0.6541)<<", "<<user_f7->Eval(4.687+0.6541)<<"]"<<endl;
    cout<<"3, 5 - f5: mPt = "<<user_f5->Eval(5.536)<<" [ "<<user_f5->Eval(5.536-1.0030)<<", "<<user_f5->Eval(5.536+1.0030)<<"]"<<endl;
    cout<<"3, 5 - f5: pPt = "<<user_f5->Eval(5.587)<<" [ "<<user_f5->Eval(5.587-0.9998)<<", "<<user_f5->Eval(5.587+0.9998)<<"]"<<endl;
    cout<<"4, 5 - f5: mPt = "<<user_f5->Eval(6.509)<<" [ "<<user_f5->Eval(6.509-1.4130)<<", "<<user_f5->Eval(6.509+1.4130)<<"]"<<endl;
    cout<<"4, 5 - f5: pPt = "<<user_f5->Eval(6.561)<<" [ "<<user_f5->Eval(6.561-1.4190)<<", "<<user_f5->Eval(6.561+1.4190)<<"]"<<endl;
    cout<<"5, 5 - f5: mPt = "<<user_f5->Eval(7.527)<<" [ "<<user_f5->Eval(7.527-1.8500)<<", "<<user_f5->Eval(7.527+1.8500)<<"]"<<endl;
    cout<<"5, 5 - f5: pPt = "<<user_f5->Eval(7.520)<<" [ "<<user_f5->Eval(7.520-1.8580)<<", "<<user_f5->Eval(7.520+1.8580)<<"]"<<endl;
*/
	int bindR = 26;
//	Draw
    TCanvas canvas("canvas");
//	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	TPaveText* paveText = new TPaveText( 0.20, 0.77, 0.30, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(0);
    paveText->AddText(Form("bin %d ", iBin));
    TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	
    if (iBin != 100) {	   
        TH2D *h2d_sigM_eff = new TH2D ("h2d_sigM_eff","title",PtBins,6.,28.,26,0.01,1.31);
	    h2d_sigM_eff->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	    h2d_sigM_eff->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	    h2d_sigM_eff->SetTitle("#mu^{-} efficiency");
        TH2D *h2d_sigP_eff = new TH2D ("h2d_sigP_eff","title",PtBins,6.,28.,26,0.01,1.31);
	    h2d_sigP_eff->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	    h2d_sigP_eff->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	    h2d_sigP_eff->SetTitle("#mu^{+} efficiency");
	    for (int i = 1; i <= PtBins; i++) {
	        for (int j = 1; j <= bindR; j++) {
		        if (h2d_2_acc->GetBinContent(i,j) == 0 || h2d_2_reco->GetBinContent(i,j) == 0) {
			        h2d_sigM_eff->SetBinContent(i,j,0.);
			        h2d_sigP_eff->SetBinContent(i,j,0.);
		        }else{
			        h2d_sigM_eff->SetBinContent(i,j,sigM[i][j]);
			        h2d_sigP_eff->SetBinContent(i,j,sigP[i][j]);
		        }
	        }
	    }
//      cout<<h2d_sigP_eff->GetBinContent(11,2)<<endl;
//      cout<<h2d_sigM_eff->GetBinContent(11,2)<<endl;
   
	    for (int i = 1; i <= PtBins; i++) {
	        for (int j = 1; j <= bindR; j++) {
		        if (h2d_2_acc->GetBinContent(i,j) == 0 || h2d_2_reco->GetBinContent(i,j) == 0 || sigM[i][j] == 0 || sigP[i][j] == 0) {
			        h2d_2_fine->SetBinContent(i,j,0.);
			        h2d_2_fine->SetBinError(i,j,0.);
		        }else{
			        h2d_2_fine->SetBinContent(i,j,h2d_2_reco->GetBinContent(i,j)/h2d_2_acc->GetBinContent(i,j)/sigM[i][j]/sigP[i][j]);
			        h2d_2_fine->SetBinError(i,j,0.);
		        }
	        }
	    }
	    for (int i = 1; i <= PtBins; i++) {
	        for (int j = 1; j <= bindR; j++) {
		        if (h2d_3_acc->GetBinContent(i,j) == 0 || h2d_3_reco->GetBinContent(i,j) == 0) {
			        h2d_3_fine->SetBinContent(i,j,0.);
			        h2d_3_fine->SetBinError(i,j,0.);
		        }else if (h2d_3_acc->GetBinContent(i,j) < h2d_3_reco->GetBinContent(i,j)) {
			        h2d_3_fine->SetBinContent(i,j,0.);
			        h2d_3_fine->SetBinError(i,j,0.);
		        }else{
			        h2d_3_fine->SetBinContent(i,j,h2d_3_reco->GetBinContent(i,j)/h2d_3_acc->GetBinContent(i,j));
			        h2d_3_fine->SetBinError(i,j,sqrt(h2d_3_fine->GetBinContent(i,j)*(1-h2d_3_fine->GetBinContent(i,j))/h2d_3_acc->GetBinContent(i,j)));
		        }
	        }
	    }
        TH2D *h2d_4_ratio = new TH2D ("h2d_4_ratio","title",PtBins,6.,28.,26,0.01,1.31);
	    h2d_4_ratio->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	    h2d_4_ratio->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	    h2d_4_ratio->SetTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)] ");
	    for (int i = 1; i <= PtBins; i++) {
	        for (int j = 1; j <= bindR; j++) {
		        if (h2d_2_fine->GetBinContent(i,j) == 0 || h2d_3_fine->GetBinContent(i,j) == 0 || sigM[i][j] == 0 || sigP[i][j] == 0 ) {
			        h2d_4_ratio->SetBinContent(i,j,0.);
			        h2d_4_ratio->SetBinError(i,j,0.);
		        }else{
			        h2d_4_ratio->SetBinContent(i,j,h2d_2_fine->GetBinContent(i,j)/h2d_3_fine->GetBinContent(i,j)/sigM[i][j]/sigP[i][j]);
			        h2d_4_ratio->SetBinError(i,j,0.);
		        }
	        }
	    }
        h2d_sigM_eff->Draw("COLZ");
        h2d_sigM_eff->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
        canvas.Update();
        canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_sigM_2D_eff_bin%d.pdf",iBin));
	    h2d_sigP_eff->Draw("COLZ");
        h2d_sigP_eff->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
        canvas.Update();
        canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_sigP_2D_eff_bin%d.pdf",iBin));

        h2d_2_acc->Draw("COLZ");
        h2d_2_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_naccL_2D_wHLT_bin%d.pdf",iBin));
//	Draw #events pass all cuts
//	    h2d_2_reco->Draw("scat=0.5");
	    h2d_2_reco->Draw("COLZ");
        h2d_2_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_nrecoL_2D_wHLT_bin%d.pdf",iBin));
	    h2d_3_acc->Draw("COLZ");
        h2d_3_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_naccL_2D_woHLT_bin%d.pdf",iBin));
//	    h2d_3_reco->Draw("scat=0.5");
	    h2d_3_reco->Draw("COLZ");
        h2d_3_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
        paveText->Draw();
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_nrecoL_2D_woHLT_bin%d.pdf",iBin));

	    h2d_2_fine->SetStats(0);
        //h2d_2_fine->SetContour((sizeof(levels2)/sizeof(Double_t)), levels2);
        //h2d_2_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
        h2d_2_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
        h2d_2_fine->Draw("COLZ");
        paveText->Draw();
	    t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	    t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	    if (iBin != 10 && iBin != 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
        } else if (iBin == 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
        }
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_recoL_2D_wHLT_bin%d.pdf",iBin));
	
        h2d_3_fine->SetStats(0);
        //h2d_3_fine->SetContour((sizeof(levels3)/sizeof(Double_t)), levels3);
        //h2d_3_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
        h2d_3_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
        h2d_3_fine->Draw("COLZZ");
        paveText->Draw();
	    t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	    t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	    if (iBin != 10 && iBin != 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
        } else if (iBin == 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
        }
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_recoL_2D_woHLT_bin%d.pdf",iBin));
	
	    h2d_4_ratio->SetStats(0);
        h2d_4_ratio->SetTitleSize(0.04,"y");
        h2d_4_ratio->GetYaxis()->SetRangeUser(0.0, 1.31);
	    h2d_4_ratio->Draw("COLZ");
        paveText->Draw();
	    t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	    t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	    if (iBin != 10 && iBin != 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
        } else if (iBin == 11) {
            t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
        }
	    canvas.Update();
	    canvas.Print(TString::Format("./plots/singleMuondRDimuonPt_recoL_2D_ratio_HLT_bin%d.pdf",iBin));
    }

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
void recoEffQ2dR(int iBin)
{//
//	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	TLorentzVector B_4vec;
	TLorentzVector *Mup_4vec;
	TLorentzVector *Mum_4vec;
    Double_t        Bmass;
    Double_t        Q2;
    Double_t        MumudR;
    Int_t           Triggers;
    Double_t        genQ2;
    Double_t        genMumudR;
    Double_t        genMupPt;
    Double_t        genMumPt;
    Double_t        genMupEta;
    Double_t        genMumEta;
    Double_t        genMupPhi;
    Double_t        genMumPhi;
 	
    TBranch        *b_Mum_4vec;   //!
    TBranch        *b_Mup_4vec;   //!
    TBranch        *b_Bmass;   //!
    TBranch        *b_Q2;   //!
    TBranch        *b_MumudR;   //!
    TBranch        *b_Triggers;   //!
    TBranch        *b_genQ2;   //!
    TBranch        *b_genMumudR;   //!
    TBranch        *b_genMupPt;   //!
    TBranch        *b_genMumPt;   //!
    TBranch        *b_genMupEta;   //!
    TBranch        *b_genMumEta;   //!
    TBranch        *b_genMupPhi;   //!
    TBranch        *b_genMumPhi;   //!
    TTree    *tree1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    if (!f1 || !f1->IsOpen()) {
        f1 = new TFile("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    }
    f1->GetObject("tree",tree1);
    Mum_4vec = 0;
    Mup_4vec = 0;
    tree1->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
    tree1->SetBranchAddress("Q2", &Q2, &b_Q2);
    tree1->SetBranchAddress("MumudR", &MumudR, &b_MumudR);
    tree1->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
    tree1->SetBranchAddress("genQ2", &genQ2, &b_genQ2);
    tree1->SetBranchAddress("genMumudR", &genMumudR, &b_genMumudR);
    tree1->SetBranchAddress("Mum_4vec", &Mum_4vec, &b_Mum_4vec);
    tree1->SetBranchAddress("Mup_4vec", &Mup_4vec, &b_Mup_4vec);
    tree1->SetBranchAddress("genMupPt", &genMupPt, &b_genMupPt);
    tree1->SetBranchAddress("genMumPt", &genMumPt, &b_genMumPt);
    tree1->SetBranchAddress("genMupEta", &genMupEta, &b_genMupEta);
    tree1->SetBranchAddress("genMumEta", &genMumEta, &b_genMumEta);
    tree1->SetBranchAddress("genMupPhi", &genMupPhi, &b_genMupPhi);
    tree1->SetBranchAddress("genMumPhi", &genMumPhi, &b_genMumPhi);
	Long64_t nentries1 = tree1->GetEntriesFast();
   
//	float thetaLBins[23]={6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
    int nLBins = 26;
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,0,1.3); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,0,1.3);
	TH1F h3_nreco_fine("h3_nreco_fine","h3_nreco_fine",nLBins,0,1.3);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.026");
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.026");
//	h2_nreco_fine.SetTitle("#Delta R distribution w/o HLT");
	h3_nreco_fine.SetMinimum(0.);
	h3_nreco_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h3_nreco_fine.SetYTitle("Reconstructed Events/0.026");
//	h3_nreco_fine.SetTitle("#Delta R distribution w/o HLT");
	
    int NLBins = 77;
	TH1F h2_nacc_q2("h2_nacc_q2" ,"h2_nacc_q2" ,NLBins,0.,23.1); 
	TH1F h2_nreco_q2("h2_nreco_q2","h2_nreco_q2",NLBins,0.,23.1);
	TH1F h3_nreco_q2("h3_nreco_q2","h3_nreco_q2",NLBins,0.,23.1);
	h2_nacc_q2.SetStats(0);
	h2_nacc_q2.SetMinimum(0.);
	h2_nacc_q2.SetXTitle("q^{2}[(GeV)^{2}]");
	h2_nacc_q2.SetYTitle("Acceptated Events/0.3");
	h2_nreco_q2.SetStats(0);
	h2_nreco_q2.SetMinimum(0.);
	h2_nreco_q2.SetXTitle("q^{2}[(GeV)^{2}]");
	h2_nreco_q2.SetYTitle("Reconstructed Events/0.3");
	h3_nreco_q2.SetStats(0);
	h3_nreco_q2.SetMinimum(0.);
	h3_nreco_q2.SetXTitle("q^{2}[(GeV)^{2}]");
	h3_nreco_q2.SetYTitle("Reconstructed Events/0.3");
   
    int Q2Bins[12]={10,23,44,14,28,22,20,20,40,22,55,60};
    TH2D *h2d_2_fine = new TH2D ("h2d_2_fine","reco Efficiency w/ HLT",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
    TH2D *h2d_2_acc = new TH2D ("h2d_2_acc","",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
    TH2D *h2d_2_reco = new TH2D ("h2d_2_reco","",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
	h2d_2_acc->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_2_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_reco->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_2_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_fine->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_2_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
    TH2D *h2d_3_fine = new TH2D ("h2d_3_fine","reco Efficiency w/o HLT",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
    TH2D *h2d_3_acc = new TH2D ("h2d_3_acc","",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
    TH2D *h2d_3_reco = new TH2D ("h2d_3_reco","",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
	h2d_3_acc->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_3_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_reco->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_3_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_fine->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_3_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");

    //double dR = 0.; double dEta = 0.; double dPhi = 0.;
    //double gdR = 0.; double gdEta = 0.; double gdPhi = 0.;
//  double dddr1 =0.2; double dddr2 =0.2;
//  double dddr3 =0.2; double dddr4 =0.2;
//  double dddr5 =0.2; double dddr6 =0.2;
    for (Long64_t entry = 0; entry < nentries1; entry++) {
        tree1->GetEntry(entry);
		if (genQ2 > Q2rangeup[iBin] || genQ2 <= Q2rangedn[iBin]) continue;
		if (iBin == 10) {  ///////////////////////////  2015-04-29
			if (genQ2 < Q2rangeup[3] && genQ2 > Q2rangedn[3]) continue;
			if (genQ2 < Q2rangeup[5] && genQ2 > Q2rangedn[5]) continue;
		}
		
        if ( fabs(genMumEta) < 2.5 && fabs(genMupEta) < 2.5 && genMumPt > 3.5 && genMumPt > 3.5 ){
            h2_nacc_fine.Fill(genMumudR);
            h2_nacc_q2.Fill(genQ2);
            h2d_2_acc->Fill(genQ2,genMumudR);
            h2d_3_acc->Fill(genQ2,genMumudR);
		}
		if (Bmass != -999 ) { // w/o HLT 2016-04-19
			h3_nreco_fine.Fill(MumudR);
			h3_nreco_q2.Fill(genQ2);
            //h2d_3_reco->Fill(Q2,MumudR);
            h2d_3_reco->Fill(genQ2,MumudR);
		}
		if (Bmass != -999 && Triggers == 1 ) { // w/ HLT 2016-04-19
            //if (MumudR > (2*TMath::Pi())) cout<<2*TMath::Pi()<<" < dR = "<<MumudR<<" ; dEta ="<<Mup_4vec->Eta()-Mum_4vec->Eta()<<" ; dPhi = "<<Mup_4vec->Phi()-Mum_4vec->Phi()<<endl;
			h2_nreco_fine.Fill(MumudR);
			h2_nreco_q2.Fill(genQ2);
            //h2d_2_reco->Fill(Q2,MumudR);
            h2d_2_reco->Fill(genQ2,MumudR);
		}
    }
//	cout<<dddr1<<"  "<<dddr2<<endl;
//	cout<<dddr3<<"  "<<dddr4<<endl;
//	cout<<dddr5<<"  "<<dddr6<<endl;
	int bindR = 26;
/*   
    TH1F h2_reco_fine("h2_reco_fine","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
		    //printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
		    h2_reco_fine.SetBinContent(i,0.);
		    //h2_reco_fine.SetBinError(i,1.);
		    h2_reco_fine.SetBinError(i,0.);
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
		    //printf("INFO: wiht HLT, recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	}
	TH1F h3_reco_fine("h3_reco_fine","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h3_nreco_fine.GetBinContent(i) == 0) {
		    //printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h3_reco_fine.SetBinContent(i,0.);
		    //h3_reco_fine.SetBinError(i,1.);
			h3_reco_fine.SetBinError(i,0.);
		}else{
			h3_reco_fine.SetBinContent(i,h3_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h3_reco_fine.SetBinError(i,sqrt(h3_reco_fine.GetBinContent(i)*(1-h3_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
		    //printf("INFO: without HLT, recoEfficiency_fine(%d)=%f +- %f.\n",i,h3_reco_fine.GetBinContent(i),h3_reco_fine.GetBinError(i));
		}
	}
	TH1F h4_reco_ratio("h4_reco_ratio","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_reco_fine.GetBinContent(i) == 0 || h3_reco_fine.GetBinContent(i) == 0) {
		    //printf("WARNING: Ratio(%d)=0, set error to be 1.\n",i);
			h4_reco_ratio.SetBinContent(i,0.);
		    //h4_reco_ratio.SetBinError(i,1.);
			h4_reco_ratio.SetBinError(i,0.);
		}else{
			h4_reco_ratio.SetBinContent(i,h2_reco_fine.GetBinContent(i)/h3_reco_fine.GetBinContent(i));
			h4_reco_ratio.SetBinError(i,h4_reco_ratio.GetBinContent(i)*sqrt(pow( h3_reco_fine.GetBinError(i)/h3_reco_fine.GetBinContent(i), 2 ) + pow( h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i), 2)));
		    //printf("INFO: Ratio(%d)=%f +- %f.\n",i,h4_reco_ratio.GetBinContent(i),h4_reco_ratio.GetBinError(i));
		}
	}
*/
//	Draw
    TCanvas canvas("canvas");
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
//	TPaveText* paveText = new TPaveText( 0.45, 0.62, 0.54, 0.72, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(0);
    paveText->AddText(Form("bin %d ", iBin));
    TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);	

/*	
if (iBin == 10 || iBin == 11) {	
	for (int i = 1; i <= Q2Bins[iBin]; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_2_acc->GetBinContent(i,j) == 0 || h2d_2_reco->GetBinContent(i,j) == 0) {
			    h2d_2_fine->SetBinContent(i,j,0.);
			    h2d_2_fine->SetBinError(i,j,0.);
		    }else{
			    h2d_2_fine->SetBinContent(i,j,h2d_2_reco->GetBinContent(i,j)/h2d_2_acc->GetBinContent(i,j));
			    h2d_2_fine->SetBinError(i,j,sqrt(h2d_2_fine->GetBinContent(i,j)*(1-h2d_2_fine->GetBinContent(i,j))/h2d_2_acc->GetBinContent(i,j)));
		    }
	    }
	}
	for (int i = 1; i <= Q2Bins[iBin]; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_3_acc->GetBinContent(i,j) == 0 || h2d_3_reco->GetBinContent(i,j) == 0) {
			    h2d_3_fine->SetBinContent(i,j,0.);
			    h2d_3_fine->SetBinError(i,j,0.);
		    }else if (h2d_3_acc->GetBinContent(i,j) < h2d_3_reco->GetBinContent(i,j)) {
			    h2d_3_fine->SetBinContent(i,j,0.);
			    h2d_3_fine->SetBinError(i,j,0.);
		    }else{
			    h2d_3_fine->SetBinContent(i,j,h2d_3_reco->GetBinContent(i,j)/h2d_3_acc->GetBinContent(i,j));
			    h2d_3_fine->SetBinError(i,j,sqrt(h2d_3_fine->GetBinContent(i,j)*(1-h2d_3_fine->GetBinContent(i,j))/h2d_3_acc->GetBinContent(i,j)));
		    }
	    }
	}
    TH2D *h2d_4_ratio = new TH2D ("h2d_4_ratio","title",Q2Bins[iBin],Q2rangedn[iBin],Q2rangeup[iBin],26,0.01,1.31);
	h2d_4_ratio->GetXaxis()->SetTitle("q^{2}[(GeV)^{2}]");
	h2d_4_ratio->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_4_ratio->SetTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)] ");
	for (int i = 1; i <= Q2Bins[iBin]; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_2_fine->GetBinContent(i,j) == 0 || h2d_3_fine->GetBinContent(i,j) == 0) {
			    h2d_4_ratio->SetBinContent(i,j,0.);
			    h2d_4_ratio->SetBinError(i,j,0.);
		    }else{
			    h2d_4_ratio->SetBinContent(i,j,h2d_2_fine->GetBinContent(i,j)/h2d_3_fine->GetBinContent(i,j));
			    h2d_4_ratio->SetBinError(i,j,h2d_4_ratio->GetBinContent(i,j)*sqrt(pow( h2d_3_fine->GetBinError(i,j)/h2d_3_fine->GetBinContent(i,j), 2 ) + pow( h2d_2_fine->GetBinError(i,j)/h2d_2_fine->GetBinContent(i,j), 2)));
		    }
	    }
	}

//	Draw #events in acceptance
//	gStyle->SetPalette("kBird");
    Int_t colors[] = { 1, 9, 4, 7, 8, 3, 5, 6, 46, 2 }; // #colors >= #levels - 1
    //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    // #levels <= #colors + 1 (notes: +-3.4e38 = +-FLT_MAX; +1.17e-38 = +FLT_MIN)
    Double_t levels[] = {1.17e-38, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 3.4e38};
    Double_t levels2[]= {1.17e-38, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.60, 3.4e38};
    Double_t levels3[]= {1.17e-38, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 3.4e38};
    //Double_t levels2[]= {1.17e-38, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 3.4e38};
    //Double_t levels3[]= {1.17e-38, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 3.4e38};

    h2d_2_acc->Draw("COLZ");
    h2d_2_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_naccL_2D_wHLT_bin%d.pdf",iBin));
//	Draw #events pass all cuts
//	h2d_2_reco->Draw("scat=0.5");
	h2d_2_reco->Draw("COLZ");
    h2d_2_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_2D_wHLT_bin%d.pdf",iBin));
	h2d_3_acc->Draw("COLZ");
    h2d_3_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_naccL_2D_woHLT_bin%d.pdf",iBin));
//	h2d_3_reco->Draw("scat=0.5");
	h2d_3_reco->Draw("COLZ");
    h2d_3_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_2D_woHLT_bin%d.pdf",iBin));

	h2d_2_fine->SetStats(0);
    //h2d_2_fine->SetContour((sizeof(levels2)/sizeof(Double_t)), levels2);
    //h2d_2_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
    h2d_2_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
    h2d_2_fine->Draw("COLZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_2D_wHLT_bin%d.pdf",iBin));
	
    h2d_3_fine->SetStats(0);
    //h2d_3_fine->SetContour((sizeof(levels3)/sizeof(Double_t)), levels3);
    //h2d_3_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
    h2d_3_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
    h2d_3_fine->Draw("COLZZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_2D_woHLT_bin%d.pdf",iBin));
	
	h2d_4_ratio->SetStats(0);
    h2d_4_ratio->SetTitleSize(0.04,"y");
    h2d_4_ratio->GetYaxis()->SetRangeUser(0.0, 1.31);
	h2d_4_ratio->Draw("COLZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_2D_ratio_HLT_bin%d.pdf",iBin));

	h2d_4_ratio->SetStats(0);
    h2d_4_ratio->SetTitleSize(0.04,"y");
    h2d_4_ratio->GetYaxis()->SetRangeUser(0.0, 1.31);
	h2d_4_ratio->Draw("SURF3");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_3D_ratio_HLT_bin%d.pdf",iBin));
}
*/

//	Draw #events in acceptance
	h2_nacc_fine.Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_naccL_fine_bin%d.pdf",iBin));
//	Draw #events pass all cuts
	h2_nreco_fine.SetStats(0);
	h3_nreco_fine.SetStats(0);
	h3_nreco_fine.Draw("PE1");
	h2_nreco_fine.Draw("PE1 SAME");
	h2_nreco_fine.SetMarkerColor(2);
	h2_nreco_fine.SetLineColor(2);
	h3_nreco_fine.SetMarkerColor(4);
	h3_nreco_fine.SetLineColor(4);
	h2_nreco_fine.SetMarkerStyle(20);
	h3_nreco_fine.SetMarkerStyle(24);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        paveText->Draw();
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    } else if (iBin == 10) {
        paveText->Draw();
    }
	TLegend *leg =new TLegend(0.73,0.75,0.87,0.87,NULL,"brNDC");
	leg->AddEntry("h2_nreco_fine"," w/ HLT "," P ");
	leg->AddEntry("h3_nreco_fine"," w/o HLT "," P ");
	leg->SetLineColor(1);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
//	canvas.Update();
//	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_fine_HLT_bin%d.pdf",iBin));
//	h3_nreco_fine.Draw();
//  paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_fine_bin%d.pdf",iBin));
/*	
	h2_reco_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h2_reco_fine.SetYTitle("Reco-Efficiency/0.05");
	h2_reco_fine.SetMinimum(0.);
	h3_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetMaximum(0.40);
	h3_reco_fine.SetMaximum(0.40);
	if (iBin == 3 || iBin == 5) {
		h2_reco_fine.SetMaximum(0.025);
		h3_reco_fine.SetMaximum(0.025);
	}
	h2_reco_fine.Draw("PE1");
	h3_reco_fine.Draw("PE1 SAME");
	h2_reco_fine.SetMarkerColor(2);
	h2_reco_fine.SetLineColor(2);
	h3_reco_fine.SetMarkerColor(4);
	h3_reco_fine.SetLineColor(4);
	h2_reco_fine.SetMarkerStyle(20);
	h3_reco_fine.SetMarkerStyle(24);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        paveText->Draw();
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    } else if (iBin == 10) {
        paveText->Draw();
    }
	TLegend *leg =new TLegend(0.73,0.75,0.87,0.87,NULL,"brNDC");
	leg->AddEntry("h2_reco_fine"," w/ HLT "," P ");
	leg->AddEntry("h3_reco_fine"," w/o HLT "," P ");
	leg->SetLineColor(1);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_fine_HLT_bin%d.pdf",iBin));

	h4_reco_ratio.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
//	h4_reco_ratio.SetYTitle("ratio/2.0[GeV]");
	h4_reco_ratio.SetYTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)]  ");
	h4_reco_ratio.SetMinimum(0.);
	h4_reco_ratio.SetMaximum(1.);
	h4_reco_ratio.SetTitleSize(0.04,"y");
	h4_reco_ratio.Draw("PE1");
	h4_reco_ratio.SetMarkerColor(1);
	h4_reco_ratio.SetLineColor(1);
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_ratio_HLT_bin%d.pdf",iBin));

if (iBin == 10 || iBin == 11) {	
    TH1F h2_reco_q2("h2_reco_q2","",NLBins,0.,23.1);
	for (int i = 1; i <= NLBins; i++) {
		if (h2_nacc_q2.GetBinContent(i) == 0 || h2_nreco_q2.GetBinContent(i) == 0) {
			h2_reco_q2.SetBinContent(i,0.);
			h2_reco_q2.SetBinError(i,0.);
		}else{
			h2_reco_q2.SetBinContent(i,h2_nreco_q2.GetBinContent(i)/h2_nacc_q2.GetBinContent(i));
			h2_reco_q2.SetBinError(i,sqrt(h2_reco_q2.GetBinContent(i)*(1-h2_reco_q2.GetBinContent(i))/h2_nacc_q2.GetBinContent(i)));
		}
	}
	TH1F h3_reco_q2("h3_reco_q2","",NLBins,0.,23.1);
	for (int i = 1; i <= NLBins; i++) {
		if (h2_nacc_q2.GetBinContent(i) == 0 || h3_nreco_q2.GetBinContent(i) == 0) {
			h3_reco_q2.SetBinContent(i,0.);
			h3_reco_q2.SetBinError(i,0.);
		}else{
			h3_reco_q2.SetBinContent(i,h3_nreco_q2.GetBinContent(i)/h2_nacc_q2.GetBinContent(i));
			h3_reco_q2.SetBinError(i,sqrt(h3_reco_q2.GetBinContent(i)*(1-h3_reco_q2.GetBinContent(i))/h2_nacc_q2.GetBinContent(i)));
		}
	}
	TH1F h4_reco_q2_ratio("h4_reco_q2_ratio","",NLBins,0.,23.1);
	for (int i = 1; i <= NLBins; i++) {
		if (h2_reco_q2.GetBinContent(i) == 0 || h3_reco_q2.GetBinContent(i) == 0) {
			h4_reco_q2_ratio.SetBinContent(i,0.);
			h4_reco_q2_ratio.SetBinError(i,0.);
		}else{
			h4_reco_q2_ratio.SetBinContent(i,h2_reco_q2.GetBinContent(i)/h3_reco_q2.GetBinContent(i));
			h4_reco_q2_ratio.SetBinError(i,h4_reco_q2_ratio.GetBinContent(i)*sqrt(pow( h3_reco_q2.GetBinError(i)/h3_reco_q2.GetBinContent(i), 2 ) + pow( h2_reco_q2.GetBinError(i)/h2_reco_q2.GetBinContent(i), 2)));
		}
	}
   
    h2_reco_q2.SetXTitle("q^{2}[(GeV)^{2}]");
	h2_reco_q2.SetYTitle("Reco-Efficiency/0.3");
	h2_reco_q2.SetMinimum(0.);
	h3_reco_q2.SetMinimum(0.);
	h2_reco_q2.SetMaximum(0.12);
	h3_reco_q2.SetMaximum(0.12);
	h2_reco_q2.Draw("PE1");
	h3_reco_q2.Draw("PE1 SAME");
	h2_reco_q2.SetMarkerColor(2);
	h2_reco_q2.SetLineColor(2);
	h3_reco_q2.SetMarkerColor(4);
	h3_reco_q2.SetLineColor(4);
	h2_reco_q2.SetMarkerStyle(20);
	h3_reco_q2.SetMarkerStyle(24);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	TLegend *leg2 =new TLegend(0.73,0.75,0.87,0.87,NULL,"brNDC");
	leg2->AddEntry("h2_reco_q2"," w/ HLT "," P ");
	leg2->AddEntry("h3_reco_q2"," w/o HLT "," P ");
	leg2->SetLineColor(1);
	leg2->SetFillColor(0);
	leg2->SetTextSize(0.02);
	leg2->Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_q2_HLT_bin%d.pdf",iBin));
	
	h4_reco_q2_ratio.SetXTitle("q^{2}[(GeV)^{2}]");
	h4_reco_q2_ratio.SetYTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)] /0.3 ");
	h4_reco_q2_ratio.SetMinimum(0.);
	h4_reco_q2_ratio.SetMaximum(1.);
	h4_reco_q2_ratio.SetTitleSize(0.04,"y");
	h4_reco_q2_ratio.Draw("PE1");
	h4_reco_q2_ratio.SetMarkerColor(1);
	h4_reco_q2_ratio.SetLineColor(1);
    paveText->Draw();
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_recoL_q2_ratio_HLT_bin%d.pdf",iBin));
}
*/
}//
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
void recoEffdRDimuonPt(int iBin)
{//{{{
//	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	TLorentzVector B_4vec;
	TLorentzVector *Mup_4vec;
	TLorentzVector *Mum_4vec;
    Double_t        Bmass;
    Double_t        Q2;
    Double_t        MupPt;
    Double_t        MumPt;
    Double_t        MupEta;
    Double_t        MumEta;
    Double_t        DimuonPt;
    Double_t        MumudR;
    Int_t           Triggers;
    Double_t        genQ2;
    Double_t        genDimuonPt;
    Double_t        genMumudR;
    Double_t        genMupPt;
    Double_t        genMumPt;
    Double_t        genMupEta;
    Double_t        genMumEta;
    Double_t        genMupPhi;
    Double_t        genMumPhi;
 	
    TBranch        *b_Mum_4vec;   //!
    TBranch        *b_Mup_4vec;   //!
    TBranch        *b_Bmass;   //!
    TBranch        *b_Q2;   //!
    TBranch        *b_MupPt;   //!
    TBranch        *b_MumPt;   //!
    TBranch        *b_MupEta;   //!
    TBranch        *b_MumEta;   //!
    TBranch        *b_DimuonPt;   //!
    TBranch        *b_MumudR;   //!
    TBranch        *b_Triggers;   //!
    TBranch        *b_genQ2;   //!
    TBranch        *b_genDimuonPt;   //!
    TBranch        *b_genMumudR;   //!
    TBranch        *b_genMupPt;   //!
    TBranch        *b_genMumPt;   //!
    TBranch        *b_genMupEta;   //!
    TBranch        *b_genMumEta;   //!
    TBranch        *b_genMupPhi;   //!
    TBranch        *b_genMumPhi;   //!
    TTree    *tree1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    if (!f1 || !f1->IsOpen()) {
        f1 = new TFile("../RootFiles/Files/MC_Signal_8TeV_v4_cut0+Q+B2+resonance-1+OPT+Anti3_HLT.root");
    }
    f1->GetObject("tree",tree1);
    Mum_4vec = 0;
    Mup_4vec = 0;
    tree1->SetBranchAddress("Bmass", &Bmass, &b_Bmass);
    tree1->SetBranchAddress("Q2", &Q2, &b_Q2);
    tree1->SetBranchAddress("MupPt", &MupPt, &b_MupPt);
    tree1->SetBranchAddress("MumPt", &MumPt, &b_MumPt);
    tree1->SetBranchAddress("MupEta", &MupEta, &b_MupEta);
    tree1->SetBranchAddress("MumEta", &MumEta, &b_MumEta);
    tree1->SetBranchAddress("DimuonPt", &DimuonPt, &b_DimuonPt);
    tree1->SetBranchAddress("MumudR", &MumudR, &b_MumudR);
    tree1->SetBranchAddress("Triggers", &Triggers, &b_Triggers);
    tree1->SetBranchAddress("genQ2", &genQ2, &b_genQ2);
    tree1->SetBranchAddress("genDimuonPt", &genDimuonPt, &b_genDimuonPt);
    tree1->SetBranchAddress("genMumudR", &genMumudR, &b_genMumudR);
    tree1->SetBranchAddress("Mum_4vec", &Mum_4vec, &b_Mum_4vec);
    tree1->SetBranchAddress("Mup_4vec", &Mup_4vec, &b_Mup_4vec);
    tree1->SetBranchAddress("genMupPt", &genMupPt, &b_genMupPt);
    tree1->SetBranchAddress("genMumPt", &genMumPt, &b_genMumPt);
    tree1->SetBranchAddress("genMupEta", &genMupEta, &b_genMupEta);
    tree1->SetBranchAddress("genMumEta", &genMumEta, &b_genMumEta);
    tree1->SetBranchAddress("genMupPhi", &genMupPhi, &b_genMupPhi);
    tree1->SetBranchAddress("genMumPhi", &genMumPhi, &b_genMumPhi);
	Long64_t nentries1 = tree1->GetEntriesFast();
   
//	float thetaLBins[23]={6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
	int PtBins = 22;// The same value as h_acc
	TH1F h2_nacc_Pt("h2_nacc_Pt" ,"h2_nacc_Pt" ,PtBins,6, 50); 
	TH1F h2_nreco_Pt("h2_nreco_Pt","h2_nreco_Pt",PtBins,6, 50);
	TH1F h3_nreco_Pt("h3_nreco_Pt","h3_nreco_Pt",PtBins,6, 50);
	h2_nacc_Pt.SetStats(0);
	h2_nacc_Pt.SetMinimum(0.);
	h2_nacc_Pt.SetXTitle(" p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2_nacc_Pt.SetYTitle("Acceptated Events/2.0[GeV]");
	h2_nreco_Pt.SetStats(0);
	h2_nreco_Pt.SetMinimum(0.);
	h2_nreco_Pt.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2_nreco_Pt.SetYTitle("Reconstructed Events/2.0[GeV]");
	h3_nreco_Pt.SetStats(0);
	h3_nreco_Pt.SetMinimum(0.);
	h3_nreco_Pt.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h3_nreco_Pt.SetYTitle("Reconstructed Events/2.0[GeV]");
	h3_nreco_Pt.SetTitle("p^{#mu^{+}#mu^{-}}_{T} distribution");
   
    int nLBins = 26;
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,0,1.3); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,0,1.3);
	TH1F h3_nreco_fine("h3_nreco_fine","h3_nreco_fine",nLBins,0,1.3);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.026");
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.026");
//	h2_nreco_fine.SetTitle("#Delta R distribution w/o HLT");
	h3_nreco_fine.SetMinimum(0.);
	h3_nreco_fine.SetXTitle("#Delta R(#mu^{+}#mu^{-})");
	h3_nreco_fine.SetYTitle("Reconstructed Events/0.026");
    h3_nreco_fine.SetTitle("#Delta R distribution");
	
    TH2D *h2d_2_fine = new TH2D ("h2d_2_fine","reco Efficiency w/ HLT",PtBins,6.,50,26,0.01,1.31);
    TH2D *h2d_2_acc = new TH2D ("h2d_2_acc","",PtBins,6.,50,26,0.01,1.31);
    TH2D *h2d_2_reco = new TH2D ("h2d_2_reco","",PtBins,6.,50,26,0.01,1.31);
	h2d_2_acc->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_reco->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_2_fine->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_2_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
    TH2D *h2d_3_fine = new TH2D ("h2d_3_fine","reco Efficiency w/o HLT",PtBins,6.,50,26,0.01,1.31);
    TH2D *h2d_3_acc = new TH2D ("h2d_3_acc","",PtBins,6.,50,26,0.01,1.31);
    TH2D *h2d_3_reco = new TH2D ("h2d_3_reco","",PtBins,6.,50,26,0.01,1.31);
	h2d_3_acc->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_acc->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_reco->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_reco->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_3_fine->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_3_fine->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");

    for (Long64_t entry = 0; entry < nentries1; entry++) {
        tree1->GetEntry(entry);
		if (genQ2 > Q2rangeup[iBin] || genQ2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
			    if (genQ2 < Q2rangeup[3] && genQ2 > Q2rangedn[3]) continue;
				if (genQ2 < Q2rangeup[5] && genQ2 > Q2rangedn[5]) continue;
			}		
        if ( fabs(genMumEta) < 2.5 && fabs(genMupEta) < 2.5 && genMumPt > 3.5 && genMumPt > 3.5 ){
            h2_nacc_fine.Fill(genMumudR);
            h2_nacc_Pt.Fill(genDimuonPt);
            h2d_2_acc->Fill(genDimuonPt,genMumudR);
            h2d_3_acc->Fill(genDimuonPt,genMumudR);
		}
		if (Bmass != -999 ) { // w/o HLT 2016-04-19
			h3_nreco_fine.Fill(MumudR);
			h3_nreco_Pt.Fill(DimuonPt);
            h2d_3_reco->Fill(DimuonPt,MumudR);
		}
		if (Bmass != -999 && Triggers == 1 ) { // w/ HLT 2016-04-19
			h2_nreco_fine.Fill(MumudR);
			h2_nreco_Pt.Fill(DimuonPt);
            h2d_2_reco->Fill(DimuonPt,MumudR);
		}
    }
	int bindR = 26;
/*   
    TH1F h2_reco_fine("h2_reco_fine","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,0.);
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
		}
	}
	TH1F h3_reco_fine("h3_reco_fine","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h3_nreco_fine.GetBinContent(i) == 0) {
			h3_reco_fine.SetBinContent(i,0.);
			h3_reco_fine.SetBinError(i,0.);
		}else{
			h3_reco_fine.SetBinContent(i,h3_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h3_reco_fine.SetBinError(i,sqrt(h3_reco_fine.GetBinContent(i)*(1-h3_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
		}
	}
	TH1F h4_reco_ratio("h4_reco_ratio","",nLBins,0,1.3);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_reco_fine.GetBinContent(i) == 0 || h3_reco_fine.GetBinContent(i) == 0) {
			h4_reco_ratio.SetBinContent(i,0.);
			h4_reco_ratio.SetBinError(i,0.);
		}else{
			h4_reco_ratio.SetBinContent(i,h2_reco_fine.GetBinContent(i)/h3_reco_fine.GetBinContent(i));
			h4_reco_ratio.SetBinError(i,h4_reco_ratio.GetBinContent(i)*sqrt(pow( h3_reco_fine.GetBinError(i)/h3_reco_fine.GetBinContent(i), 2 ) + pow( h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i), 2)));
		}
	}
*/
//	Draw
    TCanvas canvas("canvas");
//	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	TPaveText* paveText = new TPaveText( 0.20, 0.77, 0.30, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(0);
    paveText->AddText(Form("bin %d ", iBin));
    TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
/*	
//	Draw #events pass all cuts
	h2_nreco_fine.SetStats(0);
	h3_nreco_fine.SetStats(0);
	h3_nreco_fine.Draw("PE1");
	h2_nreco_fine.Draw("PE1 SAME");
	h2_nreco_fine.SetMarkerColor(2);
	h2_nreco_fine.SetLineColor(2);
	h3_nreco_fine.SetMarkerColor(4);
	h3_nreco_fine.SetLineColor(4);
	h2_nreco_fine.SetMarkerStyle(20);
	h3_nreco_fine.SetMarkerStyle(24);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        paveText->Draw();
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    } else if (iBin == 10) {
        t1->DrawLatex(.28,.81,TString::Format(" bin %d ",iBin));
        //paveText->Draw();
    }
	TLegend *leg =new TLegend(0.73,0.75,0.87,0.87,NULL,"brNDC");
	leg->AddEntry("h2_nreco_fine"," w/ HLT "," P ");
	leg->AddEntry("h3_nreco_fine"," w/o HLT "," P ");
	leg->SetLineColor(1);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_fine_bin%d.pdf",iBin));

//	Draw #events pass all cuts
	h2_nreco_Pt.SetStats(0);
	h3_nreco_Pt.SetStats(0);
	h3_nreco_Pt.Draw("PE1");
	h2_nreco_Pt.Draw("PE1 SAME");
	h2_nreco_Pt.SetMarkerColor(2);
	h2_nreco_Pt.SetLineColor(2);
	h3_nreco_Pt.SetMarkerColor(4);
	h3_nreco_Pt.SetLineColor(4);
	h2_nreco_Pt.SetMarkerStyle(20);
	h3_nreco_Pt.SetMarkerStyle(24);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        paveText->Draw();
        t1->DrawLatex(.28,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    } else if (iBin == 10) {
        t1->DrawLatex(.28,.81,TString::Format(" bin %d ",iBin));
        //paveText->Draw();
    }
	TLegend *leg1 =new TLegend(0.73,0.75,0.87,0.87,NULL,"brNDC");
	leg1->AddEntry("h2_nreco_Pt"," w/ HLT "," P ");
	leg1->AddEntry("h3_nreco_Pt"," w/o HLT "," P ");
	leg1->SetLineColor(1);
	leg1->SetFillColor(0);
	leg1->SetTextSize(0.02);
	leg1->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffQ2dR_nrecoL_Pt_bin%d.pdf",iBin));
*/
	
if (iBin != 100) {	
//if (iBin == 10 || iBin == 11) {	
	for (int i = 1; i <= PtBins; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_2_acc->GetBinContent(i,j) == 0 || h2d_2_reco->GetBinContent(i,j) == 0 ) {
			    h2d_2_fine->SetBinContent(i,j,0.);
			    h2d_2_fine->SetBinError(i,j,0.);
		    }else{
			    h2d_2_fine->SetBinContent(i,j,h2d_2_reco->GetBinContent(i,j)/h2d_2_acc->GetBinContent(i,j));
			    h2d_2_fine->SetBinError(i,j,sqrt(h2d_2_fine->GetBinContent(i,j)*(1-h2d_2_fine->GetBinContent(i,j))/h2d_2_acc->GetBinContent(i,j)));
		    }
	    }
	}
    cout<<"2, 5 wHLT: "<<h2d_2_fine->GetBinContent(2,5)<<" +- "<<h2d_2_fine->GetBinError(2,5)<<endl;
    cout<<"3, 5 wHLT: "<<h2d_2_fine->GetBinContent(3,5)<<" +- "<<h2d_2_fine->GetBinError(3,5)<<endl;
    cout<<"4, 5 wHLT: "<<h2d_2_fine->GetBinContent(4,5)<<" +- "<<h2d_2_fine->GetBinError(4,5)<<endl;
    cout<<"5, 5 wHLT: "<<h2d_2_fine->GetBinContent(5,5)<<" +- "<<h2d_2_fine->GetBinError(5,5)<<endl;
	for (int i = 1; i <= PtBins; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_3_acc->GetBinContent(i,j) == 0 || h2d_3_reco->GetBinContent(i,j) == 0) {
			    h2d_3_fine->SetBinContent(i,j,0.);
			    h2d_3_fine->SetBinError(i,j,0.);
		    }else if (h2d_3_acc->GetBinContent(i,j) < h2d_3_reco->GetBinContent(i,j)) {
			    h2d_3_fine->SetBinContent(i,j,0.);
			    h2d_3_fine->SetBinError(i,j,0.);
		    }else{
			    h2d_3_fine->SetBinContent(i,j,h2d_3_reco->GetBinContent(i,j)/h2d_3_acc->GetBinContent(i,j));
			    h2d_3_fine->SetBinError(i,j,sqrt(h2d_3_fine->GetBinContent(i,j)*(1-h2d_3_fine->GetBinContent(i,j))/h2d_3_acc->GetBinContent(i,j)));
		    }
	    }
	}
    cout<<"2, 5 woHLT: "<<h2d_3_fine->GetBinContent(2,5)<<" +- "<<h2d_3_fine->GetBinError(2,5)<<endl;
    cout<<"3, 5 woHLT: "<<h2d_3_fine->GetBinContent(3,5)<<" +- "<<h2d_3_fine->GetBinError(3,5)<<endl;
    cout<<"4, 5 woHLT: "<<h2d_3_fine->GetBinContent(4,5)<<" +- "<<h2d_3_fine->GetBinError(4,5)<<endl;
    cout<<"5, 5 woHLT: "<<h2d_3_fine->GetBinContent(5,5)<<" +- "<<h2d_3_fine->GetBinError(5,5)<<endl;
    TH2D *h2d_4_ratio = new TH2D ("h2d_4_ratio","title",PtBins,6.,50,26,0.01,1.31);
	h2d_4_ratio->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2d_4_ratio->GetYaxis()->SetTitle("#Delta R(#mu^{+}#mu^{-})");
	h2d_4_ratio->SetTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)] ");
	for (int i = 1; i <= PtBins; i++) {
	    for (int j = 1; j <= bindR; j++) {
		    if (h2d_2_fine->GetBinContent(i,j) == 0 || h2d_3_fine->GetBinContent(i,j) == 0 ) {
			    h2d_4_ratio->SetBinContent(i,j,0.);
			    h2d_4_ratio->SetBinError(i,j,0.);
		    }else{
			    h2d_4_ratio->SetBinContent(i,j,h2d_2_fine->GetBinContent(i,j)/h2d_3_fine->GetBinContent(i,j));
			    h2d_4_ratio->SetBinError(i,j,h2d_4_ratio->GetBinContent(i,j)*sqrt(pow( h2d_3_fine->GetBinError(i,j)/h2d_3_fine->GetBinContent(i,j), 2 ) + pow( h2d_2_fine->GetBinError(i,j)/h2d_2_fine->GetBinContent(i,j), 2)));
		    }
	    }
	}
    cout<<"2, 5 ratio: "<<h2d_4_ratio->GetBinContent(2,5)<<" +- "<<h2d_4_ratio->GetBinError(2,5)<<endl;
    cout<<"3, 5 ratio: "<<h2d_4_ratio->GetBinContent(3,5)<<" +- "<<h2d_4_ratio->GetBinError(3,5)<<endl;
    cout<<"4, 5 ratio: "<<h2d_4_ratio->GetBinContent(4,5)<<" +- "<<h2d_4_ratio->GetBinError(4,5)<<endl;
    cout<<"5, 5 ratio: "<<h2d_4_ratio->GetBinContent(5,5)<<" +- "<<h2d_4_ratio->GetBinError(5,5)<<endl;
   
//	Draw #events in acceptance
/*
//	gStyle->SetPalette("kBird");
    Int_t colors[] = { 1, 9, 4, 7, 8, 3, 5, 6, 46, 2 }; // #colors >= #levels - 1
    //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    // #levels <= #colors + 1 (notes: +-3.4e38 = +-FLT_MAX; +1.17e-38 = +FLT_MIN)
    Double_t levels[] = {1.17e-38, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 3.4e38};
    Double_t levels2[]= {1.17e-38, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.60, 3.4e38};
    Double_t levels3[]= {1.17e-38, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 3.4e38};
*/
    h2d_2_acc->Draw("COLZ");
    h2d_2_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_naccL_2D_wHLT_bin%d.pdf",iBin));
//	Draw #events pass all cuts
//	h2d_2_reco->Draw("scat=0.5");
	h2d_2_reco->Draw("COLZ");
    h2d_2_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_nrecoL_2D_wHLT_bin%d.pdf",iBin));
	h2d_3_acc->Draw("COLZ");
    h2d_3_acc->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_naccL_2D_woHLT_bin%d.pdf",iBin));
//	h2d_3_reco->Draw("scat=0.5");
	h2d_3_reco->Draw("COLZ");
    h2d_3_reco->GetYaxis()->SetRangeUser(0.0, 1.31);
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_nrecoL_2D_woHLT_bin%d.pdf",iBin));

	h2d_2_fine->SetStats(0);
    //h2d_2_fine->SetContour((sizeof(levels2)/sizeof(Double_t)), levels2);
    //h2d_2_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
    h2d_2_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
    h2d_2_fine->Draw("COLZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_recoL_2D_wHLT_bin%d.pdf",iBin));
	
    h2d_3_fine->SetStats(0);
    //h2d_3_fine->SetContour((sizeof(levels3)/sizeof(Double_t)), levels3);
    //h2d_3_fine->GetZaxis()->SetRangeUser(0.0, 0.60);
    h2d_3_fine->GetYaxis()->SetRangeUser(0.0, 1.31);
    h2d_3_fine->Draw("COLZZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_recoL_2D_woHLT_bin%d.pdf",iBin));
	
	h2d_4_ratio->SetStats(0);
    h2d_4_ratio->SetTitleSize(0.04,"y");
    h2d_4_ratio->GetYaxis()->SetRangeUser(0.0, 1.31);
	h2d_4_ratio->Draw("COLZ");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	if (iBin != 10 && iBin != 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < %.2f", Q2rangedn[iBin], Q2rangeup[iBin]));
    } else if (iBin == 11) {
        t1->DrawLatex(.32,.81,TString::Format(" %.2f < q^{2} < 23.04", Q2rangedn[iBin]));
    }
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_recoL_2D_ratio_HLT_bin%d.pdf",iBin));

	h2d_4_ratio->SetStats(0);
    h2d_4_ratio->SetTitleSize(0.04,"y");
    h2d_4_ratio->GetYaxis()->SetRangeUser(0.0, 1.31);
	h2d_4_ratio->Draw("SURF3");
    paveText->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffdRDimuonPt_recoL_3D_ratio_HLT_bin%d.pdf",iBin));

}

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
void recoEffDimuonPt(int iBin)
{//{{{
	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	// No HLT
//	double effUpperBound2[11] = {    0.3,   0.24,   0.16,   0.015,   0.13,  0.015,  0.12,   0.13,   0.11,   0.20,    0.13};
	// with HLT
//	double effUpperBound2[11] = {   0.10,   0.08,   0.07,   0.004,   0.06,   0.006,   0.06,   0.06,   0.06,   0.08,   0.06};
//	TLorentzVector B_4vec;
//	double Bctau = 0;
	int    Triggers = 0;
	double BMass = 0;
	double DimuonPt = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double genDimuonPt = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	
	ch->SetBranchStatus("Triggers"      , 1);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("DimuonPt"      , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genDimuonPt"   , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("Triggers"     , &Triggers);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("DimuonPt"     , &DimuonPt);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("genDimuonPt"  , &genDimuonPt);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
	
//	Fill histograms
//	float thetaLBins[21]={5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45};
//	float thetaLBins[21]={4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,300};
	float thetaLBins[23]={6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
//	TH1F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins); 
	int nLBins = 22;// The same value as h_acc
//	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,5.,1); 
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,thetaLBins); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,thetaLBins);
	TH1F h3_nreco_fine("h3_nreco_fine","h3_nreco_fine",nLBins,thetaLBins);
	h2_nacc_fine.SetStats(0);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle(" p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2_nacc_fine.SetYTitle("Acceptated Events/2.0[GeV]");
	h2_nreco_fine.SetStats(0);
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2_nreco_fine.SetYTitle("Reconstructed Events/2.0[GeV]");
	h3_nreco_fine.SetStats(0);
	h3_nreco_fine.SetMinimum(0.);
	h3_nreco_fine.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h3_nreco_fine.SetYTitle("Reconstructed Events/2.0[GeV]");
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
			h2_nacc_fine.Fill(genDimuonPt);
		}
		if (BMass != -999 ) { // w/o HLT 2016-04-19
			h3_nreco_fine.Fill(DimuonPt);
		}
		if (BMass != -999 && Triggers == 1 ) { // w/ HLT 2016-04-19
			h2_nreco_fine.Fill(DimuonPt);
		}
	}
	
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,thetaLBins);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,1.);
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: wiht HLT, recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	}
	TH1F h3_reco_fine("h3_reco_fine","",nLBins,thetaLBins);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h3_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h3_reco_fine.SetBinContent(i,0.);
			h3_reco_fine.SetBinError(i,1.);
		}else{
			h3_reco_fine.SetBinContent(i,h3_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h3_reco_fine.SetBinError(i,sqrt(h3_reco_fine.GetBinContent(i)*(1-h3_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: without HLT, recoEfficiency_fine(%d)=%f +- %f.\n",i,h3_reco_fine.GetBinContent(i),h3_reco_fine.GetBinError(i));
		}
	}
	TH1F h4_reco_ratio("h4_reco_ratio","",nLBins,thetaLBins);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_reco_fine.GetBinContent(i) == 0 || h3_reco_fine.GetBinContent(i) == 0) {
			printf("WARNING: Ratio(%d)=0, set error to be 1.\n",i);
			h4_reco_ratio.SetBinContent(i,0.);
			h4_reco_ratio.SetBinError(i,1.);
		}else{
			h4_reco_ratio.SetBinContent(i,h2_reco_fine.GetBinContent(i)/h3_reco_fine.GetBinContent(i));
			h4_reco_ratio.SetBinError(i,h4_reco_ratio.GetBinContent(i)*sqrt(pow( h3_reco_fine.GetBinError(i)/h3_reco_fine.GetBinContent(i), 2 ) + pow( h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i), 2)));
			printf("INFO: Ratio(%d)=%f +- %f.\n",i,h4_reco_ratio.GetBinContent(i),h4_reco_ratio.GetBinError(i));
		}
	}
	
	TString f1_model_format_1 ;
	f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
//	Draw #events in acceptance
	h2_nacc_fine.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffDimuonPt_naccL_fine_bin%d.pdf",iBin));
//	Draw #events pass all cuts
	h2_nreco_fine.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffDimuonPt_nrecoL_fine_HLT_bin%d.pdf",iBin));
	h3_nreco_fine.Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffDimuonPt_nrecoL_fine_bin%d.pdf",iBin));
/*
//	Draw FitResult for recoEfficiency
	const int nPar_r = 7;
	TF1 *f1_model_r = new TF1 ("f1_model_r", f1_model_format_1, -1., 1.);
	f1_model_r->SetParameter(0,0.);
	f1_model_r->SetParameter(1,0.01);
	f1_model_r->SetParameter(2,0.01);
	f1_model_r->SetParameter(3,0.01);
	f1_model_r->SetParameter(4,0.01);
	f1_model_r->SetParameter(5,0.01);
	f1_model_r->SetParameter(6,0.01);  // f1_model_format_1
		
	h2_reco_fine.Fit(f1_model_r,"R"); 
   
	h2_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetTitleOffset(1.15,"Y");
	h2_reco_fine.SetXTitle("cos#theta_{l}^{reco}");
	h2_reco_fine.SetYTitle("Reco-Efficiency / 0.1");
	h2_reco_fine.SetStats(0);
	h2_reco_fine.SetMaximum(effUpperBound2[iBin]);
	h2_reco_fine.Draw("PE1");
	h2_reco_fine.Draw();
	f1_model_r->SetTitle("");
	f1_model_r->SetLineWidth(2);
	f1_model_r->SetLineColor(2);
	f1_model_r->Draw(" SAME ");
*/	
	
	h2_reco_fine.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
	h2_reco_fine.SetYTitle("Reco-Efficiency/2.0[GeV]");
	h2_reco_fine.SetMinimum(0.);
	h3_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetMaximum(0.35);
	h3_reco_fine.SetMaximum(0.35);
	if (iBin == 3 || iBin == 5) {
		h2_reco_fine.SetMaximum(0.03);
		h3_reco_fine.SetMaximum(0.03);
	}
	h2_reco_fine.Draw("PE1");
	h3_reco_fine.Draw("PE1 SAME");
	h2_reco_fine.SetMarkerColor(2);
	h2_reco_fine.SetLineColor(2);
	h3_reco_fine.SetMarkerColor(4);
	h3_reco_fine.SetLineColor(4);
	h2_reco_fine.SetMarkerStyle(20);
	h3_reco_fine.SetMarkerStyle(24);
	TPaveText* paveText = new TPaveText( 0.16, 0.57, 0.26, 0.67, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	TLegend *leg =new TLegend(0.16,0.75,0.30,0.87,NULL,"brNDC");
	leg->AddEntry("h2_reco_fine"," w/ HLT "," P ");
	leg->AddEntry("h3_reco_fine"," w/o HLT "," P ");
	leg->SetLineColor(1);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffDimuonPt_recoL_fine_HLT_bin%d.pdf",iBin));
	
	h4_reco_ratio.SetXTitle("p^{#mu^{+}#mu^{-}}_{T} [GeV]");
//	h4_reco_ratio.SetYTitle("ratio/2.0[GeV]");
	h4_reco_ratio.SetYTitle("ratio of efficiency: [(w/ HLT) / (w/o HLT)] /2.0[GeV] ");
	h4_reco_ratio.SetMinimum(0.);
	h4_reco_ratio.SetMaximum(1.);
	h4_reco_ratio.SetTitleSize(0.04,"y");
	h4_reco_ratio.Draw("PE1");
	h4_reco_ratio.SetMarkerColor(1);
	h4_reco_ratio.SetLineColor(1);
	TPaveText* paveText1 = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText1->SetBorderSize(0);
	paveText1->SetFillColor(19);
    paveText1->AddText(Form("bin %d ", iBin));
    paveText1->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	canvas.Update();
	canvas.Print(TString::Format("./plots/recoEffDimuonPt_recoL_ratio_HLT_bin%d.pdf",iBin));
/*
//	Save fitting results
	double chi2Val_r=0;
	double arrPar_r[nPar_r], arrParErr_r[nPar_r];
	for (int iPar = 0; iPar < nPar_r; iPar++) {
		arrPar_r[iPar]    = f1_model_r->GetParameter(iPar);
		arrParErr_r[iPar] = f1_model_r->GetParError(iPar);
		chi2Val_r         = f1_model_r->GetChisquare();
	}
	std::vector<double> output_r;
	for (int iPar = 0; iPar < nPar_r; iPar++){
		output_r.push_back(arrPar_r[iPar]);
		output_r.push_back(arrParErr_r[iPar]);
		printf("%18.15f,",arrPar_r[iPar]);
		if (iPar+1 >= nPar_r) printf("\n");
	}
	for (int i = 0; i < output_r.size(); i=i+2) {
		printf("%18.15f,",output_r[i+1]);
		if (i+2 >= output_r.size()) printf("\n");
	}
*/
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEff(int iBin)
//std::vector<double> accXrecoEff(int iBin)
{//{{{
	setTDRStyle();
//	TH1::SetDefaultSumw2();
//	ROOT::Math::Minimizer::SetDefaultMaxFunctionCalls(10000);
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	// No HLT
//	double effUpperBound[11]  = { 3.0e-3, 2.8e-3, 2.4e-3, 0.15e-3, 2.2e-3, 0.3e-3, 2.8e-3, 3.5e-3, 3.5e-3, 2.6e-3, 2.2e-3};
//	double effUpperBound2[11] = {    0.3,   0.24,   0.16,   0.015,   0.13,  0.015,  0.12,   0.13,   0.11,   0.20,    0.13};
	// with HLT
//	double effUpperBound[11]  = { 1.5e-3, 1.4e-3, 1.2e-3, 0.08e-3, 1.1e-3, 0.14e-3, 1.4e-3, 1.8e-3, 1.8e-3, 1.3e-3, 1.1e-3};
//	double effUpperBound2[11] = {   0.10,   0.09,   0.07,   0.004,   0.06,   0.006,   0.06,   0.06,   0.06,   0.08,   0.06};
//	TLorentzVector B_4vec;
//	double Bctau = 0;
	int    Triggers = 0;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double gBEta  = 0;
	double gCosThetaL = 0;
	double CosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
	int    Count  = 0;
	int    Count1 = 0;
	int    Count2 = 0;
//	ch->SetBranchStatus("*",0);
//	ch->SetBranchStatus("B_4vec"        , 1);
//	ch->SetBranchStatus("Bctau"         , 1);
	ch->SetBranchStatus("Triggers"      , 1);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchStatus("genBEta"       , 1);
//	ch->SetBranchStatus("B_4vec"        , &B_4vec);
//	ch->SetBranchStatus("Bctau"         , &Bctau);
	ch->SetBranchAddress("Triggers"     , &Triggers);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	ch->SetBranchAddress("genBEta"      , &gBEta);
	
//	Load acceptance
	TFile f_acc("./RootFiles/acceptance_8TeV.root");
	TFile *fout = new TFile(TString::Format("./RootFiles/Efficiency_%d.root",iBin),"RECREATE");
//	TH1F *h1_acc       = (TH1F*)f_acc.Get(TString::Format("h1_acc_bin%d",      iBin));
	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
//	TH1F *h1_ngen      = (TH1F*)f_acc.Get(TString::Format("h1_ngen_bin%d",     iBin));
//	TH1F *h1_ngen_fine = (TH1F*)f_acc.Get(TString::Format("h1_ngen_fine_bin%d",iBin));
	
//	Fill histograms
//	float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
//	TH1F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins); 
//	TH1F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins);
	int nLBins = 20;// The same value as h_acc
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,-1,1); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,-1,1);
	h2_nacc_fine.SetStats(0);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.1");
	h2_nreco_fine.SetStats(0);
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("cos#theta_{l}^{reco}");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.1");
	
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
	//	cout<<Triggers<<" ";
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
		if (iBin == 10) {  ///////////////////////////  2015-04-29
			if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
			if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
		}
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
	//	if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 2.8 && gmuppt > 2.8 ){
			Count+=1;
            if (fabs(gBEta) > 2.5) Count1++;
//			h2_nacc.Fill(gCosThetaL);
			h2_nacc_fine.Fill(gCosThetaL);
		}
	//	if (BMass != -999 && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)) ){
	//	if (BMass != -999 ){    ////////////////////////   12-10 N.A.
		if (BMass != -999 && Triggers == 1) { // HLT 2016-04-19
//			h2_nreco.Fill(CosThetaL);
			h2_nreco_fine.Fill(CosThetaL);
			Count2+=1;
		}
	}
	cout<<"Bin "<<iBin<<" = "<<Count<<"; BEta>2.5: "<<Count1<<endl;
    cout<<"Reco_"<<iBin<<" = "<<Count2<<" / "<<Count<<" = "<<1.0*Count2/Count<<endl;
//	Calculate efficiency
//	TH1F h2_eff("h2_eff","",6,thetaLBins);
//	for (int i = 1; i <= 6; i++) {//L
//	//	Build from MC samples
//		if (h2_nacc.GetBinContent(i) == 0 || h2_nreco.GetBinContent(i) == 0) {
//			printf("WARNING: Efficiency(%d)0, set error to be 1.\n",i);
//			h2_eff.SetBinContent(i,0.);
//			h2_eff.SetBinError(i,1.);
//		}else{
//			h2_eff.SetBinContent(i,h2_nreco.GetBinContent(i)/h2_nacc.GetBinContent(i) * h1_acc->GetBinContent(i));
//			h2_eff.SetBinError(i,h2_eff.GetBinContent(i)*sqrt(-1./h2_nacc.GetBinContent(i)+1./h2_nreco.GetBinContent(i)+pow(h1_acc->GetBinError(i)/h1_acc->GetBinContent(i),2)));
//			printf("INFO: Efficiency(%d) : %f +- %f.\n",i,h2_eff.GetBinContent(i),h2_eff.GetBinError(i));
//		}
//	}
//	
//	double reco_w[20], eff_w[20];
//	double sum_reco_w_L = 0., sum_eff_w_L = 0., sum_reco_wx_L = 0., sum_eff_wx_L = 0.;
//	double sum_reco_w_R = 0., sum_eff_w_R = 0., sum_reco_wx_R = 0., sum_eff_wx_R = 0.;
//	double A_reco_L, A_eff_L, A_recoerr_L, A_efferr_L;
//	double A_reco_R, A_eff_R, A_recoerr_R, A_efferr_R;
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
	//	if (!(Selection[1])) continue;   //////////////////////////////  12-08
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,1.);
		//	h2_eff_fine.SetBinError(i,0.001e-3); ////////////////////////////////////////////////////////////////// test
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	//	if ( i == 20 && iBin == 0) { cout<<h2_reco_fine.GetBinContent(i)<<endl<<h2_nreco_fine.GetBinContent(i)<<endl<<h2_nacc_fine.GetBinContent(i)<<endl; }
//		reco_w[i-1] = 1. / ( h2_reco_fine.GetBinError(i) * h2_reco_fine.GetBinError(i) );
//		if (i <= 10) {
//			sum_reco_w_L = sum_reco_w_L + reco_w[i-1];
//			sum_reco_wx_L = sum_reco_wx_L + reco_w[i-1] * h2_reco_fine.GetBinContent(i);
//		} else {
//			sum_reco_w_R = sum_reco_w_R + reco_w[i-1];
//			sum_reco_wx_R = sum_reco_wx_R + reco_w[i-1] * h2_reco_fine.GetBinContent(i);
//		}	
	}
//	A_reco_L = sum_reco_wx_L / sum_reco_w_L;
//	A_reco_R = sum_reco_wx_R / sum_reco_w_R;
//	A_recoerr_L = 1. / sqrt( sum_reco_w_L );
//	A_recoerr_R = 1. / sqrt( sum_reco_w_R );
//	cout<<endl<<"A_reco_L["<<iBin<<"] = "<<A_reco_L<<" +- "<<A_recoerr_L<<endl;
//	cout<<"A_reco_R["<<iBin<<"] = "<<A_reco_R<<" +- "<<A_recoerr_R<<endl<<endl;
	
	TH1F h2_eff_fine("h2_eff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
	//	if (!(Selection[1])) continue;   //////////////////////////////  12-08
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_eff_fine.SetBinContent(i,0.);
			h2_eff_fine.SetBinError(i,1.);
		//	h2_eff_fine.SetBinContent(i,1.e-6);
		//	h2_eff_fine.SetBinError(i,5.e-7);
		//	h2_eff_fine.SetBinError(i,0.001e-3); ////////////////////////////////////////////////////////////////// test
		}else{
			h2_eff_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(pow(h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i),2)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
		//	h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(-1./h2_nacc_fine.GetBinContent(i)+1./h2_nreco_fine.GetBinContent(i)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
		//	h2_eff_fine.SetBinContent(i,h2_reco_fine.GetBinContent(i)*h1_acc_fine->GetBinContent(i));
		//	h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i));
			printf("INFO: Efficiency_fine(%d)=%f +- %f.\n",i,h2_eff_fine.GetBinContent(i),h2_eff_fine.GetBinError(i));
		}
//		eff_w[i-1] = 1. / ( h2_eff_fine.GetBinError(i) * h2_eff_fine.GetBinError(i) );
//		if (i <= 10) {
//			sum_eff_w_L = sum_eff_w_L + eff_w[i-1];
//			sum_eff_wx_L = sum_eff_wx_L + eff_w[i-1] * h2_eff_fine.GetBinContent(i);
//		} else {
//			sum_eff_w_R = sum_eff_w_R + eff_w[i-1];
//			sum_eff_wx_R = sum_eff_wx_R + eff_w[i-1] * h2_eff_fine.GetBinContent(i);
//		}	
	}
//	A_eff_L = sum_eff_wx_L / sum_eff_w_L;
//	A_eff_R = sum_eff_wx_R / sum_eff_w_R;
//	A_efferr_L = 1. / sqrt( sum_eff_w_L );
//	A_efferr_R = 1. / sqrt( sum_eff_w_R );
//	cout<<endl<<"A_eff_L["<<iBin<<"] = "<<A_eff_L<<" +- "<<A_efferr_L<<endl;
//	cout<<"A_eff_R["<<iBin<<"] = "<<A_eff_R<<" +- "<<A_efferr_R<<endl<<endl;
	fout->Write();
    fout->Close();
	
	TString f1_model_format_1 ;
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
//	if (iBin == 0 || iBin ==1 || iBin == 9 || iBin == 7) { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
		f1_model_format_2 = "( [0] + [1]*exp(-0.5* (((x-[2])/[3])**2)) ) * ( [4]+[5]*x+[6]*x**2+[7]*x**3+[8]*x**4+[9]*x**5+[10]*x**6 ) ";
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3 ) ";
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4 ) ";
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	}else { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0] + [1]*exp(-0.5* (((x-[2])/[3])**2)) ) * ( [4]+[5]*x+[6]*x**2+[7]*x**3+[8]*x**4+[9]*x**5+[10]*x**6 ) ";
	}
//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
//	Draw #events in acceptance
	h2_nacc_fine.Draw();
	canvas.Update();
////	canvas.Print(TString::Format("./plots/accXrecoEff_naccL_fine_bin%d.pdf",iBin));
//	Draw #events pass all cuts
	h2_nreco_fine.Draw();
	canvas.Update();
////	canvas.Print(TString::Format("./plots/accXrecoEff_nrecoL_fine_bin%d.pdf",iBin));

//	Draw FitResult for recoEfficiency
	const int nPar_r = 7;
	TF1 *f1_model_r = new TF1 ("f1_model_r", f1_model_format_1, -1., 1.);
	f1_model_r->SetParameter(0,0.);
	f1_model_r->SetParameter(1,0.01);
	f1_model_r->SetParameter(2,0.01);
	f1_model_r->SetParameter(3,0.01);
	f1_model_r->SetParameter(4,0.01);
	f1_model_r->SetParameter(5,0.01);
	f1_model_r->SetParameter(6,0.01);  // f1_model_format_1
		
	h2_reco_fine.Fit(f1_model_r,"R"); 
   
	Double_t matrix_r[7][7];
	gMinuit->mnemat(&matrix_r[0][0],7);
	cout<<"RECO ERROR MATRIX: "<<endl;
	cout<<"                   1             2            3            4            5            6            7       "<<endl;
	for (int i = 0; i<7; i++) {
		cout<<"RecoErr"<<i+1<<"  ";
		for (int j =0; j<7; j++) {
			cout<<matrix_r[i][j]<<" ";
		}
		cout<<endl;
	}
	
	h2_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetTitleOffset(1.15,"Y");
	h2_reco_fine.SetXTitle("cos#it{#theta_{l}}^{reco}");
	h2_reco_fine.SetYTitle("Reco-Efficiency / 0.1");
	h2_reco_fine.SetStats(0);
//	h2_reco_fine.SetMaximum(effUpperBound2[iBin]);
	h2_reco_fine.SetMaximum(h2_reco_fine.GetMaximum() * 1.2);
	h2_reco_fine.Draw("PE1");
//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	h2_reco_fine.Draw();
	f1_model_r->SetTitle("");
	f1_model_r->SetLineWidth(2);
	f1_model_r->SetLineColor(2);
	f1_model_r->Draw(" SAME ");
	
	
	t1->DrawLatex(.38,.90,TString::Format("Preliminary"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	t1->DrawLatex(.22,.90,TString::Format("Simulation"));
	t2->DrawLatex(.18,.82,TString::Format("#it{q^{2}}:"));
	t2->DrawLatex(.22,.82,Q2String[iBin]);
    if (iBin == 10) {
        t2->DrawLatex(.22,.77, " 14.18-22.00) GeV^{2}");
    }
	t2->DrawLatex(.84,.90,TString::Format("8TeV"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_recoL_fine_bin%d.pdf",iBin));
/*
//	Save fitting results
	double chi2Val_r=0;
	double arrPar_r[nPar_r], arrParErr_r[nPar_r];
	double recoPar1[nPar_r], recoPar2[nPar_r], recoPar3[nPar_r], recoPar4[nPar_r], recoPar5[nPar_r], recoPar6[nPar_r], recoPar7[nPar_r];
	for (int iPar = 0; iPar < nPar_r; iPar++) {
		arrPar_r[iPar]    = f1_model_r->GetParameter(iPar);
		arrParErr_r[iPar] = f1_model_r->GetParError(iPar);
		chi2Val_r         = f1_model_r->GetChisquare();
		recoPar1[iPar]   = matrix_r[0][iPar];
		recoPar2[iPar]   = matrix_r[1][iPar];
		recoPar3[iPar]   = matrix_r[2][iPar];
		recoPar4[iPar]   = matrix_r[3][iPar];
		recoPar5[iPar]   = matrix_r[4][iPar];
		recoPar6[iPar]   = matrix_r[5][iPar];
		recoPar7[iPar]   = matrix_r[6][iPar];
	}
	std::vector<double> output_r;
	for (int iPar = 0; iPar < nPar_r; iPar++){
		output_r.push_back(arrPar_r[iPar]);
		output_r.push_back(arrParErr_r[iPar]);
		printf("%18.15f,",arrPar_r[iPar]);
		if (iPar+1 >= nPar_r) printf("\n");
	}
	for (int i = 0; i < output_r.size(); i=i+2) {
		printf("%18.15f,",output_r[i+1]);
		if (i+2 >= output_r.size()) printf("\n");
	}
	writeParam(iBin,"reco",   arrPar_r,   nPar_r);  // f1_model_format_1
	writeParam(iBin,"recoErr",arrParErr_r,nPar_r);  // f1_model_format_1
	writeParam(iBin,"RecoErr1",       recoPar1,  nPar_r);  
	writeParam(iBin,"RecoErr2",       recoPar2,  nPar_r);  
	writeParam(iBin,"RecoErr3",       recoPar3,  nPar_r);  
	writeParam(iBin,"RecoErr4",       recoPar4,  nPar_r);  
	writeParam(iBin,"RecoErr5",       recoPar5,  nPar_r);  
	writeParam(iBin,"RecoErr6",       recoPar6,  nPar_r);  
	writeParam(iBin,"RecoErr7",       recoPar7,  nPar_r);  
//	Draw efficiency
//	h2_eff.SetMinimum(0.);
//	h2_eff.SetTitleOffset(1.3,"XY");
//	h2_eff.SetXTitle("genCosThetaL");
//	h2_eff.SetYTitle("Efficiency");
//	h2_eff.SetStats(0);
//	h2_eff.SetMaximum(effUpperBound[iBin]);  //03-11
//	h2_eff.Draw("PE1");
//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
///////////////////////////////////////////////////////////////////////////////////////////////

//	Draw FitResult for Total Efficiency
//	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
//	const int nPar = 11; // bin 0, 1, 9
	const int nPar = 7; // other bins
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"reco", 0));
	f1_model->SetParameter(4, readParam(iBin,"reco", 1));
	f1_model->SetParameter(5, readParam(iBin,"reco", 2));
	f1_model->SetParameter(6, readParam(iBin,"reco", 3));
//	f1_model->FixParameter(7, 0.0);
//	f1_model->FixParameter(8, 0.0);
//	f1_model->FixParameter(9, 0.0);
	if ( iBin == 2 ) {
		f1_model->SetParLimits(6, 0.0031679-0.001, 0.0031679+0.001); 
	}
	if ( iBin == 9) {
		f1_model->SetParLimits(6, 0.0006612-0.0005, 0.0006612+0.0005); 
	}
	if ( iBin == 1) {
		f1_model->SetParLimits(6, 0.0002844-0.0005, 0.0002844+0.0005); 
	}
	if ( iBin == 0) {
		f1_model->SetParLimits(6, -0.0013906-0.00090, -0.0013906-0.0004); 
	}
	if ( iBin == 0 || iBin == 1 || iBin == 9 ) {
//	if ( iBin == 12 ) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 0));
		f1_model->FixParameter(1, readParam(iBin,"acc", 1));
		f1_model->FixParameter(2, readParam(iBin,"acc", 2));
		f1_model->FixParameter(3, readParam(iBin,"acc", 3));
		f1_model->FixParameter(4, readParam(iBin,"reco", 0));
		f1_model->FixParameter(5, readParam(iBin,"reco", 1));
		f1_model->FixParameter(6, readParam(iBin,"reco", 2));
		f1_model->FixParameter(7, readParam(iBin,"reco", 3));
		f1_model->FixParameter(8, readParam(iBin,"reco", 4));
		f1_model->FixParameter(9, readParam(iBin,"reco", 5));
		f1_model->FixParameter(10, readParam(iBin,"reco", 6));		
////		f1_model->SetParError(0, readParam(iBin,"accErr", 1));
////		f1_model->SetParError(1, readParam(iBin,"accErr", 2));
////		f1_model->SetParError(2, readParam(iBin,"accErr", 3));
////		f1_model->SetParError(3, readParam(iBin,"recoErr", 0));
////		f1_model->SetParError(4, readParam(iBin,"recoErr", 1));
////		f1_model->SetParError(5, readParam(iBin,"recoErr", 2));
////		f1_model->SetParError(6, readParam(iBin,"recoErr", 3));
////		f1_model->SetParError(7, readParam(iBin,"recoErr", 4));
////		f1_model->SetParError(8, readParam(iBin,"recoErr", 5));
////		f1_model->SetParError(9, readParam(iBin,"recoErr", 6));
	} else {
//		f1_model->FixParameter(7, 0.0);
//		f1_model->FixParameter(8, 0.0);
//		f1_model->FixParameter(9, 0.0);
	}

	h2_eff_fine.SetStats(0);
//	h2_eff_fine.SetMinimum(-0.0001);
	h2_eff_fine.SetMinimum(0.);
	h2_eff_fine.SetTitleOffset(1.15,"Y");
	h2_eff_fine.SetXTitle("cos#theta_{l}^{reco}");
	h2_eff_fine.SetYTitle("Efficiency / 0.1");
//	h2_eff_fine.SetMaximum(effUpperBound[iBin]); 
	h2_eff_fine.SetMaximum(h2_eff_fine.GetMaximum() * 1.3); 
//	h2_eff_fine.Draw("TEXT");
	h2_eff_fine.Draw("PE1");
	
	h2_eff_fine.Fit(f1_model,"R"); //// 09-09
//	TFitResultPtr r = h2_eff_fine.Fit(f1_model,"WL S R");
//	r->Print();
   
	double matrix[nPar][nPar];  // bin 0,1,9
	gMinuit->mnemat(&matrix[0][0],nPar);
	cout<<"ERROR MATRIX: "<<endl;
	cout<<"            1             2            3            4            5            6            7         8         "<<endl;
	for (int i = 0; i<nPar; i++) {
		cout<<"EffErr"<<i+1<<"  ";
		for (int j =0; j<nPar; j++) {
			cout<<matrix[i][j]<<" ";
		}
		cout<<endl;
	}
	f1_model->SetTitle("");
//	f1_model->SetMaximum(effUpperBound[iBin]); 
	f1_model->SetLineWidth(2);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
	t1->DrawLatex(.23,.90,TString::Format("Simulation"));
	t1->DrawLatex(.18,.84,TString::Format("q^{2}"));
	t2->DrawLatex(.22,.84,Q2String[iBin]);
    if (iBin == 10) {
        t2->DrawLatex(.22,.79, "  14.18 - 22.00 GeV^{2}");
    }
	t2->DrawLatex(.84,.90,TString::Format("8TeV"));
  	tt->DrawLatex(.14,.90,TString::Format("CMS"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.pdf",iBin));
/////////////////////////////////////////////////////////////////////////////////////////////////
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int i = 0; i < nPar; i++) {arrPar[i]=0; arrParErr[i]=0;}
	double errPar1[nPar], errPar2[nPar], errPar3[nPar], errPar4[nPar], errPar5[nPar], errPar6[nPar], errPar7[nPar]; 
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
		errPar1[iPar]   = matrix[0][iPar];
		errPar2[iPar]   = matrix[1][iPar];
		errPar3[iPar]   = matrix[2][iPar];
		errPar4[iPar]   = matrix[3][iPar];
		errPar5[iPar]   = matrix[4][iPar];
		errPar6[iPar]   = matrix[5][iPar];
		errPar7[iPar]   = matrix[6][iPar];
//		if (iPar < 3) arrParErr[iPar] = readParam(iBin,"accErr", iPar+1);
//		else arrParErr[iPar] = readParam(iBin,"recoErr", iPar-3);
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
//	Draw compare
	TH1F h2_compFit("h2_compFit","",nLBins,-1.,1.);
	h2_compFit.SetXTitle("cos#theta_{l}^{gen}");
	TH1F h2_pullFit("h2_pullFit","",nLBins,-1.,1.);
	h2_pullFit.SetXTitle("cos#theta_{l}^{gen}");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_compFit.SetBinContent(i,f1_model->Eval(h2_eff_fine.GetXaxis()->GetBinCenter(i))/h2_eff_fine.GetBinContent(i));
			double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
			double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
			h2_pullFit.SetBinContent(i,(f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}else{
			h2_compFit.SetBinContent(i,0.);
			h2_pullFit.SetBinContent(i,0.);
		}
	}
	h2_compFit.SetMinimum(0.);
	h2_compFit.SetStats(0);
	h2_compFit.Draw("PE1");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.60,0.91,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_compFit_bin%d.pdf",iBin));
	
	h2_pullFit.SetStats(0);
	h2_pullFit.Draw(" HIST TEXT");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.50,0.91,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_pullFit_bin%d.pdf",iBin));
	
//	Draw significance of deviation
	TH1F h2_pull("PULL - Deviation/Error","",15,-3.,3.);
	h2_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_pull.Fill((f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}
	}
	h2_pull.Draw("HIST E1");
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_sigma_bin%d.pdf",iBin));
		
	delete latex;
	
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"accXrecoEff",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEffErr",arrParErr,nPar);  // f1_model_format_2
	writeParam(iBin,"EffErr1",       errPar1,  nPar);  
	writeParam(iBin,"EffErr2",       errPar2,  nPar);  
	writeParam(iBin,"EffErr3",       errPar3,  nPar);  
	writeParam(iBin,"EffErr4",       errPar4,  nPar);  
	writeParam(iBin,"EffErr5",       errPar5,  nPar);  
	writeParam(iBin,"EffErr6",       errPar6,  nPar);  
	writeParam(iBin,"EffErr7",       errPar7,  nPar); 
//	return output.c_str();	
*/
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEffSYS(int iBin)
{//{{{
	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	
	TFile f_eff(TString::Format("./RootFiles/Efficiency_%d.root",iBin));
	TH1F *h2_eff_Fine = (TH1F*)f_eff.Get("h2_eff_fine"); 
	
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) ";
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) + ( [3]*exp(-0.5* (((x-[4])/[5])**2)) ) ";
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
		f1_model_format_2 = "( [0] + [1]*exp(-0.5* (((x-[2])/[3])**2)) ) * ( [4]+[5]*x+[6]*x**2+[7]*x**3+[8]*x**4+[9]*x**5+[10]*x**6 ) ";
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
	}else { 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) + ( [3]*exp(-0.5* (((x-[4])/[5])**2)) ) ";
		f1_model_format_2 = "( [0] + [1]*exp(-0.5* (((x-[2])/[3])**2)) ) * ( [4]+[5]*x+[6]*x**2+[7]*x**3+[8]*x**4+[9]*x**5+[10]*x**6 ) ";
	}
//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
	const int nPar = 11; // other bins
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
	f1_model->SetParameter(0, readParam(iBin,"acc", 0));
	f1_model->SetParameter(1, readParam(iBin,"acc", 1));
	f1_model->SetParameter(2, readParam(iBin,"acc", 2));
	f1_model->SetParameter(3, readParam(iBin,"acc", 3));
	f1_model->SetParameter(4, readParam(iBin,"reco", 0));
	f1_model->SetParameter(5, readParam(iBin,"reco", 1));
	f1_model->SetParameter(6, readParam(iBin,"reco", 2));
	f1_model->SetParameter(7, readParam(iBin,"reco", 3));
	f1_model->SetParameter(8, readParam(iBin,"reco", 4));
//	if ( iBin == 6 || iBin == 7 || iBin == 10) {
	if ( iBin == 6 || iBin == 7) {
        f1_model->SetParameter(9, readParam(iBin,"reco", 5));      
        f1_model->SetParameter(10, readParam(iBin,"reco", 6));
		f1_model->SetParLimits(1, readParam(iBin,"acc", 1)-readParam(iBin,"accErr", 1), readParam(iBin,"acc", 1)+readParam(iBin,"accErr", 1));
		f1_model->SetParLimits(2, readParam(iBin,"acc", 2)-0.01*readParam(iBin,"accErr", 2), readParam(iBin,"acc", 2)+0.01*readParam(iBin,"accErr", 2));
		f1_model->SetParLimits(3, readParam(iBin,"acc", 3)-0.5*readParam(iBin,"accErr", 3), readParam(iBin,"acc", 3)+0.5*readParam(iBin,"accErr", 3));
	} else if ( iBin == 0 ) {
		f1_model->SetParLimits(1, readParam(iBin,"acc", 1)-readParam(iBin,"accErr", 1), readParam(iBin,"acc", 1)+readParam(iBin,"accErr", 1));
		f1_model->SetParLimits(2, readParam(iBin,"acc", 2)-0.01*readParam(iBin,"accErr", 2), readParam(iBin,"acc", 2)+0.01*readParam(iBin,"accErr", 2));
		f1_model->SetParLimits(3, readParam(iBin,"acc", 3)-0.1*readParam(iBin,"accErr", 3), readParam(iBin,"acc", 3)+0.1*readParam(iBin,"accErr", 3));
	} else if ( iBin == 1 || iBin == 9) {
      f1_model->FixParameter(7, readParam(iBin,"reco", 3));
		f1_model->FixParameter(9, 0.0);
		f1_model->FixParameter(10, 0.0);
		f1_model->SetParLimits(1, readParam(iBin,"acc", 1)-readParam(iBin,"accErr", 1), readParam(iBin,"acc", 1)+readParam(iBin,"accErr", 1));
		f1_model->SetParLimits(2, readParam(iBin,"acc", 2)-0.01*readParam(iBin,"accErr", 2), readParam(iBin,"acc", 2)+0.01*readParam(iBin,"accErr", 2));
		f1_model->SetParLimits(3, readParam(iBin,"acc", 3)-0.1*readParam(iBin,"accErr", 3), readParam(iBin,"acc", 3)+0.1*readParam(iBin,"accErr", 3));
	} else if ( iBin == 10 ) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 0));
		f1_model->FixParameter(1, readParam(iBin,"acc", 1));
		f1_model->FixParameter(2, readParam(iBin,"acc", 2));
		f1_model->FixParameter(3, readParam(iBin,"acc", 3));
		f1_model->FixParameter(4, readParam(iBin,"reco", 0));
		f1_model->FixParameter(5, readParam(iBin,"reco", 1));
		f1_model->FixParameter(6, readParam(iBin,"reco", 2));
		f1_model->FixParameter(7, readParam(iBin,"reco", 3));
		f1_model->FixParameter(8, readParam(iBin,"reco", 4));
		f1_model->FixParameter(9, readParam(iBin,"reco", 5));
		f1_model->FixParameter(10, readParam(iBin,"reco", 6));		
	} else {
		f1_model->FixParameter(9, 0.0);
		f1_model->FixParameter(10, 0.0);
		f1_model->SetParLimits(1, readParam(iBin,"acc", 1)-readParam(iBin,"accErr", 1), readParam(iBin,"acc", 1)+readParam(iBin,"accErr", 1));
		f1_model->SetParLimits(2, readParam(iBin,"acc", 2)-0.01*readParam(iBin,"accErr", 2), readParam(iBin,"acc", 2)+0.01*readParam(iBin,"accErr", 2));
		f1_model->SetParLimits(3, readParam(iBin,"acc", 3)-0.5*readParam(iBin,"accErr", 3), readParam(iBin,"acc", 3)+0.5*readParam(iBin,"accErr", 3));
	}

	h2_eff_Fine->SetStats(0);
//	h2_eff_Fine->SetMinimum(-0.0001);
	h2_eff_Fine->SetMinimum(0.);
	h2_eff_Fine->SetTitleOffset(1.15,"Y");
	h2_eff_Fine->SetXTitle("cos#theta_{l}^{reco}");
	h2_eff_Fine->SetYTitle("Efficiency / 0.1");
	h2_eff_Fine->SetMaximum(h2_eff_Fine->GetMaximum() * 1.2); 
	h2_eff_Fine->Draw("PE1");
	
	h2_eff_Fine->Fit(f1_model,"R"); //// 09-09
	f1_model->SetTitle("");
	f1_model->SetLineWidth(2);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.22,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.pdf",iBin));
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int i = 0; i < nPar; i++) {arrPar[i]=0; arrParErr[i]=0;}
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
//	return output;
	writeParam(iBin,"accXrecoEff",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEffErr",arrParErr,nPar);  // f1_model_format_2

	delete latex;
	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEffFIT(int iBin)
{//{{{
	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	
	TFile f_eff(TString::Format("./RootFiles/Efficiency_%d.root",iBin));
	TH1F *h2_eff_Fine = (TH1F*)f_eff.Get("h2_eff_fine"); 
	
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
		f1_model_format_2 = "100 * ( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
////		f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	}else { 
		f1_model_format_2 = "100* ([0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6)"; 
////		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	}
//	Draw
	TCanvas canvas("canvas");
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	tt->SetTextSize(0.09);
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->SetTextSize(0.08);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t2->SetTextSize(0.09);
//	const int nPar = 7; // other bins
	const int nPar = 10; // other bins
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
//	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -0.9, 0.9);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"reco", 0));
	f1_model->SetParameter(4, readParam(iBin,"reco", 1));
	f1_model->SetParameter(5, readParam(iBin,"reco", 2));
//	f1_model->SetParameter(6, readParam(iBin,"reco", 3));
	if ( iBin == 2) {
		f1_model->SetParLimits(6, 0.0031679-0.001, 0.0031679+0.001); 
	}
/*	if ( iBin == 9) {
		f1_model->SetParLimits(6, 0.0006612-0.0005, 0.0006612+0.0005); 
	}
	if ( iBin == 1) {
		f1_model->SetParLimits(6, 0.0002844-0.0005, 0.0002844+0.0005); 
	}
	if ( iBin == 0) {
      //f1_model->SetParameter(7, readParam(iBin,"reco", 4));
//		f1_model->SetParLimits(6, -0.0013906-0.00090, -0.0013906-0.0004); 
//	}
*/	if ( iBin == 0 || iBin == 1 || iBin == 9) {
//	if ( iBin == 20) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 1));
		f1_model->FixParameter(1, readParam(iBin,"acc", 2));
		f1_model->FixParameter(2, readParam(iBin,"acc", 3));
		f1_model->FixParameter(3, readParam(iBin,"reco", 0));
		f1_model->FixParameter(4, readParam(iBin,"reco", 1));
		f1_model->FixParameter(5, readParam(iBin,"reco", 2));
		f1_model->FixParameter(6, readParam(iBin,"reco", 3));
		f1_model->FixParameter(7, readParam(iBin,"reco", 4));
		f1_model->FixParameter(8, readParam(iBin,"reco", 5));
		f1_model->FixParameter(9, readParam(iBin,"reco", 6));		
	} else {
		f1_model->FixParameter(7, 0.0);
		f1_model->FixParameter(8, 0.0);
		f1_model->FixParameter(9, 0.0);
	}

	h2_eff_Fine->SetStats(0);
//	h2_eff_Fine->SetMinimum(-0.0001);
	h2_eff_Fine->SetTitleOffset(0.57,"X");
	h2_eff_Fine->SetTitleOffset(0.64,"Y");
////    h2_eff_Fine->GetYaxis()->CenterTitle();
////    TGaxis::SetMaxDigits(2);
////    h2_eff_Fine->GetXaxis()->CenterTitle();
	h2_eff_Fine->SetTitleSize(0.10, "Y");
	h2_eff_Fine->SetTitleSize(0.10, "X");
	h2_eff_Fine->SetXTitle("cos#theta_{#it{l}} ");
	h2_eff_Fine->SetYTitle("Efficiency (%) ");
	h2_eff_Fine->GetYaxis()->SetNdivisions(505);
	h2_eff_Fine->GetXaxis()->SetNdivisions(505);
	h2_eff_Fine->SetTickLength(0.06, "XY");
    h2_eff_Fine->SetLabelFont(22, "XY");
    h2_eff_Fine->SetLabelOffset(0.007, "XY");
    h2_eff_Fine->SetLabelSize(0.08, "XY");
////	h2_eff_Fine->SetMaximum(h2_eff_Fine->GetMaximum() * 1.6); 
	h2_eff_Fine->Scale(100);
	h2_eff_Fine->SetMaximum(0.35); 
	h2_eff_Fine->SetMinimum(0.);
    h2_eff_Fine->SetMarkerSize(1.7);
    h2_eff_Fine->SetLineWidth(3);
////	h2_eff_Fine->Draw("PE1 ");
	
	h2_eff_Fine->Fit(f1_model,"R"); //// 09-09
	
	f1_model->SetTitle("");
	f1_model->SetLineWidth(3);
	f1_model->SetLineColor(2);
	f1_model->Draw("SAME");
	if (iBin==0) h2_eff_Fine->SetBinError(1,0);
    gStyle->SetErrorX(0.000);
	h2_eff_Fine->Draw("P SAME");
	
////	t1->DrawLatex(.20,.90,TString::Format("Preliminary"));
  	tt->DrawLatex(.22,.825,TString::Format("CMS"));
	t1->DrawLatex(.37,.825,TString::Format("Simulation"));
	t2->DrawLatex(.22,.70,Q2String[iBin]);
////	tt->DrawLatex(.52,.24,FIgString[iBin]);
////	t2->DrawLatex(.18,.82,TString::Format("#it{q^{2}}:"));
////	t2->DrawLatex(.22,.82,Q2String[iBin]);
////    if (iBin == 10) {
////        t2->DrawLatex(.22,.77, " 14.18-22.00) GeV^{2}");
////    }
////	t2->DrawLatex(.80,.925,TString::Format("8 TeV"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d.pdf",iBin));
	
//	Save Fitting results
	TLatex *latex = new TLatex();
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int i = 0; i < nPar; i++) {arrPar[i]=0; arrParErr[i]=0;}
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
//	return output;
////	writeParam(iBin,"accXrecoEff",   arrPar,   nPar);  // f1_model_format_2
////	writeParam(iBin,"accXrecoEffErr",arrParErr,nPar);  // f1_model_format_2
/////////////////////////////////////////////////////////////////////////////////////////////////
	
	delete latex;
	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEffGEN(int iBin)
//std::vector<double> accXrecoEff(int iBin)
{//{{{
	setTDRStyle();
//	TH1::SetDefaultSumw2();
//	ROOT::Math::Minimizer::SetDefaultMaxFunctionCalls(10000);
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	TLorentzVector B_4vec;
//	double Bctau = 0;
	int    Triggers = 0;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double gCosThetaL = 0;
	double CosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;

//	double Nacc = 0;
//	double Nacc_gen = 0;
//	double Nreco = 0;
//	double Nreco_gen = 0;
	
	ch->SetBranchStatus("Triggers"      , 1);
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("Triggers"     , &Triggers);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Load acceptance
	TFile f_acc("./RootFiles/acceptance_8TeV.root");
	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
//	TH1F *h1_ngen_fine = (TH1F*)f_acc.Get(TString::Format("h1_ngen_fine_bin%d",iBin));
	
//	Fill histograms
	int nLBins = 20;// The same value as h_acc
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,-1,1); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,-1,1);
	TH1F h2_nRECO_fine("h2_ngen_fine","h2_nRECO_fine",nLBins,-1,1);
	
	h2_nacc_fine.SetStats(0);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.1");
	h2_nreco_fine.SetStats(0);
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.1");
	
	TH1F h2_D1_fine("h2_D1_fine","",1000,-0.1,0.1);
	double fffff = 0;
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
	//	if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 2.8 && gmuppt > 2.8 ){
			h2_nacc_fine.Fill(gCosThetaL);
		}
	//	if (BMass != -999 && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)) ){
		if (BMass != -999 && Triggers == 1){    ////////////////////////   12-10 N.A.
			h2_nreco_fine.Fill(gCosThetaL);
			h2_nRECO_fine.Fill(CosThetaL);
			
			fffff = CosThetaL - gCosThetaL;
			h2_D1_fine.Fill(fffff);
		}
	}
	
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,1.);
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	}
	TH1F h2_RECO_fine("h2_RECO_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0) {
			h2_RECO_fine.SetBinContent(i,0.);
			h2_RECO_fine.SetBinError(i,1.);
		}else{
			h2_RECO_fine.SetBinContent(i,h2_nRECO_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_RECO_fine.SetBinError(i,sqrt(h2_RECO_fine.GetBinContent(i)*(1-h2_RECO_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
		}
	}
	
	TH1F h2_eff_fine("h2_eff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_eff_fine.SetBinContent(i,0.);
			h2_eff_fine.SetBinError(i,1.);
		}else{
			h2_eff_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(pow(h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i),2)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
			printf("INFO: Efficiency_fine(%d)=%f +- %f.\n",i,h2_eff_fine.GetBinContent(i),h2_eff_fine.GetBinError(i));
		}
	}

	TH1F h2_EFF_fine("h2_EFF_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			h2_EFF_fine.SetBinContent(i,0.);
			h2_EFF_fine.SetBinError(i,1.);
		}else{
			h2_EFF_fine.SetBinContent(i,h2_nRECO_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_EFF_fine.SetBinError(i,h2_EFF_fine.GetBinContent(i)*sqrt(pow(h2_RECO_fine.GetBinError(i)/h2_RECO_fine.GetBinContent(i),2)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
		}
	}
	TH1F h2_diff_fine("h2_diff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			h2_diff_fine.SetBinContent(i,0.);
		}else{
			h2_diff_fine.SetBinContent(i, (h2_EFF_fine.GetBinContent(i) - h2_eff_fine.GetBinContent(i)) / h2_EFF_fine.GetBinContent(i) );
		}
	}
	
	TString f1_model_format_1 ;
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
		f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	}else { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
	}
//	Draw
	TCanvas canvas("canvas");
	canvas.SetGrid();
	TLatex *latex = new TLatex();
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
//	Draw #events in acceptance
	h2_nacc_fine.Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_naccL_fine_bin%d.pdf",iBin));
//	Draw #events pass all cuts
	h2_nreco_fine.Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_nrecoL_fine_bin%d.pdf",iBin));
	
	gStyle->SetOptStat();
	h2_D1_fine.SetFillColor(1);
	h2_D1_fine.Draw(" b ");
	h2_D1_fine.SetXTitle("cos#theta^{reco}_{l} - cos#theta^{gen}_{l}");
	h2_D1_fine.SetYTitle("Events/0.0002");
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_D1_fine_bin%d.pdf",iBin));
	
	double UpBound[11]  = {  0.15, 0.050, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030};
	double DoBound[11]  = { -0.15,-0.050,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030};
	h2_diff_fine.SetStats(0);
	h2_diff_fine.SetMaximum(UpBound[iBin]); 
	h2_diff_fine.SetMinimum(DoBound[iBin]); 
	h2_diff_fine.SetFillColor(7);
	h2_diff_fine.Draw(" b TEXT");
	h2_diff_fine.SetXTitle("cos#theta_{l}");
	h2_diff_fine.SetYTitle("Efficiency Resolution/0.1");
	TLatex *tdiff = new TLatex();
	tdiff->DrawLatexNDC(0.51,0.91,TString::Format("(E_{reco} - E_{gen}) / E_{reco} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_diff_fine_bin%d.pdf",iBin));

//	Draw FitResult for recoEfficiency
	const int nPar_r = 7;
	TF1 *f1_model_r = new TF1 ("f1_model_r", f1_model_format_1, -1., 1.);
	f1_model_r->SetParameter(0,0.);
	f1_model_r->SetParameter(1,0.01);
	f1_model_r->SetParameter(2,0.01);
	f1_model_r->SetParameter(3,0.01);
	f1_model_r->SetParameter(4,0.01);
	f1_model_r->SetParameter(5,0.01);
	f1_model_r->SetParameter(6,0.01);  // f1_model_format_1
		
	h2_reco_fine.Fit(f1_model_r,"R"); 
	
	h2_reco_fine.SetMinimum(0.);
	h2_reco_fine.SetTitleOffset(1.15,"Y");
	h2_reco_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_reco_fine.SetYTitle("Reco-Efficiency / 0.1");
	h2_reco_fine.SetStats(0);
	h2_reco_fine.SetMaximum(h2_reco_fine.GetMaximum() * 1.2); 
	h2_reco_fine.Draw("PE1");
//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	h2_reco_fine.Draw();
	f1_model_r->SetTitle("");
	f1_model_r->SetLineWidth(2);
	f1_model_r->SetLineColor(2);
	f1_model_r->Draw(" SAME ");
	
    paveText->Draw();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_recoL_fine_bin%d.pdf",iBin));

//	Save fitting results
	double chi2Val_r=0;
	double arrPar_r[nPar_r], arrParErr_r[nPar_r];
	for (int iPar = 0; iPar < nPar_r; iPar++) {
		arrPar_r[iPar]    = f1_model_r->GetParameter(iPar);
		arrParErr_r[iPar] = f1_model_r->GetParError(iPar);
		chi2Val_r         = f1_model_r->GetChisquare();
	}
	std::vector<double> output_r;
	for (int iPar = 0; iPar < nPar_r; iPar++){
		output_r.push_back(arrPar_r[iPar]);
		output_r.push_back(arrParErr_r[iPar]);
		printf("%18.15f,",arrPar_r[iPar]);
		if (iPar+1 >= nPar_r) printf("\n");
	}
	for (int i = 0; i < output_r.size(); i=i+2) {
		printf("%18.15f,",output_r[i+1]);
		if (i+2 >= output_r.size()) printf("\n");
	}
	writeParam(iBin,"GENreco",   arrPar_r,   nPar_r);  // f1_model_format_1
	writeParam(iBin,"GENrecoErr",arrParErr_r,nPar_r);  // f1_model_format_1
///////////////////////////////////////////////////////////////////////////////////////////////

//	Draw FitResult for Total Efficiency
//	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
	const int nPar = 10;
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"reco", 0));
	f1_model->SetParameter(4, readParam(iBin,"reco", 1));
	f1_model->SetParameter(5, readParam(iBin,"reco", 2));
	f1_model->SetParameter(6, readParam(iBin,"reco", 3));
//	f1_model->FixParameter(7, 0.0);
//	f1_model->FixParameter(8, 0.0);
//	f1_model->FixParameter(9, 0.0);
	if ( iBin == 0 || iBin == 1 || iBin == 9 ) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 1));
		f1_model->FixParameter(1, readParam(iBin,"acc", 2));
		f1_model->FixParameter(2, readParam(iBin,"acc", 3));
		f1_model->FixParameter(3, readParam(iBin,"reco", 0));
		f1_model->FixParameter(4, readParam(iBin,"reco", 1));
		f1_model->FixParameter(5, readParam(iBin,"reco", 2));
		f1_model->FixParameter(6, readParam(iBin,"reco", 3));
		f1_model->FixParameter(7, readParam(iBin,"reco", 4));
		f1_model->FixParameter(8, readParam(iBin,"reco", 5));
		f1_model->FixParameter(9, readParam(iBin,"reco", 6));		
/*		f1_model->SetParError(0, readParam(iBin,"accErr", 1));
		f1_model->SetParError(1, readParam(iBin,"accErr", 2));
		f1_model->SetParError(2, readParam(iBin,"accErr", 3));
		f1_model->SetParError(3, readParam(iBin,"recoErr", 0));
		f1_model->SetParError(4, readParam(iBin,"recoErr", 1));
		f1_model->SetParError(5, readParam(iBin,"recoErr", 2));
		f1_model->SetParError(6, readParam(iBin,"recoErr", 3));
		f1_model->SetParError(7, readParam(iBin,"recoErr", 4));
		f1_model->SetParError(8, readParam(iBin,"recoErr", 5));
		f1_model->SetParError(9, readParam(iBin,"recoErr", 6));
*/	} else {
		f1_model->FixParameter(7, 0.0);
		f1_model->FixParameter(8, 0.0);
		f1_model->FixParameter(9, 0.0);
	}

	h2_eff_fine.SetStats(0);
//	h2_eff_fine.SetMinimum(-0.00005);
	h2_eff_fine.SetMinimum(0.);
	h2_eff_fine.SetTitleOffset(1.15,"Y");
	h2_eff_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_eff_fine.SetYTitle("Efficiency / 0.1");
	h2_eff_fine.SetMaximum(h2_eff_fine.GetMaximum() * 1.2); 
//	h2_eff_fine.Draw("TEXT");
	h2_eff_fine.Draw("PE1");
	
	h2_eff_fine.Fit(f1_model,"R"); //// 09-09
//	TFitResultPtr r = h2_eff_fine.Fit(f1_model,"WL S R");
//	r->Print();
	
	f1_model->SetTitle("");
	f1_model->SetLineWidth(2);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
    paveText->Draw();
	t1->DrawLatex(.22,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_Eff_fine_bin%d.pdf",iBin));
/////////////////////////////////////////////////////////////////////////////////////////////////
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
//	Draw compare
	TH1F h2_compFit("h2_compFit","",nLBins,-1.,1.);
	h2_compFit.SetXTitle("cos#theta_{l}^{gen}");
	TH1F h2_pullFit("h2_pullFit","",nLBins,-1.,1.);
	h2_pullFit.SetXTitle("cos#theta_{l}^{gen}");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_compFit.SetBinContent(i,f1_model->Eval(h2_eff_fine.GetXaxis()->GetBinCenter(i))/h2_eff_fine.GetBinContent(i));
			double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
			double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
			h2_pullFit.SetBinContent(i,(f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}else{
			h2_compFit.SetBinContent(i,0.);
			h2_pullFit.SetBinContent(i,0.);
		}
	}
	h2_compFit.SetMinimum(0.);
	h2_compFit.SetStats(0);
	h2_compFit.Draw("PE1");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.60,0.91,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_compFit_bin%d.pdf",iBin));
	
	h2_pullFit.SetStats(0);
	h2_pullFit.Draw(" HIST TEXT");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.50,0.91,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_pullFit_bin%d.pdf",iBin));
	
//	Draw significance of deviation
	TH1F h2_pull("PULL - Deviation/Error","",15,-3.,3.);
	h2_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		double _xlo = h2_eff_fine.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h2_eff_fine.GetXaxis()->GetBinUpEdge(i);
		if (h2_eff_fine.GetBinContent(i) != 0){
			h2_pull.Fill((f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_eff_fine.GetBinContent(i))/h2_eff_fine.GetBinError(i));
		}
	}
	h2_pull.Draw("HIST E1");
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff-GEN_sigma_bin%d.pdf",iBin));
		
	delete latex;
	
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"accXrecoEff-GEN",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEff-GENErr",arrParErr,nPar);  // f1_model_format_2
//	return output.c_str();	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEffDimuon(int iBin)
//std::vector<double> accXrecoEff(int iBin)
{//{{{
	setTDRStyle();
//	TH1::SetDefaultSumw2();
//	ROOT::Math::Minimizer::SetDefaultMaxFunctionCalls(10000);
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
//	TLorentzVector B_4vec;
//	double Bctau = 0;
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double gCosThetaL = 0;
	double CosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;
//	double Nacc = 0;
//	double Nacc_gen = 0;
//	double Nreco = 0;
//	double Nreco_gen = 0;
	
	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Load acceptance
	TFile f_acc("./RootFiles/acceptance_8TeV.root");
	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
//	TH1F *h1_ngen_fine = (TH1F*)f_acc.Get(TString::Format("h1_ngen_fine_bin%d",iBin));
	
//	Fill histograms
	int nLBins = 20;// The same value as h_acc
	TH1F h2_nacc_fine("h2_nacc_fine" ,"h2_nacc_fine" ,nLBins,-1,1); 
	TH1F h2_nreco_fine("h2_nreco_fine","h2_nreco_fine",nLBins,-1,1);
	TH1F h2_nACC_fine("h2_nACC_fine" ,"h2_nACC_fine" ,nLBins,-1,1); 
	TH1F h2_nRECO_fine("h2_nRECO_fine","h2_nRECO_fine",nLBins,-1,1);
	
	h2_nacc_fine.SetStats(0);
	h2_nacc_fine.SetMinimum(0.);
	h2_nacc_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_nacc_fine.SetYTitle("Acceptated Events/0.1");
	h2_nreco_fine.SetStats(0);
	h2_nreco_fine.SetMinimum(0.);
	h2_nreco_fine.SetXTitle("cos#theta_{l}^{gen}");
	h2_nreco_fine.SetYTitle("Reconstructed Events/0.1");
	
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
			h2_nacc_fine.Fill(gCosThetaL);
			h2_nACC_fine.Fill(gCosThetaL);
		}
		if (BMass != -999 ){    ////////////////////////   12-10 N.A.
			h2_nreco_fine.Fill(CosThetaL);
		}
	}
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (Q2 > Q2rangeup[iBin] || Q2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (Q2 < Q2rangeup[3] && Q2 > Q2rangedn[3]) continue;
				if (Q2 < Q2rangeup[5] && Q2 > Q2rangedn[5]) continue;
			}
//		if ( fabs(gmumeta) < 2.5 && fabs(gmupeta) < 2.5 && gmumpt > 3.5 && gmuppt > 3.5 ){
//			h2_nACC_fine.Fill(gCosThetaL);
//		}
		if (BMass != -999 ){    ////////////////////////   12-10 N.A.
			h2_nRECO_fine.Fill(CosThetaL);
		}
	}
	
	TH1F h2_reco_fine("h2_reco_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_reco_fine.SetBinContent(i,0.);
			h2_reco_fine.SetBinError(i,1.);
		}else{
			h2_reco_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i));
			h2_reco_fine.SetBinError(i,sqrt(h2_reco_fine.GetBinContent(i)*(1-h2_reco_fine.GetBinContent(i))/h2_nacc_fine.GetBinContent(i)));
			printf("INFO: recoEfficiency_fine(%d)=%f +- %f.\n",i,h2_reco_fine.GetBinContent(i),h2_reco_fine.GetBinError(i));
		}
	}
	TH1F h2_RECO_fine("h2_RECO_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nACC_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0) {
			h2_RECO_fine.SetBinContent(i,0.);
			h2_RECO_fine.SetBinError(i,1.);
		}else{
			h2_RECO_fine.SetBinContent(i,h2_nRECO_fine.GetBinContent(i)/h2_nACC_fine.GetBinContent(i));
			h2_RECO_fine.SetBinError(i,sqrt(h2_RECO_fine.GetBinContent(i)*(1-h2_RECO_fine.GetBinContent(i))/h2_nACC_fine.GetBinContent(i)));
			printf("INFO: RECOEfficiency_fine(%d)=%f +- %f.\n",i,h2_RECO_fine.GetBinContent(i),h2_RECO_fine.GetBinError(i));
		}
	}
	
	TH1F h2_eff_fine("h2_eff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
			h2_eff_fine.SetBinContent(i,0.);
			h2_eff_fine.SetBinError(i,1.);
		}else{
			h2_eff_fine.SetBinContent(i,h2_nreco_fine.GetBinContent(i)/h2_nacc_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_eff_fine.SetBinError(i,h2_eff_fine.GetBinContent(i)*sqrt(pow(h2_reco_fine.GetBinError(i)/h2_reco_fine.GetBinContent(i),2)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
			printf("INFO: Efficiency_fine(%d)=%f +- %f.\n",i,h2_eff_fine.GetBinContent(i),h2_eff_fine.GetBinError(i));
		}
	}

	TH1F h2_EFF_fine("h2_EFF_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nACC_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			h2_EFF_fine.SetBinContent(i,0.);
			h2_EFF_fine.SetBinError(i,1.);
		}else{
			h2_EFF_fine.SetBinContent(i,h2_nRECO_fine.GetBinContent(i)/h2_nACC_fine.GetBinContent(i) * h1_acc_fine->GetBinContent(i));
			h2_EFF_fine.SetBinError(i,h2_EFF_fine.GetBinContent(i)*sqrt(pow(h2_RECO_fine.GetBinError(i)/h2_RECO_fine.GetBinContent(i),2)+pow(h1_acc_fine->GetBinError(i)/h1_acc_fine->GetBinContent(i),2)));
			printf("INFO: EFFICIENCY_fine(%d)=%f +- %f.\n",i,h2_EFF_fine.GetBinContent(i),h2_EFF_fine.GetBinError(i));
		}
	}
	TH1F h2_diff_fine("h2_diff_fine","",nLBins,-1,1);
	for (int i = 1; i <= nLBins; i++) {
		if (h2_nacc_fine.GetBinContent(i) == 0 || h2_nACC_fine.GetBinContent(i) == 0 || h2_nRECO_fine.GetBinContent(i) == 0 || h2_nreco_fine.GetBinContent(i) == 0 || h1_acc_fine->GetBinContent(i) == 0 ) {
			h2_diff_fine.SetBinContent(i,0.);
		}else{
			h2_diff_fine.SetBinContent(i, (h2_EFF_fine.GetBinContent(i) - h2_eff_fine.GetBinContent(i)) / h2_EFF_fine.GetBinContent(i) );
		}
	}
	
	TString f1_model_format_1 ;
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4 ) ";
		f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
	}else { 
		f1_model_format_1 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6"; 
	//	f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
	//	f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]"; 
	}
//	Draw
	TCanvas canvas("canvas");
	canvas.SetGrid();
	TLatex *latex = new TLatex();
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
//	Draw #events in acceptance
	h2_nACC_fine.Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_nACCL_fine_bin%d.pdf",iBin));
//	Draw #events pass all cuts
	h2_nRECO_fine.Draw();
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_nRECOL_fine_bin%d.pdf",iBin));
	
	double UpBound[11]  = {  0.10, 0.050, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030};
	double DoBound[11]  = { -0.03,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030,-0.030};
	h2_diff_fine.SetStats(0);
	h2_diff_fine.SetMaximum(UpBound[iBin]); 
	h2_diff_fine.SetMinimum(DoBound[iBin]); 
	h2_diff_fine.SetFillColor(7);
	h2_diff_fine.Draw(" b TEXT");
	h2_diff_fine.SetXTitle("cos#theta^{reco}_{l}");
	h2_diff_fine.SetYTitle("Efficiency Resolution/0.1");
	TLatex *tdiff = new TLatex();
	tdiff->DrawLatexNDC(0.51,0.92,TString::Format("(E_{reco} - E_{gen}) / E_{reco} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_diff_fine_bin%d.pdf",iBin));

//	Draw FitResult for recoEfficiency
	const int nPar_r = 7;
	TF1 *f1_model_r = new TF1 ("f1_model_r", f1_model_format_1, -1., 1.);
	f1_model_r->SetParameter(0,0.);
	f1_model_r->SetParameter(1,0.01);
	f1_model_r->SetParameter(2,0.01);
	f1_model_r->SetParameter(3,0.01);
	f1_model_r->SetParameter(4,0.01);
	f1_model_r->SetParameter(5,0.01);
	f1_model_r->SetParameter(6,0.01);  // f1_model_format_1
		
	h2_RECO_fine.Fit(f1_model_r,"R"); 
	
	h2_RECO_fine.SetMinimum(0.);
	h2_RECO_fine.SetTitleOffset(1.15,"Y");
	h2_RECO_fine.SetXTitle("cos#theta_{l}^{reco}");
	h2_RECO_fine.SetYTitle("Reco-Efficiency / 0.1");
	h2_RECO_fine.SetStats(0);
	h2_RECO_fine.SetMaximum(h2_RECO_fine.GetMaximum() * 1.2); 
	h2_RECO_fine.Draw("PE1");
//	latex->DrawLatexNDC(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
	h2_RECO_fine.Draw();
	f1_model_r->SetTitle("");
	f1_model_r->SetLineWidth(2);
	f1_model_r->SetLineColor(2);
	f1_model_r->Draw(" SAME ");
	
    paveText->Draw();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_RECOL_fine_bin%d.pdf",iBin));

//	Save fitting results
	double chi2Val_r=0;
	double arrPar_r[nPar_r], arrParErr_r[nPar_r];
	for (int iPar = 0; iPar < nPar_r; iPar++) {
		arrPar_r[iPar]    = f1_model_r->GetParameter(iPar);
		arrParErr_r[iPar] = f1_model_r->GetParError(iPar);
		chi2Val_r         = f1_model_r->GetChisquare();
	}
	std::vector<double> output_r;
	for (int iPar = 0; iPar < nPar_r; iPar++){
		output_r.push_back(arrPar_r[iPar]);
		output_r.push_back(arrParErr_r[iPar]);
		printf("%18.15f,",arrPar_r[iPar]);
		if (iPar+1 >= nPar_r) printf("\n");
	}
	for (int i = 0; i < output_r.size(); i=i+2) {
		printf("%18.15f,",output_r[i+1]);
		if (i+2 >= output_r.size()) printf("\n");
	}
	writeParam(iBin,"Dimuonreco",   arrPar_r,   nPar_r);  // f1_model_format_1
	writeParam(iBin,"DimuonrecoErr",arrParErr_r,nPar_r);  // f1_model_format_1
///////////////////////////////////////////////////////////////////////////////////////////////

//	Draw FitResult for Total Efficiency
//	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
//	const int nPar = 10;
	const int nPar = 10;
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
	f1_model->SetParameter(0, readParam(iBin,"acc", 1));
	f1_model->SetParameter(1, readParam(iBin,"acc", 2));
	f1_model->SetParameter(2, readParam(iBin,"acc", 3));
	f1_model->SetParameter(3, readParam(iBin,"Dimuonreco", 0));
	f1_model->SetParameter(4, readParam(iBin,"Dimuonreco", 1));
	f1_model->SetParameter(5, readParam(iBin,"Dimuonreco", 2));
	f1_model->SetParameter(6, readParam(iBin,"Dimuonreco", 3));
//	f1_model->FixParameter(7, 0.0);
//	f1_model->FixParameter(8, 0.0);
//	f1_model->FixParameter(9, 0.0);
	if ( iBin == 0 || iBin == 9 ) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 1));
		f1_model->FixParameter(1, readParam(iBin,"acc", 2));
		f1_model->FixParameter(2, readParam(iBin,"acc", 3));
		f1_model->FixParameter(3, readParam(iBin,"Dimuonreco", 0));
		f1_model->FixParameter(4, readParam(iBin,"Dimuonreco", 1));
		f1_model->FixParameter(5, readParam(iBin,"Dimuonreco", 2));
		f1_model->FixParameter(6, readParam(iBin,"Dimuonreco", 3));
		f1_model->FixParameter(7, readParam(iBin,"Dimuonreco", 4));
		f1_model->FixParameter(8, readParam(iBin,"Dimuonreco", 5));
		f1_model->FixParameter(9, readParam(iBin,"Dimuonreco", 6));		
/*		f1_model->SetParError(0, readParam(iBin,"accErr", 1));
		f1_model->SetParError(1, readParam(iBin,"accErr", 2));
		f1_model->SetParError(2, readParam(iBin,"accErr", 3));
		f1_model->SetParError(3, readParam(iBin,"DimuonrecoErr", 0));
		f1_model->SetParError(4, readParam(iBin,"DimuonrecoErr", 1));
		f1_model->SetParError(5, readParam(iBin,"DimuonrecoErr", 2));
		f1_model->SetParError(6, readParam(iBin,"DimuonrecoErr", 3));
		f1_model->SetParError(7, readParam(iBin,"DimuonrecoErr", 4));
		f1_model->SetParError(8, readParam(iBin,"DimuonrecoErr", 5));
		f1_model->SetParError(9, readParam(iBin,"DimuonrecoErr", 6));
*/	} else if (iBin == 1) {
		f1_model->FixParameter(0, readParam(iBin,"acc", 1));
		f1_model->FixParameter(1, readParam(iBin,"acc", 2));
		f1_model->FixParameter(2, readParam(iBin,"acc", 3));
		f1_model->FixParameter(3, readParam(iBin,"Dimuonreco", 0));
		f1_model->FixParameter(4, readParam(iBin,"Dimuonreco", 1));
		f1_model->FixParameter(5, readParam(iBin,"Dimuonreco", 2));
		f1_model->FixParameter(6, readParam(iBin,"Dimuonreco", 3));
		f1_model->FixParameter(7, readParam(iBin,"Dimuonreco", 4));
		f1_model->FixParameter(8, readParam(iBin,"Dimuonreco", 5));
		f1_model->FixParameter(9, readParam(iBin,"Dimuonreco", 6));		
	}else {
		f1_model->FixParameter(7, 0.0);
		f1_model->FixParameter(8, 0.0);
		f1_model->FixParameter(9, 0.0);
	}

	h2_EFF_fine.SetStats(0);
	h2_EFF_fine.SetMinimum(-0.0001);
//	h2_EFF_fine.SetMinimum(0.);
	h2_EFF_fine.SetTitleOffset(1.15,"Y");
	h2_EFF_fine.SetXTitle("cos#theta_{l}^{reco}");
	h2_EFF_fine.SetYTitle("Efficiency / 0.1");
	h2_EFF_fine.SetMaximum(h2_EFF_fine.GetMaximum() * 1.2); 
//	h2_EFF_fine.Draw("TEXT");
	h2_EFF_fine.Draw("PE1");
	
	h2_EFF_fine.Fit(f1_model,"R"); //// 09-09
//	TFitResultPtr r = h2_EFF_fine.Fit(f1_model,"WL S R");
//	r->Print();
	
	f1_model->SetTitle("");
	f1_model->SetLineWidth(2);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
    paveText->Draw();
	t1->DrawLatex(.22,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_EFF_fine_bin%d.pdf",iBin));
/////////////////////////////////////////////////////////////////////////////////////////////////
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
//	Draw compare
	TH1F h2_compFit("h2_compFit","",nLBins,-1.,1.);
	h2_compFit.SetXTitle("cos#theta_{l}^{reco}");
	TH1F h2_pullFit("h2_pullFit","",nLBins,-1.,1.);
	h2_pullFit.SetXTitle("cos#theta_{l}^{reco}");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		if (h2_EFF_fine.GetBinContent(i) != 0){
			h2_compFit.SetBinContent(i,f1_model->Eval(h2_EFF_fine.GetXaxis()->GetBinCenter(i))/h2_EFF_fine.GetBinContent(i));
			double _xlo = h2_EFF_fine.GetXaxis()->GetBinLowEdge(i);
			double _xhi = h2_EFF_fine.GetXaxis()->GetBinUpEdge(i);
			h2_pullFit.SetBinContent(i,(f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_EFF_fine.GetBinContent(i))/h2_EFF_fine.GetBinError(i));
		}else{
			h2_compFit.SetBinContent(i,0.);
			h2_pullFit.SetBinContent(i,0.);
		}
	}
	h2_compFit.SetMinimum(0.);
	h2_compFit.SetStats(0);
	h2_compFit.Draw("PE1");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.60,0.91,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_compFit_bin%d.pdf",iBin));
	
	h2_pullFit.SetStats(0);
	h2_pullFit.Draw(" HIST TEXT");
	latex->DrawLatexNDC(0.01,0.91,TString::Format("#chi^{2} = %f",chi2Val));
	latex->DrawLatexNDC(0.50,0.91,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_pullFit_bin%d.pdf",iBin));
	
//	Draw significance of deviation
	TH1F h2_pull("PULL - Deviation/Error","",15,-3.,3.);
	h2_pull.SetXTitle("Significance of deviation");
	for (int i = 1; i <= nLBins; i++) {//thetaL
		double _xlo = h2_EFF_fine.GetXaxis()->GetBinLowEdge(i);
		double _xhi = h2_EFF_fine.GetXaxis()->GetBinUpEdge(i);
		if (h2_EFF_fine.GetBinContent(i) != 0){
			h2_pull.Fill((f1_model->Integral(_xlo,_xhi)/(_xhi-_xlo)-h2_EFF_fine.GetBinContent(i))/h2_EFF_fine.GetBinError(i));
		}
	}
	h2_pull.Draw("HIST E1");
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEffDimuon_sigma_bin%d.pdf",iBin));
		
	delete latex;
	
	std::vector<double> output;
	for (int iPar = 0; iPar < nPar; iPar++){
		output.push_back(arrPar[iPar]);
		output.push_back(arrParErr[iPar]);
		
		printf("%18.15f,",arrPar[iPar]);
		if (iPar+1 >= nPar) printf("\n");
	}
	for (int i = 0; i < output.size(); i=i+2) {
		printf("%18.15f,",output[i+1]);
		if (i+2 >= output.size()) printf("\n");
	}
//	return output;
	writeParam(iBin,"accXrecoEffDimuon",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEffDimuonErr",arrParErr,nPar);  // f1_model_format_2
//	return output.c_str();	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void accXrecoEffLimited(int iBin, double EffP1,  double EffP2,  double EffP3,  double EffP4,  double EffP5,  double EffP6,  double EffP7,  double EffP8,  double EffP9,  double EffP10, int Index)
//void accXrecoEffLimited(int iBin, double EffP1,  double EffP2,  double EffP3,  double EffP4,  double EffP5,  double EffP6,  double EffP7,  double EffP8,  double EffP9,  double EffP10)
{//{{{
	setTDRStyle();
	printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
	
	TFile f_eff(TString::Format("./RootFiles/Efficiency_%d.root",iBin));
//	TH1F *h1_acc_fine  = (TH1F*)f_acc.Get(TString::Format("h1_acc_fine_bin%d", iBin));
	TH1F *h2_eff_Fine = (TH1F*)f_eff.Get("h2_eff_fine"); 
	
	TString f1_model_format_2 ;
	if (iBin == 0 || iBin ==1 || iBin == 9 ) { 
		f1_model_format_2 = "( [0]*exp(-0.5* (((x-[1])/[2])**2)) ) * ( [3]+[4]*x+[5]*x**2+[6]*x**3+[7]*x**4+[8]*x**5+[9]*x**6 ) ";
	}else { 
		f1_model_format_2 = "[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]+[8]+[9]"; 
	}
//	Draw
	TCanvas canvas("canvas");
	TLatex *latex = new TLatex();
	const int nPar = 10;
	TF1 *f1_model = new TF1("f1_model", f1_model_format_2, -1., 1.);
		f1_model->FixParameter(0, EffP1);
		f1_model->FixParameter(1, EffP2);
		f1_model->FixParameter(2, EffP3);
		f1_model->FixParameter(3, EffP4);
		f1_model->FixParameter(4, EffP5);
		f1_model->FixParameter(5, EffP6);
		f1_model->FixParameter(6, EffP7);
		f1_model->FixParameter(7, EffP8);
		f1_model->FixParameter(8, EffP9);
		f1_model->FixParameter(9, EffP10);
	if ( iBin == 2) {
		f1_model->SetParLimits(6, 0.0031679-0.001, 0.0031679+0.001); 
	}

	h2_eff_Fine->SetStats(0);
//	h2_eff_Fine->SetMinimum(-0.00005);
	h2_eff_Fine->SetMinimum(0.);
	h2_eff_Fine->SetTitleOffset(1.15,"Y");
	h2_eff_Fine->SetXTitle("cos#theta_{l}^{reco}");
	h2_eff_Fine->SetYTitle("Efficiency / 0.1");
	h2_eff_Fine->SetMaximum(h2_eff_Fine->GetMaximum() * 1.2); 
	h2_eff_Fine->Draw("PE1");
	
//	h2_eff_Fine.Fit(f1_model,"R"); //// 09-09
	
	f1_model->SetTitle("");
	f1_model->SetLineWidth(2);
	f1_model->SetLineColor(2);
	f1_model->Draw(" SAME ");
	
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.22,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	
	canvas.Update();
	canvas.Print(TString::Format("./plots/accXrecoEff_Eff_fine_bin%d_Index%d.pdf",iBin,Index));
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
	}
	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
//	return output;
	writeEffP(iBin,TString::Format("EffP_%d",Index),   arrPar,   nPar);  // f1_model_format_2
/////////////////////////////////////////////////////////////////////////////////////////////////
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

void DimuonResolution(int iBin)
//std::vector<double> accXrecoEff(int iBin)
{//{{{
	setTDRStyle();
	double BMass = 0;
	double Mumumass = 0;
	double Mumumasserr = 0;
	double gQ2 = 0;
	double Q2 = 0;
	double gCosThetaL = 0;
	double CosThetaL = 0;
	double gmuppt = 0;
	double gmupeta= 0;
	double gmumpt = 0;
	double gmumeta= 0;

	ch->SetBranchStatus("Bmass"         , 1);
	ch->SetBranchStatus("Mumumass"      , 1);
	ch->SetBranchStatus("Mumumasserr"   , 1);
	ch->SetBranchStatus("genQ2"         , 1);
	ch->SetBranchStatus("Q2"            , 1);
	ch->SetBranchStatus("genCosThetaL"  , 1);
	ch->SetBranchStatus("CosThetaL"     , 1);
	ch->SetBranchStatus("genMu*"        , 1);
	ch->SetBranchAddress("Bmass"        , &BMass);
	ch->SetBranchAddress("Mumumass"     , &Mumumass);
	ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
	ch->SetBranchAddress("Q2"           , &Q2);
	ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
	ch->SetBranchAddress("genQ2"        , &gQ2);
	ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
	ch->SetBranchAddress("genMupPt"     , &gmuppt);
	ch->SetBranchAddress("genMupEta"    , &gmupeta);
	ch->SetBranchAddress("genMumPt"     , &gmumpt);
	ch->SetBranchAddress("genMumEta"    , &gmumeta);
	
//	Fill histograms
	TH1F h2_Dimuon_fine("h2_Dimuon_fine","",1000,-2.0,2.0);
	double fffff = 0;
	for (int entry = 0; entry < ch->GetEntries(); entry++) {
		ch->GetEntry(entry);
		if (gQ2 > Q2rangeup[iBin] || gQ2 <= Q2rangedn[iBin]) continue;
			if (iBin == 10) {  ///////////////////////////  2015-04-29
				if (gQ2 < Q2rangeup[3] && gQ2 > Q2rangedn[3]) continue;
				if (gQ2 < Q2rangeup[5] && gQ2 > Q2rangedn[5]) continue;
			}
		if (BMass != -999 ){    ////////////////////////   12-10 N.A.
			fffff = Q2 - gQ2;
			h2_Dimuon_fine.Fill(fffff);
		}
	}
	
//	Draw
	TCanvas canvas("canvas");
	canvas.SetGrid();
	TLatex *latex = new TLatex();
	TPaveText* paveText = new TPaveText( 0.16, 0.77, 0.26, 0.87, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
	
	gStyle->SetOptStat();
	h2_Dimuon_fine.SetFillColor(4);
	h2_Dimuon_fine.Draw(" b ");
//	h2_Dimuon_fine.SetXTitle("M^{reco}_{#mu^{+}#mu^{-}} - M^{gen}_{#mu^{+}#mu^{-}}");
	h2_Dimuon_fine.SetXTitle("q^{2}_{reco} - q^{2}_{gen}");
	h2_Dimuon_fine.SetYTitle("Events/0.004");
    paveText->Draw();
	canvas.Update();
	canvas.Print(TString::Format("./plots/DimuonResolution_fine_bin%d.pdf",iBin));
	
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
