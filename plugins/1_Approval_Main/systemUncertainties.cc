// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular_gen_R_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular_gen_R")
{//{{{
	setTDRStyle();
//	double effUpperBound[11]  = {3.e4, 6.e4, 10.e4, 1., 6.e4, 1., 3.5e4, 3.5e4, 4.e4, 14.e4, 35.e4};
	RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{l}^{gen}", -1., 1.);
	RooRealVar genQ2("genQ2","q^{2}",1.0,22.);
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3.);
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................

	RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));
//	accP0.setConstant(kTRUE);
//	accP1.setConstant(kTRUE);
//	accP2.setConstant(kTRUE);
//	accP3.setConstant(kTRUE);
	accP0.setError(readParam(iBin,"accErr", 0));
	accP1.setError(readParam(iBin,"accErr", 1));
	accP2.setError(readParam(iBin,"accErr", 2));
	accP3.setError(readParam(iBin,"accErr", 3));  
	RooArgSet f_sigA_argset(genCosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
	TString f_sigA_format;
	TString f_ang_format;
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-genCosThetaL*genCosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*genCosThetaL )";
	}
	TString f_acc_format = "( accP0 + accP1 *exp(-0.5*(((genCosThetaL-accP2)/accP3)**2)) ) ";
	f_sigA_argset.add(RooArgSet(f_acc_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_acc_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf f("f","",f_sig,nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(RooArgSet(afb, fh)),PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) {
	    std::vector<double> output;
	    output.push_back(fh.getVal());
	    output.push_back(fh.getError());
	    output.push_back(afb.getVal());
	    output.push_back(afb.getError());
	return output;
	}

//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	RooPlot* framecosl = genCosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetTitleOffset(1.15,"Y");
	framecosl->SetMinimum(0);
//	framecosl->SetMaximum(effUpperBound[iBin]);
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.2);
	framecosl->Draw();
	
////////////////////////////////////// scan   //////////////////////////////////.........................
//	t1->DrawLatex(.30,.83,TString::Format("F_{H}  =%4.4f ", 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi()  ));
//	t1->DrawLatex(.60,.83,TString::Format("A_{FB} =%4.4f ", (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )   ));
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f #pm %4.4f ",fh.getVal(),fh.getError()));
//	t1->DrawLatex(.55,.83,TString::Format("A_{FB} =%4.4f #pm %4.4f ",afb.getVal(),afb.getError()));
	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f %4.4f",fh.getVal(),fh.getAsymErrorHi(),fh.getAsymErrorLo()));
	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f %4.4f",afb.getVal(),afb.getAsymErrorHi(),afb.getAsymErrorLo()));
////////////////////////////////////// refit  //////////////////////////////////.........................
	t1->SetTextFont(12);
	t1->DrawLatex(.19,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
//	TPaveText* paveText = new TPaveText( 0.80, 0.76, 0.90, 0.86, "NDC" );
	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.27, 0.80, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.pdf",outfile,iBin,Index));
	
	delete t1;
	delete c;
	delete data;
	
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
	//	val[0] = fh.getVal();val[1] = fh.getError(); val[2] = fh.getError();
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, "genfh_R", val, 3);
		val[1]=0; val[2]=0;
	//	val[0] = afb.getVal();val[1] = afb.getError();val[2] = afb.getError();
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "genafb_R",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "genFCN_R", val);
		printf("genAfb_R[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"genafb_R",0), fabs(readParam(iBin,"genafb_R",2)), fabs(readParam(iBin,"genafb_R",1)));
		printf("genFh_R[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"genfh_R",0), fabs(readParam(iBin,"genfh_R",2)), fabs(readParam(iBin,"genfh_R",1)));
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeOutput(outfile,iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(outfile,iBin, Index, "fh", val);
		val[1]=0; val[2]=0;
		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );
		writeOutput(outfile,iBin, Index, "F_afb", val);
		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
		writeOutput(outfile,iBin, Index, "F_fh", val);
		val[0] = f_fitresult->minNll();
		writeOutput(outfile,iBin, Index, "FCN", val);
	}
//	write output
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
}//}}}

void angular_gen_R(const char outfile[] = "angular_gen_R")
{//{{{
	setTDRStyle();
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = { 0.00,   0.00};
	double conYAerr[2]= { 0.20,   0.20};
	double conYF[2]   = { 0.20,   0.20};
	double conYFerr[2]= { 0.30,   0.30};
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yderrfh[9], yuerrfh[9], yafb[9], yderrafb[9], yuerrafb[9]; 
	double yrecofh[9], yrecoderrfh[9], yrecouerrfh[9], yrecoafb[9], yrecoderrafb[9], yrecouerrafb[9]; 
//	Check input 
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	reco
		yrecoafb[i]      = readParam(ibin,"recoafb",0);
		yrecoderrafb[i]   = fabs(readParam(ibin,"recoafb",1));
		yrecouerrafb[i]   = fabs(readParam(ibin,"recoafb",2));
		yrecofh[i]       = readParam(ibin,"recofh",0);
		yrecoderrfh[i]   = fabs(readParam(ibin,"recofh",1));
		yrecouerrfh[i]   = fabs(readParam(ibin,"recofh",2));
		printf("recoAfb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecoafb[i],yrecouerrafb[i],yrecoderrafb[i]);
		printf("recoFh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecofh[i],yrecouerrfh[i],yrecoderrfh[i]);
	//	gen_R
		yafb[i]      = readParam(ibin,"genafb_R",0);
		yderrafb[i]   = fabs(readParam(ibin,"genafb_R",1));
		yuerrafb[i]   = fabs(readParam(ibin,"genafb_R",2));
		yfh[i]       = readParam(ibin,"genfh_R",0);
		yderrfh[i]   = fabs(readParam(ibin,"genfh_R",1));
		yuerrfh[i]   = fabs(readParam(ibin,"genfh_R",2));
		printf("genAfb_R[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("genFh_R [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i],yuerrfh[i],yderrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitle("");
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->Draw();
//	FH
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,0.5,"Y");  
	TGraphAsymmErrors *r_fh  = new TGraphAsymmErrors(7,x,yrecofh,xerr,xerr,yrecoderrfh,yrecouerrfh);
	r_fh->SetMarkerColor(4);
	r_fh->SetMarkerStyle(20);	
	r_fh->SetFillColor(4);
	r_fh->SetLineColor(4);
	r_fh->SetFillStyle(3001);
	r_fh->Draw("2");
	r_fh->Draw("P");
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(3);
	g_fh->SetMarkerStyle(24);	
	g_fh->SetFillColor(3);
	g_fh->SetLineColor(3);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();
	
	TLegend *leg =new TLegend(0.64,0.72,0.90,0.86);
	leg->AddEntry(r_fh," reco level"," lePf ");
	leg->AddEntry(g_fh," gen level(acceptance corrected) "," lePf ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.19,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/Angular_gen_R_reco_fh.pdf"));
	c->Clear();
	
//	AFB
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.2,0.2,"Y"); 
	frame->Draw();
	TGraphAsymmErrors *r_afb = new TGraphAsymmErrors(7,x,yrecoafb,xerr,xerr,yrecoderrafb,yrecouerrafb);
	r_afb->SetMarkerColor(4);
	r_afb->SetMarkerStyle(20);
	r_afb->SetFillColor(4);
	r_afb->SetLineColor(4);
	r_afb->SetFillStyle(3001);
	r_afb->Draw("2");
	r_afb->Draw("P");
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	g_afb->SetMarkerColor(3);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(3);
	g_afb->SetLineColor(3);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2"); 
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	
	leg->Draw();
	t1->DrawLatex(.19,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.54,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/Angular_gen_R_reco_afb.pdf"));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular_reco_R_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular_reco_R")
{//{{{
	setTDRStyle();
	RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{l}^{gen}", -1., 1.);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................
	
	RooRealVar nsig("nsig","nsig",1E6,1E1,1E9);
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));
//	accP0.setConstant(kTRUE);
//	accP1.setConstant(kTRUE);
//	accP2.setConstant(kTRUE);
//	accP3.setConstant(kTRUE);
	accP0.setError(readParam(iBin,"accErr", 0));
	accP1.setError(readParam(iBin,"accErr", 1));
	accP2.setError(readParam(iBin,"accErr", 2));
	accP3.setError(readParam(iBin,"accErr", 3));  
//	reco Efficiency
	RooRealVar recoP0("recoP0","recoP0",readParam(iBin,"reco", 0));
	RooRealVar recoP1("recoP1","recoP1",readParam(iBin,"reco", 1));
	RooRealVar recoP2("recoP2","recoP2",readParam(iBin,"reco", 2));
	RooRealVar recoP3("recoP3","recoP3",readParam(iBin,"reco", 3));   
	RooRealVar recoP4("recoP4","recoP4",readParam(iBin,"reco", 4));  
	RooRealVar recoP5("recoP5","recoP5",readParam(iBin,"reco", 5));   
	RooRealVar recoP6("recoP6","recoP6",readParam(iBin,"reco", 6)); 
//	recoP0.setConstant(kTRUE);
//	recoP1.setConstant(kTRUE);
//	recoP2.setConstant(kTRUE);
//	recoP3.setConstant(kTRUE);
//	recoP4.setConstant(kTRUE);
//	recoP5.setConstant(kTRUE);
//	recoP6.setConstant(kTRUE);
	recoP0.setError(readParam(iBin,"recoErr", 0));
	recoP1.setError(readParam(iBin,"recoErr", 1));
	recoP2.setError(readParam(iBin,"recoErr", 2));
	recoP3.setError(readParam(iBin,"recoErr", 3));  
	recoP4.setError(readParam(iBin,"recoErr", 4)); 
	recoP5.setError(readParam(iBin,"recoErr", 5)); 
	recoP6.setError(readParam(iBin,"recoErr", 6)); 
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));  
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5)); 
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6)); 
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7)); 
//	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8)); 
//	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9)); 
//	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10));   
//	effP0.setConstant(kTRUE);
//	effP1.setConstant(kTRUE);
//	effP2.setConstant(kTRUE);
//	effP3.setConstant(kTRUE);
//	effP4.setConstant(kTRUE);
//	effP5.setConstant(kTRUE);
//	effP6.setConstant(kTRUE);
//	effP7.setConstant(kTRUE);
//	effP8.setConstant(kTRUE);
//	effP9.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
//	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
//	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
//	effP10.setError(readParam(iBin,"accXrecoEffErr", 10));

//	RooArgSet f_effA_argset(CosThetaL);
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
//	f_effA_argset.add(RooArgSet(effP8, effP9));      
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	RooArgSet f_accA_argset(CosThetaL);
//	f_accA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
//	RooArgSet f_recoA_argset(CosThetaL);
//	f_recoA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
	
//	RooGenericPdf f_acc("f_acc", "accP0 + accP1 *exp(-0.5*((CosThetaL-accP2)/accP3)^2) ", f_accA_argset);     
//	RooGenericPdf f_reco("f_reco", "recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6", f_recoA_argset);     
//	RooGenericPdf f_sigA("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb));
//	//RooGenericPdf f_sig("f_sig", "0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL", RooArgSet(CosThetaL,fh,afb),TString::Format("fabs(%s) <= (%s)/2.",afb,fh);
	
//	RooProdPdf    f_eff("f_eff","", f_acc, f_reco);
//	//RooProdPdf    f_effXsig("f_effXsig","", f_acc, f_reco, f_sig);
//	RooProdPdf    f_sig("f_effXsig","", f_eff, f_sigA); 
//	RooExtendPdf  f("f","", f_sig, nsig);
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

	RooArgSet f_sigA_argset(genCosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format;
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-genCosThetaL*genCosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*genCosThetaL )";
	}
	if (iBin != 0 && iBin !=1 && iBin != 9 ) { 
//	if (iBin != 0 && iBin !=1 && iBin != 9 && iBin != 7) { 
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*genCosThetaL+effP2*genCosThetaL**2+effP3*genCosThetaL**3+effP4*genCosThetaL**4+effP5*genCosThetaL**5+effP6*genCosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((genCosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * genCosThetaL + recoP2 * genCosThetaL**2 + recoP3 * genCosThetaL**3 + recoP4 * genCosThetaL**4 + recoP5 * genCosThetaL**5 + recoP6 * genCosThetaL**6  )";
	//	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
	//	f_rec_format = "( effP0 *exp(-0.5*(((genCosThetaL-effP1)/effP2)**2)) ) * (effP3+effP4*genCosThetaL+effP5*genCosThetaL**2+effP6*genCosThetaL**3+effP7*genCosThetaL**4 )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf  f("f","", f_sig, nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,Q2),Q2range[iBin],0);    // 12-08
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(RooArgSet(afb, fh)),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Warnings(-1), PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}

//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
		
	RooPlot* framecosl = genCosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.3);
	framecosl->Draw();
		
////////////////////////////////////// scan  //////////////////////////////////.........................
//	t1->DrawLatex(.30,.83,TString::Format("F_{H}  = %6.4f ", 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi()  ));
//	t1->DrawLatex(.60,.83,TString::Format("A_{FB} = %6.4f ", (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )   ));
////////////////////////////////////// scan  //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f #pm %4.4f ",fh.getVal(),fh.getError()));
//	t1->DrawLatex(.55,.83,TString::Format("A_{FB} =%4.4f #pm %4.4f ",afb.getVal(),afb.getError()));
	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f -%4.4f",fh.getVal(),fh.getError(),fh.getError()));
	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f -%4.4f",afb.getVal(),afb.getError(),afb.getError()));
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f %4.4f",fh.getVal(),fh.getAsymErrorHi(),fh.getAsymErrorLo()));
//	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f %4.4f",afb.getVal(),afb.getAsymErrorHi(),afb.getAsymErrorLo()));
////////////////////////////////////// refit  //////////////////////////////////.........................
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
//	TPaveText* paveText = new TPaveText( 0.80, 0.76, 0.90, 0.86, "NDC" );
	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.27, 0.80, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.pdf",outfile,iBin,Index));
//	clear
	delete t1;
	delete c;
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {  // ReFit
		val[0] = fh.getVal();val[1] = fh.getError();val[2] = fh.getError();
	//	val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, "recoRfh", val, 3);
		val[1]=0; val[2]=0;
		val[0] = afb.getVal();val[1] = afb.getError(); val[2] = afb.getError();
	//	val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "recoRafb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "recoFCN", val);
		printf("reco_R_Afb[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recoRafb",0), fabs(readParam(iBin,"recoRafb",2)), fabs(readParam(iBin,"recoRafb",1)));
		printf(" reco_R_Fh[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recoRfh",0), fabs(readParam(iBin,"recoRfh",2)), fabs(readParam(iBin,"recoRfh",1)));
/*	} else if (Index == -2) {
		cout<<endl<<endl<<"This is for a data fitting test!!!"<<endl;
		cout<<"Fit Status  = "<<f_fitresult->status()<<endl;
		cout<<"Afb = "<<afb.getVal()<<" +- "<<afb.getError()<<endl;
		cout<<"Fh  = "<<fh.getVal()<<" +- "<<fh.getError()<<endl;
*/	} else {  // ScanFit
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeOutput(outfile,iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(outfile,iBin, Index, "fh", val);
		val[1]=0; val[2]=0;
		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );
		writeOutput(outfile,iBin, Index, "F_afb", val);
		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
		writeOutput(outfile,iBin, Index, "F_fh", val);
		val[0] = f_fitresult->minNll();
		writeOutput(outfile,iBin, Index, "FCN", val);
	}
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
}//}}}

void angular_reco_R(const char outfile[] = "angular_reco_R")
{//{{{
	setTDRStyle();
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = { 0.00,   0.00};
	double conYAerr[2]= { 0.20,   0.20};
	double conYF[2]   = { 0.20,   0.20};
	double conYFerr[2]= { 0.30,   0.30};
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yrecoRfh[9], yrecoRderrfh[9], yrecoRuerrfh[9], yrecoRafb[9], yrecoRderrafb[9], yrecoRuerrafb[9]; 
	double yfh[9], yderrfh[9], yuerrfh[9], yafb[9], yderrafb[9], yuerrafb[9]; 
//	Check input 
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	reco
		yafb[i]      = readParam(ibin,"recoafb",0);
		yderrafb[i]  = fabs(readParam(ibin,"recoafb",1));
		yuerrafb[i]  = fabs(readParam(ibin,"recoafb",2));
		yfh[i]       = readParam(ibin,"recofh",0);
		yderrfh[i]   = fabs(readParam(ibin,"recofh",1));
		yuerrfh[i]   = fabs(readParam(ibin,"recofh",2));
		printf("recoAfb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("recoFh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i],yuerrfh[i],yderrfh[i]);
	//	recoR
		yrecoRafb[i]      = readParam(ibin,"recoRafb",0);
		yrecoRderrafb[i]  = fabs(readParam(ibin,"recoRafb",1));
		yrecoRuerrafb[i]  = fabs(readParam(ibin,"recoRafb",2));
		yrecoRfh[i]       = readParam(ibin,"recoRfh",0);
		yrecoRderrfh[i]   = fabs(readParam(ibin,"recoRfh",1));
		yrecoRuerrfh[i]   = fabs(readParam(ibin,"recoRfh",2));
		printf("reco_R_Afb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecoRafb[i],yrecoRuerrafb[i],yrecoRderrafb[i]);
		printf("reco_R_Fh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecoRfh[i],yrecoRuerrfh[i],yrecoRderrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitle("");
	frame->Draw();
//	FH
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,0.5,"Y");  
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);	
	g_fh->SetFillColor(4);
	g_fh->SetLineColor(4);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	TGraphAsymmErrors *recoR_fh  = new TGraphAsymmErrors(7,x,yrecoRfh,xerr,xerr,yrecoRderrfh,yrecoRuerrfh);
	recoR_fh->SetMarkerColor(3);
	recoR_fh->SetMarkerStyle(24);
	recoR_fh->SetFillColor(3);
	recoR_fh->SetLineColor(3);
	recoR_fh->SetFillStyle(3001);
	recoR_fh->Draw("2");
	recoR_fh->Draw("P");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();
	
	TLegend *leg =new TLegend(0.64,0.72,0.90,0.86);
	leg->AddEntry(g_fh," reco level"," lePf ");
	leg->AddEntry(recoR_fh," reco level(gen angular) "," lePf ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/angular_recoR_reco_fh.pdf"));
	c->Update();
	c->Clear();
//	AFB
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.2,0.2,"Y"); 
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(4);
	g_afb->SetLineColor(4);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	TGraphAsymmErrors *recoR_afb = new TGraphAsymmErrors(7,x,yrecoRafb,xerr,xerr,yrecoRderrafb,yrecoRuerrafb);
	recoR_afb->SetMarkerColor(3);
	recoR_afb->SetMarkerStyle(24);
	recoR_afb->SetFillColor(3);
	recoR_afb->SetLineColor(3);
	recoR_afb->SetFillStyle(3001);
	recoR_afb->Draw("2");
	recoR_afb->Draw("P");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2"); 
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	
	leg->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Update();
	c->Print(TString::Format("./plots/angular_recoR_reco_afb.pdf"));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular_reco_D_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular_reco_D")
{//{{{
	setTDRStyle();
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1., 1.);
//	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar genQ2("genQ2","q^{2}",1.0,22.);
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................
	
	RooRealVar nsig("nsig","nsig",1E6,1E1,1E9);
	
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
//	Acceptance
	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));
//	accP0.setConstant(kTRUE);
//	accP1.setConstant(kTRUE);
//	accP2.setConstant(kTRUE);
//	accP3.setConstant(kTRUE);
	accP0.setError(readParam(iBin,"accErr", 0));
	accP1.setError(readParam(iBin,"accErr", 1));
	accP2.setError(readParam(iBin,"accErr", 2));
	accP3.setError(readParam(iBin,"accErr", 3));  
//	reco Efficiency
	RooRealVar recoP0("recoP0","recoP0",readParam(iBin,"reco", 0));
	RooRealVar recoP1("recoP1","recoP1",readParam(iBin,"reco", 1));
	RooRealVar recoP2("recoP2","recoP2",readParam(iBin,"reco", 2));
	RooRealVar recoP3("recoP3","recoP3",readParam(iBin,"reco", 3));   
	RooRealVar recoP4("recoP4","recoP4",readParam(iBin,"reco", 4));  
	RooRealVar recoP5("recoP5","recoP5",readParam(iBin,"reco", 5));   
	RooRealVar recoP6("recoP6","recoP6",readParam(iBin,"reco", 6)); 
//	recoP0.setConstant(kTRUE);
//	recoP1.setConstant(kTRUE);
//	recoP2.setConstant(kTRUE);
//	recoP3.setConstant(kTRUE);
//	recoP4.setConstant(kTRUE);
//	recoP5.setConstant(kTRUE);
//	recoP6.setConstant(kTRUE);
	recoP0.setError(readParam(iBin,"recoErr", 0));
	recoP1.setError(readParam(iBin,"recoErr", 1));
	recoP2.setError(readParam(iBin,"recoErr", 2));
	recoP3.setError(readParam(iBin,"recoErr", 3));  
	recoP4.setError(readParam(iBin,"recoErr", 4)); 
	recoP5.setError(readParam(iBin,"recoErr", 5)); 
	recoP6.setError(readParam(iBin,"recoErr", 6)); 
/////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////

/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",readParam(iBin,"accXrecoEff", 0));
	RooRealVar effP1("effP1","effP1",readParam(iBin,"accXrecoEff", 1));
	RooRealVar effP2("effP2","effP2",readParam(iBin,"accXrecoEff", 2));
	RooRealVar effP3("effP3","effP3",readParam(iBin,"accXrecoEff", 3));   
	RooRealVar effP4("effP4","effP4",readParam(iBin,"accXrecoEff", 4));  
	RooRealVar effP5("effP5","effP5",readParam(iBin,"accXrecoEff", 5)); 
	RooRealVar effP6("effP6","effP6",readParam(iBin,"accXrecoEff", 6)); 
	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7)); 
	RooRealVar effP8("effP8","effP8",readParam(iBin,"accXrecoEff", 8)); 
	RooRealVar effP9("effP9","effP9",readParam(iBin,"accXrecoEff", 9)); 
//	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10));   
//	effP0.setConstant(kTRUE);
//	effP1.setConstant(kTRUE);
//	effP2.setConstant(kTRUE);
//	effP3.setConstant(kTRUE);
//	effP4.setConstant(kTRUE);
//	effP5.setConstant(kTRUE);
//	effP6.setConstant(kTRUE);
//	effP7.setConstant(kTRUE);
//	effP8.setConstant(kTRUE);
//	effP9.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
//	effP10.setError(readParam(iBin,"accXrecoEffErr", 10));

//	RooArgSet f_effA_argset(CosThetaL);
//	f_effA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
//	f_effA_argset.add(RooArgSet(effP8, effP9));      
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format;
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
	}
	if (iBin != 0 && iBin !=1 && iBin != 9 ) { 
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
	//	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
	//	f_sigA_argset.add(RooArgSet(effP7, effP8, effP9));
	//	f_rec_format = "( effP0 *exp(-0.5*(((CosThetaL-effP1)/effP2)**2)) ) * ( effP3 + effP4 * CosThetaL + effP5 * CosThetaL**2 + effP6 * CosThetaL**3 + effP7 * CosThetaL**4 + effP8 * CosThetaL**5 + effP9 * CosThetaL**6)";
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf  f("f","", f_sig, nsig);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);    // 12-08
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,genQ2),TString::Format("(%s) && (%s)",genQ2range[iBin], CTL[0]),0);    // 12-08
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(RooArgSet(afb, fh)),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Warnings(-1), PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) {
	    std::vector<double> output;
	    output.push_back(fh.getVal());
	    output.push_back(fh.getError());
	    output.push_back(afb.getVal());
	    output.push_back(afb.getError());
	    return output;
	}

//	Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
		
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(100)); 
	f.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.25);
	framecosl->Draw();
		
////////////////////////////////////// scan  //////////////////////////////////.........................
//	t1->DrawLatex(.30,.83,TString::Format("F_{H}  = %6.4f ", 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi()  ));
//	t1->DrawLatex(.60,.83,TString::Format("A_{FB} = %6.4f ", (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )   ));
////////////////////////////////////// scan  //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f #pm %4.4f ",fh.getVal(),fh.getError()));
//	t1->DrawLatex(.55,.83,TString::Format("A_{FB} =%4.4f #pm %4.4f ",afb.getVal(),afb.getError()));
	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f -%4.4f",fh.getVal(),fh.getError(),fh.getError()));
	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f -%4.4f",afb.getVal(),afb.getError(),afb.getError()));
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f %4.4f",fh.getVal(),fh.getAsymErrorHi(),fh.getAsymErrorLo()));
//	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f %4.4f",afb.getVal(),afb.getAsymErrorHi(),afb.getAsymErrorLo()));
////////////////////////////////////// refit  //////////////////////////////////.........................
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
//	TPaveText* paveText = new TPaveText( 0.80, 0.76, 0.90, 0.86, "NDC" );
	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.27, 0.80, "NDC" );
	paveText->SetBorderSize(0);
	paveText->SetFillColor(19);
    paveText->AddText(Form("bin %d ", iBin));
    paveText->Draw();
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.pdf",outfile,iBin,Index));
//	clear
	delete t1;
	delete c;
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {  // ReFit
		val[0] = fh.getVal();val[1] = fh.getError();val[2] = fh.getError();
	//	val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, "recoDfh", val, 3);
		val[1]=0; val[2]=0;
		val[0] = afb.getVal();val[1] = afb.getError(); val[2] = afb.getError();
	//	val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "recoDafb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "recoDFCN", val);
		printf("recoDAfb[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recoDafb",0), fabs(readParam(iBin,"recoDafb",2)), fabs(readParam(iBin,"recoDafb",1)));
		printf(" recoDFh[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recoDfh",0), fabs(readParam(iBin,"recoDfh",2)), fabs(readParam(iBin,"recoDfh",1)));
/*	} else if (Index == -2) {
		cout<<endl<<endl<<"This is for a data fitting test!!!"<<endl;
		cout<<"Fit Status  = "<<f_fitresult->status()<<endl;
		cout<<"Afb = "<<afb.getVal()<<" +- "<<afb.getError()<<endl;
		cout<<"Fh  = "<<fh.getVal()<<" +- "<<fh.getError()<<endl;
*/	} else {  // ScanFit
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeOutput(outfile,iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(outfile,iBin, Index, "fh", val);
		val[1]=0; val[2]=0;
		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );
		writeOutput(outfile,iBin, Index, "F_afb", val);
		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
		writeOutput(outfile,iBin, Index, "F_fh", val);
		val[0] = f_fitresult->minNll();
		writeOutput(outfile,iBin, Index, "FCN", val);
	}
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
}//}}}

void angular_reco_D(const char outfile[] = "angular_reco_D")
{//{{{
	setTDRStyle();
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = { 0.00,   0.00};
	double conYAerr[2]= { 0.20,   0.20};
	double conYF[2]   = { 0.20,   0.20};
	double conYFerr[2]= { 0.30,   0.30};
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yrecoDfh[9], yrecoDderrfh[9], yrecoDuerrfh[9], yrecoDafb[9], yrecoDderrafb[9], yrecoDuerrafb[9]; 
	double yfh[9], yderrfh[9], yuerrfh[9], yafb[9], yderrafb[9], yuerrafb[9]; 
//	Check input 
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	reco
		yafb[i]      = readParam(ibin,"recoafb",0);
		yderrafb[i]  = fabs(readParam(ibin,"recoafb",1));
		yuerrafb[i]  = fabs(readParam(ibin,"recoafb",2));
		yfh[i]       = readParam(ibin,"recofh",0);
		yderrfh[i]   = fabs(readParam(ibin,"recofh",1));
		yuerrfh[i]   = fabs(readParam(ibin,"recofh",2));
	//	if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
	//	else { yderrfh[i] = yuerrfh[i]; }
		printf("recoAfb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("recoFh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i],yuerrfh[i],yderrfh[i]);
	//	recoD
		yrecoDafb[i]      = readParam(ibin,"recoDafb",0);
		yrecoDderrafb[i]  = fabs(readParam(ibin,"recoDafb",1));
		yrecoDuerrafb[i]  = fabs(readParam(ibin,"recoDafb",2));
		yrecoDfh[i]       = readParam(ibin,"recoDfh",0);
		yrecoDderrfh[i]   = fabs(readParam(ibin,"recoDfh",1));
		yrecoDuerrfh[i]   = fabs(readParam(ibin,"recoDfh",2));
	//	if (yrecoDuerrfh[i] > fabs(yrecoDfh[i])) { yrecoDderrfh[i] = fabs(yrecoDfh[i]);}
	//	else { yrecoDderrfh[i] = yrecoDuerrfh[i]; }
		printf("recoDAfb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecoDafb[i],yrecoDuerrafb[i],yrecoDderrafb[i]);
		printf("recoDFh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yrecoDfh[i],yrecoDuerrfh[i],yrecoDderrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitle("");
	frame->Draw();
//	FH
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,0.5,"Y");  
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(4);
	g_fh->SetMarkerStyle(20);	
	g_fh->SetLineColor(4);
	g_fh->SetFillColor(4);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	TGraphAsymmErrors *recoD_fh  = new TGraphAsymmErrors(7,x,yrecoDfh,xerr,xerr,yrecoDderrfh,yrecoDuerrfh);
	recoD_fh->SetMarkerColor(3);
	recoD_fh->SetMarkerStyle(24);
	recoD_fh->SetFillColor(3);
	recoD_fh->SetLineColor(3);
	recoD_fh->SetFillStyle(3002);
	recoD_fh->Draw("2");
	recoD_fh->Draw("P");
	
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();
	
	TLegend *leg =new TLegend(0.64,0.72,0.90,0.86);
	leg->AddEntry(g_fh," reco level"," lePf ");
	leg->AddEntry(recoD_fh," gen level(Q2 bins) "," lePf ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/angular_recoD_reco_fh.pdf"));
	c->Update();
	c->Clear();
//	AFB
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetAxisRange(-0.2,0.2,"Y"); 
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	g_afb->SetMarkerColor(4);
	g_afb->SetMarkerStyle(20);
	g_afb->SetFillColor(4);
	g_afb->SetLineColor(4);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	TGraphAsymmErrors *recoD_afb = new TGraphAsymmErrors(7,x,yrecoDafb,xerr,xerr,yrecoDderrafb,yrecoDuerrafb);
	recoD_afb->SetMarkerColor(3);
	recoD_afb->SetMarkerStyle(24);
	recoD_afb->SetFillColor(3);
	recoD_afb->SetLineColor(3);
	recoD_afb->SetFillStyle(3002);
	recoD_afb->Draw("2");
	recoD_afb->Draw("P");
	
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2"); 
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	
	leg->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Update();
	c->Print(TString::Format("./plots/angular_recoD_reco_afb.pdf"));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////

//std::vector<double> angular2D_Limited_bin(int iBin, float Iafb, float Ifh, int Index, double EffP1,  double EffP2,  double EffP3,  double EffP4,  double EffP5,  double EffP6,  double EffP7,  double EffP8,  double EffP9,  double EffP10,  const char outfile[] = "angular2D_Limited")
std::vector<double> angular2D_Limited_bin(int iBin, float Iafb, float Ifh, int Index, double EffP1,  double EffP2,  double EffP3,  double EffP4,  double EffP5,  double EffP6,  double EffP7,  double EffP8,  double EffP9,  double EffP10,  int idex, const char outfile[] = "angular2D_Limited")
//std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
	setTDRStyle();
	// Remark: You must use RooFit!! It's better in unbinned fit.
	//         Extended ML fit is adopted by Mauro, just follow!!
	//         Need some modification for accXrecoEff.
	cout<<endl<<"iBin = "<<iBin<<"   idex = "<<idex<<endl<<endl; 
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................
		
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
//	Total Efficiency
	RooRealVar effP0("effP0","effP0",EffP1);
	RooRealVar effP1("effP1","effP1",EffP2);
	RooRealVar effP2("effP2","effP2",EffP3);
	RooRealVar effP3("effP3","effP3",EffP4);   
	RooRealVar effP4("effP4","effP4",EffP5);  
	RooRealVar effP5("effP5","effP5",EffP6); 
	RooRealVar effP6("effP6","effP6",EffP7); 
	RooRealVar effP7("effP7","effP7",EffP8); 
	RooRealVar effP8("effP8","effP8",EffP9); 
	RooRealVar effP9("effP9","effP9",EffP10); 
	effP0.setConstant(kTRUE);
	effP1.setConstant(kTRUE);
	effP2.setConstant(kTRUE);
	effP3.setConstant(kTRUE);
	effP4.setConstant(kTRUE);
	effP5.setConstant(kTRUE);
	effP6.setConstant(kTRUE);
	effP7.setConstant(kTRUE);
	effP8.setConstant(kTRUE);
	effP9.setConstant(kTRUE);
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
	// // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.26,5.30);
//	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.269,5.289);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
//	sigGauss1_sigma.setConstant(kTRUE);
	sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
//	sigGauss2_sigma.setConstant(kTRUE);
	sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
	RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
//	sigM_frac.setConstant(kTRUE);
	sigM_frac.setError(readParam(iBin,"sigM_frac",1));
	// // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
	
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format;
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
	}	
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_sigA_argset.add(RooArgSet(effP7, effP8, effP9));
		f_rec_format = "( effP0 *exp(-0.5*(((CosThetaL-effP1)/effP2)**2)) ) * ( effP3 + effP4 * CosThetaL + effP5 * CosThetaL**2 + effP6 * CosThetaL**3 + effP7 * CosThetaL**4 + effP8 * CosThetaL**5 + effP9 * CosThetaL**6  )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);
	
	// Create signal distribution
	RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_sig prepared. <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	
	// Create combinatorial background distribution

//	Build Chebychev polynomial p.d.f.  
	RooRealVar bkgCombM_c("bkgCombM_c","c",-0.5,-20.,1.);
	RooRealVar offset("offset","offset",-5.);
	RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
	RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
	
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",readParam(iBin,"bkgCombL_c0",0));
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0));
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0));
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0)); 
//	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0));
//	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",readParam(iBin,"bkgCombL_c5",0),-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
//			bkgCombL_c0.setConstant(kTRUE);
//			bkgCombL_c1.setConstant(kTRUE);
//			bkgCombL_c2.setConstant(kTRUE);
			//if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
			if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
//				bkgCombL_c3.setConstant(kTRUE);
			}
			//bkgCombL_c4.setConstant(kTRUE);
			//bkgCombL_c5.setConstant(kTRUE);
			break;
	}
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	// Create peak background distribution
	RooRealVar bkgGauss_mean("bkgGauss_mean","cos#theta_{l}", readParam(iBin,"bkgGauss_mean",0));
	RooRealVar bkgGauss_sigma("bkgGauss_sigma","#sigma",  readParam(iBin,"bkgGauss_sigma",0));
	RooRealVar bkg_frac("bkg_frac","bkg_frac",readParam(iBin,"bkg_frac",0));
	
	RooGaussian f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", CosThetaL, bkgGauss_mean, bkgGauss_sigma);
	RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombL_P), bkg_frac);
	
	RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL, f_bkgCombM);
	cout<<">>>>>>>>>>>>>>>> INFO: f_bkgComb prepared. <<<<<<<<<<<<<<<<<<<<<<"<<endl;

	RooRealVar nsig("nsig","nsig",500,0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",500,0,6E3);
//	nbkgComb.setConstant(kTRUE);
	
	RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_sig),RooArgList(nbkgComb,nsig));// no penalty term
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	

///////////////////////////////////////////////////////////// Gaussian constraints ///////////////////////////////////////////////////	
	// Gaussian constraints
	RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",sigGauss1_sigma,RooConst(readParam(iBin,"sigGauss1_sigma",0)),RooConst(readParam(iBin,"sigGauss1_sigma",1)));
	RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",sigGauss2_sigma,RooConst(readParam(iBin,"sigGauss2_sigma",0)),RooConst(readParam(iBin,"sigGauss2_sigma",1)));
	RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",sigM_frac,RooConst(readParam(iBin,"sigM_frac",0)),RooConst(readParam(iBin,"sigM_frac",1)));
	
	RooArgSet gausConstraints(gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
///////////////////////////////////////////////////////////// Gaussian constraints ///////////////////////////////////////////////////	
	// Get data and apply unbinned fit
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
//	ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Minos(kTRUE));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit2","Simplex"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit2"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3));	
	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Strategy(2), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2","Simplex"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3) {
		if (Index == -1) {
			double val_1[4]={0,0,0,0};
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("fh_L_%d",idex), val_1);
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("afb_L_%d",idex),val_1);
			val_1[1]= 0;
			writeEffP(iBin, TString::Format("FCN_L_%d",idex), val_1);
		}
//	if (f_fitresult->status() != 0) {
	    std::vector<double> output;
	    output.push_back(fh.getVal());
	    output.push_back(fh.getError());
	    output.push_back(afb.getVal());
	    output.push_back(afb.getError());
	    return output;
	}
/*

	double UpperBound_mass[11] ={ 50,100,190,0.,110,0., 55, 70, 50,220,600};
	double UpperBound_cosl[11] ={ 50,100,160,0., 80,0., 40, 40, 35,200,400};


	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	TPad *tp1 = new TPad("tp1","",0.0,0.15,1.0,1.00);
	TPad *tp2 = new TPad("tp2","",0.0,0.00,1.0,0.15);
	tp1->Draw();
	tp2->Draw();
	tp1->cd();
	RooPlot* framemass = Bmass.frame();
	data->plotOn(framemass,RooFit::Name("data"), Binning(25)); 
	f.plotOn(framemass,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullmass = framemass->pullHist();
	//f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
//	f.plotOn(framemass, Components(f_sig),FillStyle(3001),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framemass, Components(f_sig),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(2),LineWidth(2));
	f.plotOn(framemass,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(5),LineWidth(4),LineStyle(4));
	
	framemass->SetTitle("");
	framemass->SetTitleOffset(1.1,"Y");
	framemass->SetMinimum(0);
	framemass->SetMaximum(UpperBound_mass[iBin]);
	framemass->Draw();
    
	TLegend *leg =new TLegend(0.69,0.69,0.90,0.86,NULL,"brNDC");
	leg->AddEntry("data"," Data "," PE ");
	leg->AddEntry("pdf"," Total P.d.f. "," L ");
	leg->AddEntry("sig"," Signal "," L ");
	leg->AddEntry("bkgComb"," Comb. bkg. "," L");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->Draw();

	TPaveText* paveText = new TPaveText( 0.17, 0.65, 0.38, 0.88, "NDC" ); 
//	TPaveText* paveText = new TPaveText( 0.17, 0.70, 0.41, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
////////////////////////////////////// scan  //////////////////////////////////.........................
//	paveText->AddText(Form("F_{H}  =%4.4f ", 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )); 
//	paveText->AddText(Form("A_{FB} =%4.4f ", (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )  )); 
////////////////////////////////////// scan  //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	paveText->AddText(Form("F_{H}  =%6.4f #pm %6.4f ", fh.getVal(),fh.getError())); 
	paveText->AddText(Form("A_{FB} =%6.4f #pm %6.4f ", afb.getVal(),afb.getError())); 
//	paveText->AddText(Form("F_{H}  =%6.4f + %6.4f  %6.4f", fh.getVal(),fh.getAsymErrorHi(),fh.getAsymErrorLo())); 
//	paveText->AddText(Form("A_{FB} =%6.4f + %6.4f  %6.4f", afb.getVal(),afb.getAsymErrorHi(),afb.getAsymErrorLo())); 
////////////////////////////////////// refit  //////////////////////////////////.........................
	paveText->AddText(Form("   N_{sig}   = %.1f #pm %.1f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("   N_{bkg}   = %.1f #pm %.1f ", nbkgComb.getVal(), nbkgComb.getError())); 
//	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
	paveText->AddText(Form("#Chi^{2}_{Bmass}   = %.2f  ", framemass->chiSquare())); 
	paveText->AddText(Form("Fit Status   = %d  ", f_fitresult->status() )); 
	paveText->Draw(); 
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.66,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
//	TPaveText* paveText_l = new TPaveText( 0.20, 0.78, 0.26, 0.84, "NDC" );
//	TPaveText* paveText_l = new TPaveText( 0.80, 0.76, 0.90, 0.86, "NDC" );
	TPaveText* paveText_l = new TPaveText( 0.81, 0.58, 0.89, 0.68, "NDC" );
	paveText_l->SetBorderSize(0);
	paveText_l->SetFillColor(19);
    paveText_l->AddText(Form("bin %d ", iBin));
    paveText_l->Draw();
	
	tp2->cd();
	RooPlot *framepullmass = Bmass.frame(Title("  "));
	framepullmass->addPlotable(pullmass,"P");
	framepullmass->GetYaxis()->SetRangeUser(-3.,3.);
	framepullmass->GetYaxis()->SetNdivisions(6);
	framepullmass->GetXaxis()->SetTitle("  ");
	framepullmass->SetLabelFont(22,"XY");
	framepullmass->SetLabelSize(0.12,"XY");
	TLine *tl1 =new TLine(5.10, 0., 5.60, 0.);
	tl1->SetLineColor(2);
	framepullmass->Draw();
	tl1->Draw();

	c->Update();
	c->Print(TString::Format("./plots/%s_bin%d_idex_%d_Index_%d.pdf",outfile,iBin,idex,Index));
	
	// Draw projection to CosThetaL
	tp1->cd();
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,RooFit::Name("data"), Binning(25)); 
	f.plotOn(framecosl,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullcosl = framecosl->pullHist();
	//f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
	//f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
//	f.plotOn(framecosl, Components(f_sig),FillStyle(3001),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framecosl, Components(f_sig),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(2),LineWidth(2));
	f.plotOn(framecosl,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(5),LineWidth(4),LineStyle(4));
	
	framecosl->SetTitle("");
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(UpperBound_cosl[iBin]);
	framecosl->Draw();
	leg->Draw();
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.66,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	paveText->AddText(Form("#Chi^{2}_{cos#theta}   = %.2f  ", framecosl->chiSquare())); 
	paveText->Draw(); 
    paveText_l->Draw();

	tp2->cd();
	RooPlot *framepullcosl = CosThetaL.frame(Title("  "));
	framepullcosl->addPlotable(pullcosl,"P");
	framepullcosl->GetYaxis()->SetRangeUser(-5.,5.);
	framepullcosl->GetYaxis()->SetNdivisions(10);
	framepullcosl->GetXaxis()->SetTitle("  ");
	framepullcosl->SetLabelFont(22,"XY");
	framepullcosl->SetLabelSize(0.12,"XY");
	TLine *tl2 =new TLine(-1., 0., 1., 0.);
	tl2->SetLineColor(2);
	framepullcosl->Draw();
	tl2->Draw();      
	
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d_idex_%d_Index_%d.pdf",outfile,iBin,idex,Index));
	delete c;
	delete t1;
*/
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = fh.getVal();val[1] = fh.getError();val[2] = fh.getError();
	//	val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeEffP(iBin, TString::Format("fh_L_%d",idex), val, 3);
		val[0] = afb.getVal();val[1] = afb.getError();val[2] = afb.getError();
	//	val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeEffP(iBin, TString::Format("afb_L_%d",idex),val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeEffP(iBin, TString::Format("FCN_L_%d",idex), val);
/*		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
		val[1]= 3./2. + 3. * atan( fh.getVal() + fh.getError() ) / TMath::Pi() - val[0]; 
		val[2]= 3./2. + 3. * atan( fh.getVal() - fh.getError() ) / TMath::Pi() - val[0]; 
		writeParam(iBin, "fh", val, 3);
		double FH, FHu, FHd;
		FH = val[0]; FHu = val[1]; FHd = val[2];
		val[1]=0; val[2]=0;
		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );
		val[1]=   sqrt( pow( ( 1. * FH * atan( afb.getVal() + afb.getError() ) / TMath::Pi() - val[0] ), 2) + pow( ( 1. * atan( afb.getVal() ) * FHu / TMath::Pi() ), 2)  + 2. * ( 1. * FH * atan( afb.getVal() + afb.getError() ) / TMath::Pi() - val[0] ) / afb.getError() * ( 1. * atan( afb.getVal() ) * FHu / TMath::Pi() ) / fh.getError() * f_fitresult->correlation(afb, fh) );
		val[2]= - sqrt( pow( ( 1. * FH * atan( afb.getVal() - afb.getError() ) / TMath::Pi() - val[0] ), 2) + pow( ( 1. * atan( afb.getVal() ) * FHd / TMath::Pi() ), 2)  + 2. * ( 1. * FH * atan( afb.getVal() - afb.getError() ) / TMath::Pi() - val[0] ) / afb.getError() * ( 1. * atan( afb.getVal() ) * FHd / TMath::Pi() ) / fh.getError() * f_fitresult->correlation(afb, fh) );
		writeParam(iBin, "afb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "FCN", val);
		val[0] = sigGauss_mean.getVal(); val[1] = sigGauss_mean.getError();
		writeParam(iBin, "sigGauss_mean", val);
		val[0] = bkgCombM_c.getVal(); val[1] = bkgCombM_c.getError();
		writeParam(iBin, "bkgCombM_c", val);
*/	} else if (Index == -2) {
		cout<<endl<<endl<<"This is for a data fitting test!!!"<<endl;
		cout<<"Fit Status  = "<<f_fitresult->status()<<endl;
		cout<<"Afb = "<<afb.getVal()<<" +- "<<afb.getError()<<endl;
		cout<<"Fh  = "<<fh.getVal()<<" +- "<<fh.getError()<<endl;
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeEffPOutput(outfile,iBin, idex, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeEffPOutput(outfile,iBin, idex, Index, "fh", val);
		val[1]=0; val[2]=0;
		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrained
	//	val[0] = afb.getVal();   // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_afb", val);
		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrained
	//	val[0] = fh.getVal();    // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_fh", val); 
		val[0] = f_fitresult->minNll();
		writeEffPOutput(outfile,iBin, idex, Index, "FCN", val);
		val[0] = idex;
		writeEffPOutput(outfile,iBin, idex, Index, "idex", val);
	}
	
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}

void PlotFCN_Limited( int iBin, const char outfile[] = "FCN")
{
	setTDRStyle();
//for (int idex = 1; idex <= 200; idex++) {	
    for (int idex = 1; idex <= 200; idex++) {	
	    TGraph *gr_afb = new TGraph();
	    TGraph *gr_fh  = new TGraph();
	    double FCN = 999, fcn = 999;
	    int Index = 0, index = 0, n = -1, NIndex = 1000;
	    double Iafb, Ifh, afb, fh;
	    do {
		    index+=1;
		    while ( ! (fcn = readEffPOutput(outfile, iBin, idex, index, "FCN", 0)) && index < NIndex) index+=1;
		    afb = readEffPOutput(outfile, iBin, idex, index, "F_afb", 0);  
		    fh  = readEffPOutput(outfile, iBin, idex, index, "F_fh", 0);
	    //	if (fcn < 0 && (afberr > 0.001) && (fherr > 0.001) ) {
		    if (fcn < 0 ) {
			    n+=1;
			    gr_afb->SetPoint(n, afb, fcn);
			    gr_fh->SetPoint(n, fh, fcn);
		    }
		    if ( FCN > fcn ) { FCN = fcn; Index = index; }
	    //	if ( FCN > fcn && (fabs(afb) < (fh / 2.))) { FCN = fcn; Index = index; }
	    } while ( index < NIndex);
	    Iafb   = readEffPOutput(outfile, iBin, idex, Index, "afb", 0);
	    Ifh   = readEffPOutput(outfile, iBin, idex, Index, "fh", 0);
	    afb    = readEffPOutput(outfile, iBin, idex, Index, "F_afb", 0);  
	    fh    = readEffPOutput(outfile, iBin, idex, Index, "F_fh", 0);  
	    TCanvas *c = new TCanvas("c","c",800,600);
//	    c->SetTitle("FCN distribution of A_{FB}");
	    gr_afb->GetXaxis()->SetTitle("A_{FB}");
	    gr_afb->GetYaxis()->SetTitle("NLL");
	    gr_afb->Draw("AP");
//	    c->Print(TString::Format("./plots/%s_FCN_afb_bin%d_EffP_%d.pdf",outfile,iBin,idex));
	    c->Clear();
//	    c->SetTitle("FCN distribution of F_{H}");
	    gr_fh->GetXaxis()->SetTitle("F_{H}");
	    gr_fh->GetYaxis()->SetTitle("NLL");
	    gr_fh->Draw("AP");
//	    c->Print(TString::Format("./plots/%s_FCN_fh_bin%d_EffP_%d.pdf",outfile,iBin,idex));
	    c->Clear();
	    delete c;
	    cout<<"EffP_"<<idex<<endl;
	    cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
//	    cout<<"Iafb  = "<<Iafb<<"   Ifh = "<<Ifh<<endl;
	    cout<<"afb = "<<afb<<endl;
	    cout<<"fh  = "<<fh<<endl;
	    double val[3]={0,0,0};
	    val[0] = Iafb; val[1] = readEffPOutput(outfile, iBin, idex, Index, "afb", 1);
	    writeEffP(iBin, TString::Format("Iafb_%s_EffP_%d",outfile,idex),val);
	    val[1]=0; val[2]=0;
	    val[0] = Ifh; val[1] = readEffPOutput(outfile, iBin, idex, Index, "fh", 1);
	    writeEffP(iBin, TString::Format("Ifh_%s_EffP_%d",outfile,idex),val);
	    val[1]=0; val[2]=0;
	    val[0] = afb;
	    writeEffP(iBin, TString::Format("afb_%s_EffP_%d",outfile,idex),val);
	    val[1]=0; val[2]=0;
	    val[0] = fh;
	    writeEffP(iBin, TString::Format("fh_%s_EffP_%d",outfile,idex),val);
    }
}

void angular2D_Limited_FCN( int iBin, int Index, const char outfile[] = "angular2D_Limited_FCN")
{//{{{
	setTDRStyle();
	TGraphAsymmErrors *gr   = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_0 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_1 = new TGraphAsymmErrors();
	TCanvas *c = new TCanvas("c","c",800,600);
	TLatex *t1 = new TLatex();
	for (float n = 0, x =-1., y = 0; x <= 1.; n+=1, x+=0.0005) {
		y = 2 * fabs(x);
		gr_0->SetPoint(n, x, 0);
		gr_1->SetPoint(n, x, 3);
		if (x > -1.0 && x < 1.0) gr->SetPoint(n, x, y);
	}
	gr->SetMarkerColor(2);
	gr->SetMarkerSize(0.1);
	gr->SetMarkerStyle(7);
	gr->GetYaxis()->SetTitleOffset(1.15);
	gr_0->SetLineColor(4);
	gr_1->SetLineColor(4);
	gr_0->GetXaxis()->SetTitle("A_{FB}");
	gr_0->GetYaxis()->SetTitle("F_{H}");
	gr_0->GetXaxis()->SetRangeUser(-0.40,0.40);  // data fitted
	gr_0->GetYaxis()->SetRangeUser(-0.05,0.80);
	if (iBin == 1 || iBin == 999 || iBin == 998) {
        gr_0->GetXaxis()->SetRangeUser(-0.60,0.60);  // bin1
	    gr_0->GetYaxis()->SetRangeUser(-0.05,1.20);
    }
	if (iBin == 2) {
        gr_0->GetXaxis()->SetRangeUser(-0.005,0.005);  // bin2
	    gr_0->GetYaxis()->SetRangeUser(-0.001,0.01);
    }
	if (iBin == 6 || iBin == 7 || iBin == 8) {
        gr_0->GetXaxis()->SetRangeUser(-0.15,0.15);  // bin4
	    gr_0->GetYaxis()->SetRangeUser(-0.01,0.30);
    }
	if (iBin == 4 || iBin == 10) {
        gr_0->GetXaxis()->SetRangeUser(-0.05,0.05);  // bin4
	    gr_0->GetYaxis()->SetRangeUser(-0.01,0.10);
    }
//	gr_0->GetXaxis()->SetRangeUser(-0.008,0.008);  // data fitted
//	gr_0->GetYaxis()->SetRangeUser(-0.001,0.016);
	gr_0->Draw("AL");	
	gr_1->Draw("L");	
	gr->Draw("P");	
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.66,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
//	t1->DrawLatex(-0.4,.81,TString::Format("CMS Preliminary"));
//	t1->DrawLatex( 0.1,.81,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	
	TGraphAsymmErrors *gr_afb_fh = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_0 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_1 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_2 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_4 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_6 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_7 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_8 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_9 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_afb_fh_10 = new TGraphAsymmErrors();
	double fcn = 999;
	int index = 0, n = -1, NIndex = 100;
	double afb[200], fh[200];
	double rms_afb[11], rms_fh[11], SQ_afb[11], SQ_fh[11];
	double sd_afb[11], sd_fh[11], sum_afb[11], sum_fh[11];
	double avg_afb[11], avg_fh[11], sq_afb[11], sq_fh[11];
	double aa, bb;
	if (Index == -999) {
        for (int idex = 1; idex <= 200; idex++) {
    	    index = 0;
    		do {
    			index+=1;
    			while ( ! (fcn =  readEffPOutput(outfile, iBin, idex, index, "FCN", 0) )  && index < NIndex) index+=1;
    			n+=1;
    			aa = readEffPOutput(outfile, iBin, idex, index, "F_afb", 0);  // best fitted values
    			bb = readEffPOutput(outfile, iBin, idex, index, "F_fh", 0);
    			if (fcn < 0 ) {
    				gr_afb_fh->SetPoint(n, aa, bb);
    			}
    		} while ( index < NIndex);
        }
	} else if (Index == -1) {
		n = -1;
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
		for (int idex = 1; idex <= 200; idex++) {
			n+=1;
			afb[idex] = readEffP(ibin, TString::Format("afb_%s_EffP_%d",outfile,idex), 0);  // best fitted values!
			fh[idex]  = readEffP(ibin, TString::Format("fh_%s_EffP_%d",outfile,idex), 0);
			gr_afb_fh->SetPoint(n, afb[idex], fh[idex]);
		}
			t1->SetTextColor(4);
			t1->DrawLatex(afb[199]+0.1, fh[199]-0.1, TString::Format("%d", ibin));
		}
	} else if (Index == -2) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
			n = -1;
			SQ_afb[ibin] = 0;
			SQ_fh[ibin]  = 0;
			sum_afb[ibin] = 0;
			sum_fh[ibin]  = 0;
			for (int idex = 1; idex <= 200; idex++) {
				n+=1;
				aa = readEffP(ibin, TString::Format("afb_L_%d",idex), 0);  // best fitted values!
				if ( aa == -999) n=n-1;
				else {
					afb[idex] = readEffP(ibin, TString::Format("afb_L_%d",idex), 0);
					fh[idex]  = readEffP(ibin, TString::Format("fh_L_%d",idex), 0);
					SQ_afb[ibin] = SQ_afb[ibin] + afb[idex] * afb[idex];
					SQ_fh[ibin]  = SQ_fh[ibin]  + fh[idex]  * fh[idex];
					sum_afb[ibin] = sum_afb[ibin] + afb[idex];
					sum_fh[ibin]  = sum_fh[ibin]  + fh[idex];
					if (ibin == 0 ) {
						gr_afb_fh_0->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 1) {
						gr_afb_fh_1->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 2) {
						gr_afb_fh_2->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 4) {
						gr_afb_fh_4->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 6) {
						gr_afb_fh_6->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 7) {
						gr_afb_fh_7->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 8) {
						gr_afb_fh_8->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 9) {
						gr_afb_fh_9->SetPoint(n, afb[idex], fh[idex]);
					} else if (ibin == 10) {
						gr_afb_fh_10->SetPoint(n, afb[idex], fh[idex]);
					}
				}
			}
			avg_afb[ibin] = sum_afb[ibin] / (n+1);
			avg_fh[ibin] = sum_fh[ibin] / (n+1);
			rms_afb[ibin] = sqrt( SQ_afb[ibin] / (n+1)  );
			rms_fh[ibin]  = sqrt( SQ_fh[ibin]  / (n+1)  );
			cout<<"N_"<<ibin<<" = "<<n+1<<endl;
		//	cout<<"rms_afb"<<ibin<<" = "<<rms_afb[ibin]<<endl;
		//	cout<<"rms_fh"<<ibin<<" = "<<rms_fh[ibin]<<endl;
			n = -1;
			sq_afb[ibin] = 0;
			sq_fh[ibin]  = 0;
			for (int idex = 1; idex <= 200; idex++) {
				n+=1;
				aa = readEffP(ibin, TString::Format("afb_L_%d",idex), 0);  // best fitted values!
				if ( aa == -999) n=n-1;
				else {
					afb[idex] = readEffP(ibin, TString::Format("afb_L_%d",idex), 0);
					fh[idex]  = readEffP(ibin, TString::Format("fh_L_%d",idex), 0);
					sq_afb[ibin] = sq_afb[ibin] + (afb[idex]-avg_afb[ibin]) * (afb[idex]-avg_afb[ibin]);
					sq_fh[ibin]  = sq_fh[ibin] + (fh[idex]-avg_fh[ibin]) * (fh[idex]-avg_fh[ibin]);
				}
			}
			sd_afb[ibin] = sqrt( sq_afb[ibin] / n  );
			sd_fh[ibin]  = sqrt( sq_fh[ibin]  / n  );
			cout<<"sd_afb"<<ibin<<" = "<<sd_afb[ibin]<<"  avg_afb"<<ibin<<" = "<<avg_afb[ibin]<<endl;
			cout<<"sd_fh"<<ibin<<" = "<<sd_fh[ibin]<<"  avg_fh"<<ibin<<" = "<<avg_fh[ibin]<<endl;
		}
	}
	
	if (Index == -999) {
		TGraphAsymmErrors *gr_afb_fh_f = new TGraphAsymmErrors();
		for (int idex = 1; idex <= 200; idex++) {
			aa = readEffP(iBin, TString::Format("afb_%s_EffP_%d",outfile,idex), 0);  // best fitted values!
			bb  = readEffP(iBin, TString::Format("fh_%s_EffP_%d",outfile,idex), 0);
			gr_afb_fh_f->SetPoint(idex-1, aa, bb);
		}
	    gr_afb_fh->Draw("P");
	    gr_afb_fh->SetMarkerStyle(20);
		gr_afb_fh_f->SetMarkerColor(3);
		gr_afb_fh_f->SetMarkerStyle(24);
		gr_afb_fh_f->SetMarkerSize(1);
		gr_afb_fh_f->Draw(" P ");
		
	//	TLegend *leg =new TLegend(0.16,0.75,0.38,0.88,NULL,"brNDC");
		TLegend *leg =new TLegend(0.69,0.75,0.91,0.88,NULL,"brNDC");
		leg->AddEntry(gr_afb_fh,  " Succeeded fitted results ","P");
		leg->AddEntry(gr_afb_fh_f," Best fitted results ",    "P");
		leg->SetLineColor(1);
		leg->SetFillColor(0);
		leg->SetTextSize(0.02);
		leg->Draw();					
		c->Print(TString::Format("./plots/%s_afb_fh_bin%d_f.pdf",outfile,iBin));  // fitted
	} else if (Index == -1 || Index == -2) {
		if (Index == -1) {
	    gr_afb_fh->Draw("P");
	    gr_afb_fh->SetMarkerStyle(20);
			TLegend *leg_1 =new TLegend(0.69,0.75,0.91,0.88,NULL,"brNDC");
		//	TLegend *leg_1 =new TLegend(0.16,0.79,0.38,0.88,NULL,"brNDC");
			leg_1->AddEntry(gr_afb_fh," Best fitted result ", "P");
			leg_1->SetLineColor(1);
			leg_1->SetFillColor(0);
			leg_1->SetTextSize(0.02);
			leg_1->Draw();					
			c->Print(TString::Format("./plots/%s_afb_fh_bin%d_f.pdf",outfile,iBin));
		} else if (Index == -2){
    		gr_afb_fh_0->SetMarkerStyle(24);
    		gr_afb_fh_1->SetMarkerStyle(25);
    		gr_afb_fh_2->SetMarkerStyle(2);
    		gr_afb_fh_4->SetMarkerStyle(28);
    		gr_afb_fh_6->SetMarkerStyle(30);
    		gr_afb_fh_7->SetMarkerStyle(5);
    		gr_afb_fh_8->SetMarkerStyle(32);
    		gr_afb_fh_9->SetMarkerStyle(30);
    		gr_afb_fh_10->SetMarkerStyle(29);
    		gr_afb_fh_0->SetMarkerColor(3);
    		gr_afb_fh_1->SetMarkerColor(2);
    		gr_afb_fh_2->SetMarkerColor(5);
    		gr_afb_fh_4->SetMarkerColor(6);
    		gr_afb_fh_6->SetMarkerColor(7);
    		gr_afb_fh_7->SetMarkerColor(8);
    		gr_afb_fh_8->SetMarkerColor(9);
    		gr_afb_fh_9->SetMarkerColor(4);
    		gr_afb_fh_10->SetMarkerColor(1);
    		gr_afb_fh_0->Draw("P");
    		gr_afb_fh_1->Draw("P");
    		gr_afb_fh_2->Draw("P");
    		gr_afb_fh_4->Draw("P");
    		gr_afb_fh_6->Draw("P");
    		gr_afb_fh_7->Draw("P");
    		gr_afb_fh_8->Draw("P");
    		gr_afb_fh_9->Draw("P");
    		gr_afb_fh_10->Draw("P");
    		//	TLegend *leg_2 =new TLegend(0.16,0.79,0.40,0.88,NULL,"brNDC");
    		TLegend *leg_2 =new TLegend(0.79,0.60,0.91,0.88,NULL,"brNDC");
    		leg_2->AddEntry(gr_afb_fh_0," Bin 0 ", "P");
    		leg_2->AddEntry(gr_afb_fh_1," Bin 1 ", "P");
    		leg_2->AddEntry(gr_afb_fh_2," Bin 2 ", "P");
    		leg_2->AddEntry(gr_afb_fh_4," Bin 4 ", "P");
            leg_2->AddEntry(gr_afb_fh_6," Bin 6 ", "P");
    		leg_2->AddEntry(gr_afb_fh_7," Bin 7 ", "P");
    		leg_2->AddEntry(gr_afb_fh_8," Bin 8 ", "P");
    		leg_2->AddEntry(gr_afb_fh_9," Bin 9 ", "P");
    		leg_2->AddEntry(gr_afb_fh_10," Bin 10 ", "P");
    		leg_2->SetLineColor(1);
    		leg_2->SetFillColor(0);
    		leg_2->SetTextSize(0.02);
    		leg_2->Draw();					
    		c->Print(TString::Format("./plots/%s_afb_fh_bin%d.pdf",outfile,iBin));
    	}
	}
	
	c->Clear();
	t1->Clear();
/*	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = {-0.05,  -0.05};
	double conYAerr[2]= { 0.45,   0.45};
	double conYF[2]   = { 0.70,   0.70};
	double conYFerr[2]= { 0.80,   0.80};

	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9],yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9];
	for (int i = 0; i < 9; i++) {
		yfh[i]  = 0; yuerrfh[i] = 0; yderrfh[i] = 0;
		yafb[i] = 0; yuerrafb[i]= 0; yderrafb[i]= 0;
	}
////Checkout input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]      = readParam(ibin,"afb_L",0);
		yuerrafb[i]  = fabs(readParam(ibin,"afb_L",2));
		yderrafb[i]  = fabs(readParam(ibin,"afb_L",1));
		yfh[i]       = readParam(ibin,"fh_L",0);
		yuerrfh[i]   = fabs(readParam(ibin,"fh_L",2));
		yderrfh[i]   = fabs(readParam(ibin,"fh_L",1));
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("Afb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("Fh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i], yuerrfh[i], yderrfh[i]);
	}
//	Draw
//	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitleOffset(1.1,"Y");
	frame->SetTitle("");
	frame->Draw("2");	
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,1.5,"Y");
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	d_fh->SetMarkerColor(1);
	d_fh->SetMarkerStyle(20);
	d_fh->Draw("P");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();

//	TLegend *leg =new TLegend(0.69,0.69,0.90,0.86);
//	leg->AddEntry("d_fh"," Data ","lep");
//	leg->AddEntry("sys_fh"," Systematic error ","f");
//	leg->AddEntry("s_fh"," SM prediction ","f");
//	leg->Draw();
	
//	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.62,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Clear();
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("A_{FB}");
	frame->SetAxisRange(-0.5,0.4,"Y");
	frame->Draw();
	TGraphAsymmErrors *d_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	d_afb->SetMarkerColor(1);
	d_afb->SetMarkerStyle(20);
	d_afb->Draw("P");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2");
   
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.62,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
*/
}//}}}
void angular2D_priorSYS(int iBin, const char outfile[] = "angular2D_priorSYS", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1.1, 1.1);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -0.9, 0.9);
	
/*	// Create combinatorial background distribution
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",0.,-3.,3.);
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-3.,3.);
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-3.,3.);
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-3.,3.);
	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-2.,2.);
//	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
        case 1:
		case 9:
		    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
	        f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3,bkgCombL_c4));
		    break;
		case 2:
		case 10:
		    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
		    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
			break;
	    default:
		    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
	        //f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5,bkgCombL_c4));
		    break;
    }
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
//	RooChebychev f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
*/	// Create peak background distribution
	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}",  0.2, -0.1,  0.3);
	RooRealVar bkgGauss_mean2("bkgGauss_mean2","cos#theta_{l}",  0.2, -0.2,  0.5);
	RooRealVar bkgGauss_mean3("bkgGauss_mean3","cos#theta_{l}", -0.2, -0.9, -0.1);
//	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.2, -0.4, -0.1);
	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.2, -0.25, -0.16);
	RooRealVar bkgGauss_mean5("bkgGauss_mean5","cos#theta_{l}",  0.2,  0.0,  0.4);
	RooRealVar bkgGauss_mean6("bkgGauss_mean6","cos#theta_{l}",  0.1, -0.3,  0.3);
	RooRealVar bkgGauss_mean7("bkgGauss_mean7","cos#theta_{l}",  0.1,-0.20,  0.3);
	RooRealVar bkgGauss_mean8("bkgGauss_mean8","cos#theta_{l}",  0.1, -0.7,  0.5);
	
	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",  .50,  .01,  3);
	RooRealVar bkgGauss_sigma2("bkgGauss_sigma2","#sigma_{2}",  .20,  .01,  2);
	RooRealVar bkgGauss_sigma3("bkgGauss_sigma3","#sigma_{3}",  .80,  .01,  5);
//	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  .50,  .01,  1);
	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  .12,  .01,0.5);
	RooRealVar bkgGauss_sigma5("bkgGauss_sigma5","#sigma_{5}",  .40,  .01,  2);
	RooRealVar bkgGauss_sigma6("bkgGauss_sigma6","#sigma_{6}",  .12,  .10,0.3);
	RooRealVar bkgGauss_sigma7("bkgGauss_sigma7","#sigma_{7}",  .40,  .01,  2);
	RooRealVar bkgGauss_sigma8("bkgGauss_sigma8","#sigma_{8}",  .50,  .01,  3);
	
	RooGaussian f_bkgCombLGauss11("f_bkgCombLGauss11","f_bkgCombLGauss11", CosThetaL, bkgGauss_mean1, bkgGauss_sigma1);
	RooGaussian f_bkgCombLGauss22("f_bkgCombLGauss22","f_bkgCombLGauss22", CosThetaL, bkgGauss_mean2, bkgGauss_sigma2);
//	RooGaussian f_bkgCombLGauss33("f_bkgCombLGauss33","f_bkgCombLGauss33", CosThetaL, bkgGauss_mean3, bkgGauss_sigma3);
	RooGaussian f_bkgCombLGauss44("f_bkgCombLGauss44","f_bkgCombLGauss44", CosThetaL, bkgGauss_mean4, bkgGauss_sigma4);
	RooGaussian f_bkgCombLGauss55("f_bkgCombLGauss55","f_bkgCombLGauss55", CosThetaL, bkgGauss_mean5, bkgGauss_sigma5);
	RooGaussian f_bkgCombLGauss66("f_bkgCombLGauss66","f_bkgCombLGauss66", CosThetaL, bkgGauss_mean6, bkgGauss_sigma6);
	RooGaussian f_bkgCombLGauss77("f_bkgCombLGauss77","f_bkgCombLGauss77", CosThetaL, bkgGauss_mean7, bkgGauss_sigma7);
	RooGaussian f_bkgCombLGauss88("f_bkgCombLGauss88","f_bkgCombLGauss88", CosThetaL, bkgGauss_mean8, bkgGauss_sigma8);
	
//	RooGaussian f_bkgCombLGauss62("f_bkgCombLGauss62","f_bkgCombLGauss62", CosThetaL, bkgGauss_mean6, bkgGauss_sigma2);
//	RooGaussian f_bkgCombLGauss42("f_bkgCombLGauss42","f_bkgCombLGauss42", CosThetaL, bkgGauss_mean4, bkgGauss_sigma2);
	RooGaussian f_bkgCombLGauss31("f_bkgCombLGauss31","f_bkgCombLGauss31", CosThetaL, bkgGauss_mean3, bkgGauss_sigma1);
//	RooRealVar bkg_frac_G("bkg_frac_G","bkg_frac_G",1.,0.,1.);
//	RooAddPdf f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", RooArgList(f_bkgCombLGauss55,f_bkgCombLGauss66), bkg_frac_G);
//	
	RooRealVar bkg_frac("bkg_frac","bkg_frac",1.,0.,1.);

	RooAddPdf *f_bkgCombL = 0;
	switch (iBin) {
		case 9:
		    //f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss11, f_bkgCombL_P), bkg_frac);
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss11, f_bkgCombLGauss66), bkg_frac);
			break;
		case 0:
		    //f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss22, f_bkgCombL_P), bkg_frac);
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss22, f_bkgCombLGauss55), bkg_frac);
			break;
		case 10:
			//f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss77, f_bkgCombL_P), bkg_frac);
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss77, f_bkgCombLGauss55), bkg_frac);
			break;
		case 2:
			//f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss66, f_bkgCombL_P), bkg_frac);
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss66, f_bkgCombLGauss77), bkg_frac);
			break;
		case 1:
			//f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss77, f_bkgCombL_P), bkg_frac);
			//f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss55, f_bkgCombL_P), bkg_frac);
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss55, f_bkgCombLGauss66), bkg_frac);
			//f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombL_P), bkg_frac);
			break;
		case 6:
		    //f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss88, f_bkgCombL_P), bkg_frac);
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss88, f_bkgCombLGauss31), bkg_frac);
			break;
		case 4:
		    //f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss44, f_bkgCombL_P), bkg_frac);
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss44, f_bkgCombLGauss77), bkg_frac);
			break;
		case 7:
		case 8:
		    //f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss31, f_bkgCombL_P), bkg_frac);
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss31, f_bkgCombLGauss77), bkg_frac);
			break;
	}

	// Get data and apply unbinned fit
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.384 || Bmass < 5.174)",Q2range[iBin]),0);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.349 || Bmass < 5.209)",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.344 || Bmass < 5.214)",Q2range[iBin]),0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
/*	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	if ( iBin == 1 ) {
	    for (int i = 0; i < 2; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
*/
/*	if ( iBin == 2 || iBin == 9) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
*/
/*
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
    if ( iBin == 6 || iBin == 4) {
    	for (int i = 0; i < 2; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
*/

//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Extended(),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended(),Strategy(2),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));

//bin 1	
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
//bin 4,6	
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
//bin 0,2,7,8,9,10, 1, 4,6	
	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) return;
	
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
//	data->plotOn(framecosl,Binning(25)); 
//	data->plotOn(framecosl,Binning(22)); 
	data->plotOn(framecosl,Binning(20)); 
//	data->plotOn(framecosl,Binning(18)); 
	f_bkgCombL->plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
//	framecosl->SetAxisRange(-1.,1.,"X"); 
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.63,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	TPaveText* paveText_l = new TPaveText( 0.17, 0.76, 0.27, 0.86, "NDC" );
	paveText_l->SetBorderSize(0);
	paveText_l->SetFillColor(19);
    paveText_l->AddText(Form("bin %d ", iBin));
    paveText_l->Draw();
	
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d.pdf",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
	
	// Prepare datacard
	if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
/*        
        bkgCombL_c0.setConstant(kTRUE);
        bkgCombL_c1.setConstant(kTRUE);
        bkgCombL_c2.setConstant(kTRUE);
        bkgCombL_c3.setConstant(kTRUE);
*/        
        bkg_frac.setConstant(kTRUE);
//      wspace->import(f_bkgCombLGauss);
//      wspace->import(f_bkgCombL_P);
        wspace->import(*f_bkgCombL);
        wspace->writeToFile(TString::Format("%s/wspace_prior_bin%d.root",owspacepath.Data(),iBin),true);
		double val[3] = {0,0,0};
/*		val[0] = bkgCombL_c0.getVal();val[1] = bkgCombL_c0.getError();
		writeParam(iBin, "bkgCombL_c0", val);
		val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
		writeParam(iBin, "bkgCombL_c1", val);
		val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
		writeParam(iBin, "bkgCombL_c2", val);
		val[0] = 0; val[1] = 0;
		writeParam(iBin, "bkgCombL_c3", val);
	//	if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
		if (iBin == 2 || iBin ==9 || iBin == 10) {
			val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
			writeParam(iBin, "bkgCombL_c3", val);
	    }
		val[0] = 0; val[1] = 0;
		writeParam(iBin, "bkgCombL_c4", val);
		if (iBin ==9 ) {
			val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
			writeParam(iBin, "bkgCombL_c4", val);
	    }
*/		val[0] = 0; val[1] = 0;
		if ( iBin == 9) {
			val[0] = bkgGauss_mean1.getVal();val[1] = bkgGauss_mean1.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean6.getVal();val[1] = bkgGauss_mean6.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma6.getVal();val[1] = bkgGauss_sigma6.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 0 ) {
			val[0] = bkgGauss_mean2.getVal();val[1] = bkgGauss_mean2.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma2.getVal();val[1] = bkgGauss_sigma2.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean5.getVal();val[1] = bkgGauss_mean5.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma5.getVal();val[1] = bkgGauss_sigma5.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
		} else if (iBin == 10 ) {
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean5.getVal();val[1] = bkgGauss_mean5.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma5.getVal();val[1] = bkgGauss_sigma5.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 2 ) {
			val[0] = bkgGauss_mean6.getVal();val[1] = bkgGauss_mean6.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma6.getVal();val[1] = bkgGauss_sigma6.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 1) {
			val[0] = bkgGauss_mean5.getVal();val[1] = bkgGauss_mean5.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma5.getVal();val[1] = bkgGauss_sigma5.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean6.getVal();val[1] = bkgGauss_mean6.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma6.getVal();val[1] = bkgGauss_sigma6.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 4) {
			val[0] = bkgGauss_mean4.getVal();val[1] = bkgGauss_mean4.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma4.getVal();val[1] = bkgGauss_sigma4.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 6 ) {
			val[0] = bkgGauss_mean8.getVal();val[1] = bkgGauss_mean8.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma8.getVal();val[1] = bkgGauss_sigma8.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean3.getVal();val[1] = bkgGauss_mean3.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
	    } else if (iBin == 7 || iBin == 8) {
			val[0] = bkgGauss_mean3.getVal();val[1] = bkgGauss_mean3.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean_1", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma_1", val);
		}
		val[0] = bkg_frac.getVal();val[1] = bkg_frac.getError();
		writeParam(iBin, "bkg_frac", val);
    }
	
}//}}}


