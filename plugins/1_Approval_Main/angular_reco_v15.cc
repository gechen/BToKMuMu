// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
std::vector<double> angular_reco_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular_reco")
{//{{{
	setTDRStyle();
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1., 1.);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -0.97, 0.97);
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
//	if (iBin != 0 ) { 
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
	//	f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	//	f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	//	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
	//	f_rec_format = "( effP0 *exp(-0.5*(((CosThetaL-effP1)/effP2)**2)) ) * (effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3 )";
	//	f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	//	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6, effP7));
	//	f_rec_format = "( effP0 *exp(-0.5*(((CosThetaL-effP1)/effP2)**2)) ) * (effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**5 )";
	//	f_rec_format = "( effP0 *exp(-0.5*(((CosThetaL-effP1)/effP2)**2)) ) * (effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4 )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sig("f_sig", f_sigA_format,f_sigA_argset);
	RooExtendPdf  f("f","", f_sig, nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);    // 12-08
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(RooArgSet(afb, fh)),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Warnings(-1), PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Range(-0.9,0.9),Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), PrintEvalErrors(-1));
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
//	data->plotOn(framecosl,Binning(95)); 
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
		writeParam(iBin, "recofh", val, 3);
		val[1]=0; val[2]=0;
		val[0] = afb.getVal();val[1] = afb.getError(); val[2] = afb.getError();
	//	val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "recoafb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "recoFCN", val);
		printf("recoAfb[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recoafb",0), fabs(readParam(iBin,"recoafb",2)), fabs(readParam(iBin,"recoafb",1)));
		printf(" recoFh[%d]=%6.4f + %6.4f - %6.4f\n", iBin, readParam(iBin,"recofh",0), fabs(readParam(iBin,"recofh",2)), fabs(readParam(iBin,"recofh",1)));
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

void angular_reco(const char outfile[] = "angular_reco")
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
	double ygenfh[9], ygenderrfh[9], ygenuerrfh[9], ygenafb[9], ygenderrafb[9], ygenuerrafb[9]; 
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
	//	gen
		ygenafb[i]      = readParam(ibin,"genafb",0);
		ygenderrafb[i]  = fabs(readParam(ibin,"genafb",1));
		ygenuerrafb[i]  = fabs(readParam(ibin,"genafb",2));
		ygenfh[i]       = readParam(ibin,"genfh",0);
		ygenderrfh[i]   = fabs(readParam(ibin,"genfh",1));
		ygenuerrfh[i]   = fabs(readParam(ibin,"genfh",2));
		printf("genAfb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,ygenafb[i],ygenuerrafb[i],ygenderrafb[i]);
		printf("genFh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,ygenfh[i],ygenuerrfh[i],ygenderrfh[i]);
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
	TGraphAsymmErrors *gen_fh  = new TGraphAsymmErrors(7,x,ygenfh,xerr,xerr,ygenderrfh,ygenuerrfh);
	gen_fh->SetMarkerColor(3);
	gen_fh->SetMarkerStyle(24);
	gen_fh->SetFillColor(3);
	gen_fh->SetLineColor(3);
	gen_fh->SetFillStyle(3001);
	gen_fh->Draw("2");
	gen_fh->Draw("P");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();
	
	TLegend *leg =new TLegend(0.64,0.72,0.90,0.86);
	leg->AddEntry(g_fh," reco level"," lePf ");
	leg->AddEntry(gen_fh," gen level(unfiltered MC) "," lePf ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/angular_gen_reco_fh.pdf"));
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
	TGraphAsymmErrors *gen_afb = new TGraphAsymmErrors(7,x,ygenafb,xerr,xerr,ygenderrafb,ygenuerrafb);
	gen_afb->SetMarkerColor(3);
	gen_afb->SetMarkerStyle(24);
	gen_afb->SetFillColor(3);
	gen_afb->SetLineColor(3);
	gen_afb->SetFillStyle(3001);
	gen_afb->Draw("2");
	gen_afb->Draw("P");
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
	c->Print(TString::Format("./plots/angular_gen_reco_afb.pdf"));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
