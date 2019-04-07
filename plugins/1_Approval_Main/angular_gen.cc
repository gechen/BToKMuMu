// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular_gen_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular_gen")
{//{{{
	setTDRStyle();
	RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{l}^{gen}", -1., 1.);
	RooRealVar genQ2("genQ2","q^{2}",1.0,22.);
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................
	
	RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
	
	RooArgSet f_sigA_argset(genCosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_ang_format;
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-genCosThetaL*genCosThetaL) + 0.5*fh + afb*genCosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-genCosThetaL*genCosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*genCosThetaL )";
	}
	f_sigA_argset.add(RooArgSet(f_ang_format));
	RooGenericPdf f_sig("f_sig", f_ang_format,f_sigA_argset);
	RooExtendPdf f("f","",f_sig,nsig);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaL,genQ2),genQ2range[iBin],0);
//	RooFitResult *f_fitresult = f.fitTo(*data, Extended(kTRUE), Save(kTRUE), Minimizer("Minuit"), Strategy(2),Warnings(-1), PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f.fitTo(*data, Extended(kTRUE), Save(kTRUE), Minimizer("Minuit"), Strategy(2),Warnings(-1), Minos(RooArgSet(afb, fh)), PrintEvalErrors(-1));
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
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.25);
	framecosl->Draw();
////////////////////////////////////// scan   //////////////////////////////////.........................
//	t1->DrawLatex(.30,.83,TString::Format("F_{H}  =%4.4f ", 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi()  ));
//	t1->DrawLatex(.60,.83,TString::Format("A_{FB} =%4.4f ", (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() )   ));
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
//	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%6.4f #pm %6.4f ", fh.getVal(),fh.getError()));
//	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%6.4f #pm %6.4f ", afb.getVal(),afb.getError()));
	t1->DrawLatex(.19,.83,TString::Format("F_{H}  =%4.4f + %4.4f %4.4f",fh.getVal(),fh.getAsymErrorHi(),fh.getAsymErrorLo()));
	t1->DrawLatex(.45,.78,TString::Format("A_{FB} =%4.4f + %4.4f %4.4f",afb.getVal(),afb.getAsymErrorHi(),afb.getAsymErrorLo()));
////////////////////////////////////// refit  //////////////////////////////////.........................
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.64,.90,TString::Format("private unfiltered MC"));
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
	double val[4]={0,0,0,0};
	if (Index == -1) {
	//	val[0] = fh.getVal();val[1] = fh.getError();val[2] = fh.getError();
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, "genfh", val, 3);
	//	val[0] = afb.getVal();val[1] = afb.getError();val[2] = afb.getError();
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "genafb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "genFCN", val);
		printf("genAfb[%d]=%6.4f + %6.4f - %6.4f \n", iBin, readParam(iBin,"genafb",0), fabs(readParam(iBin,"genafb",2)), fabs(readParam(iBin,"genafb",1)));
		printf("genFh[%d]=%6.4f + %6.4f - %6.4f \n", iBin, readParam(iBin,"genfh",0), fabs(readParam(iBin,"genfh",2)), fabs(readParam(iBin,"genfh",1)));
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

void angular_gen( const char outfile[] = "angular_gen")
{//{{{
	setTDRStyle();
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.09,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9], yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9]; 
//	Check input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]     = readParam(ibin,"genafb",0);
		yderrafb[i]  = fabs(readParam(ibin,"genafb",1));
		yuerrafb[i]  = fabs(readParam(ibin,"genafb",2));
		yfh[i]      = readParam(ibin,"genfh",0);
		yderrfh[i]   = fabs(readParam(ibin,"genfh",1));
		yuerrfh[i]   = fabs(readParam(ibin,"genfh",2));
		printf("genAfb[%d] =%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("genFh[%d]  =%6.4f + %6.4f - %6.4f\n",ibin,yfh[i],yuerrfh[i],yderrfh[i]);
	}
//	plotting
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->SetTitle("");
//	Fh
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,0.5,"Y");  
	frame->Draw();
	TGraphAsymmErrors *g_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	g_fh->SetMarkerColor(3);
	g_fh->SetMarkerStyle(24);
	g_fh->SetFillColor(3);
	g_fh->SetFillStyle(3001);
	g_fh->Draw("2");
	g_fh->Draw("P");
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.64,.90,TString::Format("private unfiltered MC"));
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Clear();
//	Afb	
	frame->SetYTitle("A_{FB}");
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
//	frame->SetAxisRange(-0.02,0.02,"Y");
	frame->SetAxisRange(-0.2,0.2,"Y");  
	frame->Draw();
	TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	g_afb->SetMarkerColor(3);
	g_afb->SetMarkerStyle(24);
	g_afb->SetFillColor(3);
	g_afb->SetFillStyle(3001);
	g_afb->Draw("2");
	g_afb->Draw("P");
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.64,.90,TString::Format("private unfiltered MC"));
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
	c->Clear();
	c->Close();
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
