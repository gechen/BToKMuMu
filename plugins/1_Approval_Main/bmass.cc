// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
///////////////////////////////////////////////////////////////////////////////////////////////////////
void bmass(int iBin, const char outfile[] = "bmass")
{//{{{
	bool test = false; 
	
	RooRealVar     DimuonPt("DimuonPt","DimuonPt", 0., 400.);
	RooRealVar     MupEta("MupEta","MupEta", -10., 10.);
	RooRealVar     MumEta("MumEta","MumEta", -10., 10.);
//	RooRealVar     Tk_4vec("Tk_4vec","Tk_4vec", -10., 10.);
	
    RooRealVar     Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.10, 5.60);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
	//TString::Format("(%s) && (%s) ",Q2range[iBin],mumuMassWindow[0]),0);	// data_v1
	TString::Format("(%s)",Q2range[iBin]),0);	// data_v1
	//RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr, DimuonPt, MupEta, MumEta),
	//TString::Format("(%s) && (DimuonPt >7.0) && (abs(MupEta)<=2.0) && (abs(MumEta)<=2.0)",Q2range[iBin]),0);	// data_v1
	//RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr, DimuonPt, MupEta, MumEta, Tk_4vec),
	//TString::Format("(%s) && (DimuonPt >7.0) && (fabs(MupEta)<=2.0) && (fabs(MumEta)<=2.0) && (fabs(Tk_4vec.Eta())<=2.4)",Q2range[iBin]),0);	// data_v1

	data->Print();
	
//	Create model and dataset
//	-------------------------------------------------------------------------
//	Gaussian signal 
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.25, 5.32);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
//	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.35);
//	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.005, 0.00, 0.02);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.030, 0.01, 0.1);
//	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.100, 0.06, 0.2);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
	
//	Build Chebychev polynomial p.d.f.  
	RooRealVar     a0("a0", "constant", 0.4, -1, 1);
	RooRealVar     a1("a1", "linear", 0.1, -1, 1);
	RooRealVar     a2("a2", "quadratic", 0.1, -1, 1);
	RooChebychev   bkg("bkg", "Background", Bmass, RooArgSet(a0, a1, a2));
	
//	Construct signal+background PDF
	RooRealVar     nsig("nsig", "number of signal events", 4648, 0, 1E8); 
	RooRealVar     nbkg("nbkg", "number of background events", 21472, 0, 1E8);
	RooAddPdf      model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig));
	
//	Print structure of composite p.d.f.
	model.Print("t");
	
//	Fit model to data, save fitresult 
//	------------------------------------------------------------------------
	RooFitResult* fitres; 
	if (! test) {
		fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE));
		fitres->Print("v"); 
	}
	
//	Plot model 
//	---------------------------------------------------------
	TString title = "B^{+/-} mass";
//	int nbins = 25; 
	int nbins = 20; 
	RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
	data->plotOn(frame);
	model.plotOn(frame, LineColor(1));
	
//	Overlay the background component of model with a dashed line
	model.plotOn(frame,Components("bkg"), LineStyle(kDashed), LineColor(2));
//	Overlay the signal component of model with a blue line
	model.plotOn(frame,Components("sig"), LineStyle(1), LineColor(4));
	
//	Draw the frame on the canvas
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
//	double chi2Val=0;
//	chi2Val = model.GetChisquare();
//	latex->DrawLatexNDC(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
	
	TPaveText* paveText = new TPaveText( 0.62, 0.70, 0.89, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("nbkg   = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
	paveText->AddText(Form("mean   = %.5f #pm %.5f ", mean.getVal(), mean.getError())); 
	paveText->AddText(Form("sigma1 = %.5f #pm %.5f ", sigma1.getVal(), sigma1.getError())); 
	paveText->AddText(Form("sigma2 = %.5f #pm %.5f ", sigma2.getVal(), sigma2.getError())); 
	paveText->AddText(Form("frac   = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->Draw(); 
	
	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
//	delete paveText; 
	delete c;

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void bmass_JpsiKGEN(int iBin, const char outfile[] = "bmass_JpsiKGEN")
{//{{{
	bool test = false; 
	
	RooRealVar     Bmass("Bmass", "J/#psi K mass(GeV/c^{2})", 5.10, 5.60);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
	TString::Format("(%s)",Q2range[iBin]),0);	// JpsiK
//	TString::Format("(%s) && (%s) && (%s)",Q2range[iBin],mumuMassWindow[3],bmassWindow[0]),0);	// JpsiK

	data->Print();
	
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.23, 5.32);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.35);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
//	Construct signal PDF
	RooRealVar     nsig("nsig", "number of signal events", 4648, 0, 1E8); 
	RooAddPdf      model("model", "g", RooArgList(sig), RooArgList(nsig));
//	Print structure of composite p.d.f.
	model.Print("t");
//	Fit model to data, save fitresult 
	RooFitResult* fitres; 
	if (! test) {
		fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE));
		fitres->Print("v"); 
	}
//	Plot model 
	TString title = "B^{+/-} mass of J/#psi K MC Control Channel ";
//	TString title = "B^{+/-} mass of J/#psi K Data Control Channel ";
	int nbins = 25; 
	RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
	data->plotOn(frame);
	model.plotOn(frame, LineColor(1));
	model.plotOn(frame,Components("sig"), LineStyle(1), LineColor(4));
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
	TPaveText* paveText = new TPaveText( 0.62, 0.70, 0.89, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("mean   = %.5f #pm %.5f ", mean.getVal(), mean.getError())); 
	paveText->AddText(Form("sigma1 = %.5f #pm %.5f ", sigma1.getVal(), sigma1.getError())); 
	paveText->AddText(Form("sigma2 = %.5f #pm %.5f ", sigma2.getVal(), sigma2.getError())); 
	paveText->AddText(Form("frac   = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->Draw(); 
	
//	c->Print(TString::Format("./plots/%s_bin%d_JpsiKMC.pdf",outfile,iBin));
//	c->Print(TString::Format("./plots/%s_bin%d_Data.pdf",outfile,iBin));
	c->Print(TString::Format("./plots/%s_bin%d_Unfilter_SignalMC.pdf",outfile,iBin));
	delete paveText; 
	delete c;

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void bmass_JpsiK(int iBin, const char outfile[] = "bmass_JpsiK")
{//{{{
	bool test = false; 
	
	RooRealVar     Bmass("Bmass", "J/#psi K mass(GeV/c^{2})", 5.10, 5.60);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
	TString::Format("(%s)",Q2range[iBin]),0);	// JpsiK
//	TString::Format("(%s) && (%s) && (%s)",Q2range[iBin],mumuMassWindow[3],bmassWindow[0]),0);	// JpsiK

	data->Print();
	
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.23, 5.32);
//	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.08);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.35);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
//	Construct signal PDF
	RooRealVar     nsig("nsig", "number of signal events", 4648, 0, 1E8); 
	RooAddPdf      model("model", "g", RooArgList(sig), RooArgList(nsig));
//	Print structure of composite p.d.f.
	model.Print("t");
//	Fit model to data, save fitresult 
	RooFitResult* fitres; 
	if (! test) {
		fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE));
		fitres->Print("v"); 
	}
//	Plot model 
	TString title = "B^{+/-} mass of J/#psi K MC Control Channel ";
//	TString title = "B^{+/-} mass of J/#psi K Data Control Channel ";
	int nbins = 25; 
	RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
	data->plotOn(frame);
	model.plotOn(frame, LineColor(1));
	model.plotOn(frame,Components("sig"), LineStyle(1), LineColor(4));
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
	TPaveText* paveText = new TPaveText( 0.62, 0.70, 0.89, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("mean   = %.5f #pm %.5f ", mean.getVal(), mean.getError())); 
	paveText->AddText(Form("sigma1 = %.5f #pm %.5f ", sigma1.getVal(), sigma1.getError())); 
	paveText->AddText(Form("sigma2 = %.5f #pm %.5f ", sigma2.getVal(), sigma2.getError())); 
	paveText->AddText(Form("frac   = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->Draw(); 
	
	c->Print(TString::Format("./plots/%s_bin%d_JpsiKMC.pdf",outfile,iBin));
//	c->Print(TString::Format("./plots/%s_bin%d_Data.pdf",outfile,iBin));
//	c->Print(TString::Format("./plots/%s_bin%d_SignalMC.pdf",outfile,iBin));
	delete paveText; 
	delete c;

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void bmass_Psi2SK(int iBin, const char outfile[] = "bmass_Psi2SK")
{//{{{
	bool test = false; 
	
	RooRealVar     Bmass("Bmass", "#psi' K mass(GeV/c^{2})", 5.10, 5.60);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
	TString::Format("(%s)",Q2range[iBin]),0);	// Psi2SK
//	TString::Format("(%s) && (%s) && (%s)",Q2range[iBin],mumuMassWindow[4],bmassWindow[0]),0);	// Psi2SK

	data->Print();
	
	RooRealVar     mean("mean","mean of gaussians", 5.27925, 5.23, 5.32);
	RooRealVar     sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05);
	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.03, 0.35);
//	RooRealVar     sigma2("sigma2","width of Gaussian2", 0.065, 0.05, 0.35);
	RooRealVar     sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
	RooGaussian    sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1);
	RooGaussian    sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2);
	RooAddPdf      sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));
//	Construct signal PDF
	RooRealVar     nsig("nsig", "number of signal events", 4648, 0, 1E8); 
	RooAddPdf      model("model", "g", RooArgList(sig), RooArgList(nsig));
//	Print structure of composite p.d.f.
	model.Print("t");
//	Fit model to data, save fitresult 
	RooFitResult* fitres; 
	if (! test) {
		fitres = model.fitTo(*data, Extended(kTRUE), Save(kTRUE));
		fitres->Print("v"); 
	}
//	Plot model 
	TString title = "B^{+/-} mass of #psi' K MC Control Channel ";
//	TString title = "B^{+/-} mass of #psi' K Data Control Channel ";
	int nbins = 25; 
	RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
	data->plotOn(frame);
	model.plotOn(frame, LineColor(1));
	model.plotOn(frame,Components("sig"), LineStyle(1), LineColor(4));
	TCanvas *c = new TCanvas("c", "c", 800, 600); 
	set_root_style(); 
	c->UseCurrentStyle();
	
	gPad->SetLeftMargin(0.15);
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->Draw();
	
	TPaveText* paveText = new TPaveText( 0.62, 0.70, 0.89, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
	paveText->AddText(Form("mean   = %.5f #pm %.5f ", mean.getVal(), mean.getError())); 
	paveText->AddText(Form("sigma1 = %.5f #pm %.5f ", sigma1.getVal(), sigma1.getError())); 
	paveText->AddText(Form("sigma2 = %.5f #pm %.5f ", sigma2.getVal(), sigma2.getError())); 
	paveText->AddText(Form("frac   = %.5f #pm %.5f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->Draw(); 
	
	c->Print(TString::Format("./plots/%s_bin%d_Psi2SKMC.pdf",outfile,iBin));
//	c->Print(TString::Format("./plots/%s_bin%d_Data.pdf",outfile,iBin));
//	c->Print(TString::Format("./plots/%s_bin%d_SignalMC.pdf",outfile,iBin));
	delete paveText; 
	delete c;

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
