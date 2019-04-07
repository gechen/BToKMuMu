// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_1a_Sm(int iBin, const char outfile[] = "angular2D_1a_Sm", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10, 5.60);
//	RooRealVar Bmass("Bmass","M_{K^{+/-}#Mu#Mu}",5.27925-0.28, 5.27925+0.28);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	
	// Create parameters and PDFs
	   // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K#Mu#Mu}",5.27925,5.24,5.31);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",.025,0.001,0.06);
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.075,0.060,0.15);
//	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",.065,0.01,0.15);
	RooRealVar sigM_frac("sigM_frac","sigM_frac",0.5,0.,1.);
	
	// Create signal distribution
	   // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", f_sigMGauss1, f_sigMGauss2, sigM_frac);
   
	RooRealVar nsig("nsig","nsig",1E5,0,1E8);
	RooAddPdf f("f", "f",RooArgList(f_sigM),RooArgList(nsig));
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE),Strategy(2),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(kTRUE),PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) return;

	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	RooPlot* frame = Bmass.frame(); 
	data->plotOn(frame,Binning(20)); 
	f.plotOn(frame,LineColor(1)); 
	f.plotOn(frame,Components(f_sigM),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
	f.plotOn(frame,Components(f_sigM),LineStyle(2),LineColor(4),LineWidth(2));

	frame->SetTitle("");
	frame->SetMinimum(0);
	frame->SetTitleOffset(1.1,"Y");
	frame->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.51,.90,TString::Format("signal MC: 3296.81 fb^{-1}(8TeV)"));
	TPaveText* paveText_l = new TPaveText( 0.17, 0.76, 0.27, 0.86, "NDC" );
	paveText_l->SetBorderSize(0);
	paveText_l->SetFillColor(19);
    paveText_l->AddText(Form("bin %d ", iBin));
    paveText_l->Draw();
	
	TPaveText* paveText = new TPaveText( 0.65, 0.70, 0.90, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
	paveText->AddText(Form("   nsig   = %.0f #pm %.0f ", nsig.getVal(), nsig.getError())); 
//	paveText->AddText(Form(" nbkg   = %.0f #pm %.0f ", nbkg.getVal(), nbkg.getError())); 
//	paveText->AddText(Form(" sigmean   = %.3f #pm %.3f ", sigGauss_mean.getVal(), sigGauss_mean.getError())); 
//	paveText->AddText(Form("sig_sigma1 = %.3f #pm %.3f ", sigGauss1_sigma.getVal(), sigGauss1_sigma.getError())); 
//	paveText->AddText(Form("sig_sigma2 = %.3f #pm %.3f ", sigGauss2_sigma.getVal(), sigGauss2_sigma.getError())); 
//	paveText->AddText(Form("    frac   = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
	paveText->AddText(Form("#Chi^{2}   = %.2f  ", frame->chiSquare())); 
	paveText->Draw(); 
	
	c->Update();
	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	// clear
	delete t1;
	delete c;
	delete data;
	// Prepare datacard
	if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        nsig.setConstant(kTRUE);
        sigGauss_mean.setConstant(kTRUE);
        sigGauss1_sigma.setConstant(kTRUE);
        sigGauss2_sigma.setConstant(kTRUE);
        sigM_frac.setConstant(kTRUE);
        wspace->import(nsig);
        wspace->import(f_sigM);
        wspace->import(f);
        wspace->writeToFile(TString::Format("%s/wspace_Sm_bin%d.root",owspacepath.Data(),iBin),true);
		double val[3]={0,0,0};
		writeParam(iBin, "iBin", new double((double)iBin), 1);
		if (is7TeVCheck){
			writeParam(iBin, "mode", new double(2011), 1);
		}else{
			writeParam(iBin, "mode", new double(2012), 1);
		}
		val[0]=sigGauss1_sigma.getVal();val[1]=sigGauss1_sigma.getError();
		writeParam(iBin, "sigGauss1_sigma", val);
		val[0]=sigGauss2_sigma.getVal();val[1]=sigGauss2_sigma.getError();
		writeParam(iBin, "sigGauss2_sigma", val);
		val[0]=sigM_frac.getVal();val[1]=sigM_frac.getError();
		writeParam(iBin, "sigM_frac", val);
	}
}//}}}

void angular2D_1b_YpPm(int iBin, const char outfile[] = "angular2D_1b_YpPm", bool keepParam = true)
{//{{{
	setTDRStyle();
	bool mctype = false;
	TString func = outfile;
	if (func == "angular2D_1b_YpPm_Jpsi") {
		if (iBin != 10 && iBin != 4 && iBin != 2){    
	//	if (iBin != 10 && iBin != 2){    
			mctype = true;
		}
	}else if (func == "angular2D_1b_YpPm_Psi") {
		if (iBin != 10 && iBin != 4 && iBin != 6){
	//	if (iBin != 10 && iBin != 4 ){
			mctype = true;
		}
	}else if (func == "angular2D_1b_YpPm_JX") {
		if (iBin != 10 && iBin != 4 && iBin != 6 && iBin != 2){
	//	if (iBin != 10 && iBin != 2){
			mctype = true;
		}
	}
/*	if (keepParam && mctype){
		double val[3]={-1,0,0};
		writeParam(iBin, TString::Format("bkgPeakM_c_%s",outfile), val);
		val[0]=0;
		writeParam(iBin, TString::Format("nbkgPeak_%s",outfile), val);
		return;
	}
*/	// Fit to control channel simulations by YpPm to determine Yp,Pm.
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	
	// Create peak background distribution
	RooRealVar bkgPeakM_c("bkgPeakM_c","c",-0.5,-20.,1.);
	RooRealVar bkgoffset("bkgoffset","bkgoffset",-5.);
	RooAddition bkgBmass_offset("bkgBmass_offset","bkgBmass_offset",RooArgList(Bmass,bkgoffset));
	RooExponential f_bkgPeakM_P("f_bkgPeakM_P","f_bkgPeakM_P",bkgBmass_offset,bkgPeakM_c);// exponential decay

	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",5E2,1,3E3);
	RooExtendPdf *f = new RooExtendPdf("f","f",f_bkgPeakM_P,nbkgPeak);
	
	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass),Q2range[iBin],0);
	RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended(),Strategy(2),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Strategy(2),Warnings(-1), Minos(kTRUE),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f->fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended());
	f_fitresult->Print();
//	if (f_fitresult->status() != 0) return;
	
	// Draw the frame on the canvas
	TCanvas* c = new TCanvas("c");
	RooPlot* frame = Bmass.frame(); 
	data->plotOn(frame,Binning(20)); 
	f->plotOn(frame); 
	frame->SetTitle("");
	frame->SetMinimum(0);
	frame->SetTitleOffset(1.1,"Y");
	frame->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	if (func == "angular2D_1b_YpPm_Jpsi") {
		t1->DrawLatex(.56,.90,TString::Format("J/#Psi K MC: 18.58 fb^{-1}(8TeV)"));
	}else if (func == "angular2D_1b_YpPm_Psi") {
		t1->DrawLatex(.56,.90,TString::Format("#Psi(2S) K MC: 212.50 fb^{-1}(8TeV)"));
	}else if (func == "angular2D_1b_YpPm_JX") {
		t1->DrawLatex(.56,.90,TString::Format("#Psi(#mu^{+}#mu^{-}) X MC: 9.82 fb^{-1}(8TeV)"));
	}
	TPaveText* paveText_l = new TPaveText( 0.17, 0.76, 0.27, 0.86, "NDC" );
	paveText_l->SetBorderSize(0);
	paveText_l->SetFillColor(19);
    paveText_l->AddText(Form("bin %d ", iBin));
    paveText_l->Draw();
	c->Update();
	c->Print(TString::Format("./plots/%s_bin%d.pdf",outfile,iBin));
	
	// clear
	delete t1;
	delete c;
	delete data;
/*
	if (keepParam){
		double val[3]={0,0,0};
		val[0] = bkgPeakM_c.getVal();val[1] = bkgPeakM_c.getError();
		writeParam(iBin, TString::Format("bkgPeakM_c_%s",outfile), val);
		val[0]=0;
		val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
		writeParam(iBin, TString::Format("nbkgPeak_%s",outfile), val);
	}
*/	
}//}}}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_2a_PkPl_Jpsi(int iBin, const char outfile[] = "angular2D_2a_PkPl_Jpsi", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Gaussian constraint on yields and mass is needed.
	bool mctype = false;
	TString func = outfile;
	const char read1b[] = "angular2D_1b_YpPm_Jpsi";     
	if (iBin != 10 && iBin != 4 && iBin != 2){
//	if (iBin != 10 && iBin != 2){
		mctype = true;
	}
/*	if (keepParam && mctype){
		double val[3]={0,0,0};
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
		return;
	}
*/	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.0, 1.0);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.1, 1.1);
	RooArgSet f_bkgPeakL_argset;
	RooRealVar bkgPeakL_c0("bkgPeakL_c0","c0",0.1,-5,5);
	RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0.1,-5,5);
	RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0.2,-5,5);
	RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0.2,-5,5);
//	RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c0,bkgPeakL_c1,bkgPeakL_c2));
	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c3));
	//	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c4));
	RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
//	RooChebychev f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
//	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",100,0.,2E3);
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",400,0.,3E3);
	RooExtendPdf f_bkgPeakL_ext("f_bkgPeakL_ext","f_bkgPeakL_ext",f_bkgPeakL,nbkgPeak);
	
	// Gaussian Constraint
	RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b), 0)),RooConst(readParam(iBin, TString::Format("nbkgPeak_%s",read1b), 1) ) );
	// Get data
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	for (int i = 0; i < 2; i++) {
		CosThetaL = 1 + (i+1) * .02;
		d->add(RooArgSet(CosThetaL));
		CosThetaL = -1 - (i+1) * .02;
		d->add(RooArgSet(CosThetaL));
	}
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
	
//	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
//	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*data,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak),Strategy(2),Warnings(-1),PrintEvalErrors(-1));
	f_fitresult->Print();
//	if (f_fitresult->status() != 0) return;
	
//	Draw CosThetaL
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
//	data->plotOn(framecosl,Binning(22)); 
	data->plotOn(framecosl,Binning(20)); 
//	f_bkgPeakL.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.56,.90,TString::Format("J/#Psi K MC: 18.58 fb^{-1}(8TeV)"));
	
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
/*
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgPeakL_c0.getVal();val[1] = bkgPeakL_c0.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
	}
*/	
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_2a_PkPl_Psi(int iBin, const char outfile[] = "angular2D_2a_PkPl_Psi", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Gaussian constraint on yields and mass is needed.
	bool mctype = false;
	TString func = outfile;
	const char read1b[] = "angular2D_1b_YpPm_Psi";
	if (iBin != 10 && iBin != 4 && iBin != 6){
//	if (iBin != 10 && iBin != 4 ){
		mctype = true;
	}
/*	if (keepParam && mctype){
		double val[3]={0,0,0};
		writeParam(iBin, TString::Format("bkgPeakL_mean_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_sigma_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_frac_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
		return;
	}
*/	
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.0, 1.0);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.1, 1.1);
	RooArgSet f_bkgPeakL_argset;
	RooRealVar bkgPeakL_c0("bkgPeakL_c0","c0",0.2,-5.,5.);
	RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0.1,-5.,5.);
	RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0.2,-5.,5.);
	RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0.1,-5.,5.);
//	RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0.2,-5.,5.);
	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c0,bkgPeakL_c1,bkgPeakL_c2));
	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c3));
//	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c4));
	RooPolynomial f_bkgPeakL_P("f_bkgPeakL_P","f_bkgPeakL_P",CosThetaL,f_bkgPeakL_argset);
	
	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}", -0.01, -0.11, 0.09);
	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",   .02,  .0, 0.15);
	RooGaussian f_bkgPeakLGauss11("f_bkgPeakLGauss11","f_bkgPeakLGauss11", CosThetaL, bkgGauss_mean1, bkgGauss_sigma1);
	RooRealVar bkg_frac_G("bkg_frac_G","bkg_frac_G",1.,0.,1.);
	RooAddPdf f_bkgPeakL("f_bkgPeakL","f_bkgPeakL", RooArgList(f_bkgPeakLGauss11,f_bkgPeakL_P), bkg_frac_G);
	
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",450,0.,3E3);
	RooExtendPdf f_bkgPeakL_ext("f_bkgPeakL_ext","f_bkgPeakL_ext",f_bkgPeakL,nbkgPeak);
	
	// Gaussian Constraint
	RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b), 0)),RooConst(readParam(iBin, TString::Format("nbkgPeak_%s",read1b), 1) ) );
	// Get data
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	for (int i = 0; i < 4; i++) {
	//	CosThetaL = 1 + (i+1) * .02;
		CosThetaL = 1 + 1 * .05;
		d->add(RooArgSet(CosThetaL));
	//	CosThetaL = -1 - (i+1) * .02;
		CosThetaL = -1 - 1 * .05;
		d->add(RooArgSet(CosThetaL));
	}
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
//	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"));  // 2015-09-11
	
	f_fitresult->Print();
//	if (f_fitresult->status() != 0) return;
	
//	Draw CosThetaL
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(20)); 
//	data->plotOn(framecosl,Binning(22)); 
//	f_bkgPeakL.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.56,.90,TString::Format("#Psi(2S) K MC: 212.50 fb^{-1}(8TeV)"));
	
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
/*	
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgGauss_mean1.getVal();val[1] = bkgGauss_mean1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_mean_%s",outfile), val);
		val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_sigma_%s",outfile), val);
		val[0] = bkg_frac_G.getVal();val[1] = bkg_frac_G.getError();
		writeParam(iBin, TString::Format("bkgPeakL_frac_%s",outfile), val);

		val[0] = bkgPeakL_c0.getVal();val[1] = bkgPeakL_c0.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
		val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
	}	
*/	
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_2a_PkPl_JX(int iBin, const char outfile[] = "angular2D_2a_PkPl_JX", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Gaussian constraint on yields and mass is needed.
	bool mctype = false;
	TString func = outfile;
	const char read1b[] = "angular2D_1b_YpPm_JX";     
	if (iBin != 10 && iBin != 4 && iBin != 2 && iBin != 6){
//	if (iBin != 10 && iBin != 2 ){
		mctype = true;
	}
/*	if (keepParam && mctype){
		double val[3]={0,0,0};
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
	//	writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
		return;
	}
*/	
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.0, 1.0);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}^{reco}", -1.1, 1.1);
	RooArgSet f_bkgPeakL_argset;
	RooRealVar bkgPeakL_c0("bkgPeakL_c0","c0",0.1,-5,5);
	RooRealVar bkgPeakL_c1("bkgPeakL_c1","c1",0.1,-5,5);
	RooRealVar bkgPeakL_c2("bkgPeakL_c2","c2",0.2,-5,5);
//	RooRealVar bkgPeakL_c3("bkgPeakL_c3","c3",0,-5,5);
//	RooRealVar bkgPeakL_c4("bkgPeakL_c4","c4",0,-5,5);
	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c0,bkgPeakL_c1,bkgPeakL_c2));
	//	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c3));
	//	f_bkgPeakL_argset.add(RooArgSet(bkgPeakL_c4));
	
	RooPolynomial f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
//	RooChebychev f_bkgPeakL("f_bkgPeakL","f_bkgPeakL",CosThetaL,f_bkgPeakL_argset);
	RooRealVar nbkgPeak("nbkgPeak","nbkgPeak",200,0.,2E3);
	RooExtendPdf f_bkgPeakL_ext("f_bkgPeakL_ext","f_bkgPeakL_ext",f_bkgPeakL,nbkgPeak);
	
	// Gaussian Constraint
	RooGaussian gaus_nbkgPeak("gaus_nbkgPeak","gaus_nbkgPeak",nbkgPeak,RooConst(readParam(iBin,TString::Format("nbkgPeak_%s",read1b), 0)),RooConst(readParam(iBin, TString::Format("nbkgPeak_%s",read1b), 1) ) );
	// Get data
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaL,Q2),Q2range[iBin],0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	for (int i = 0; i < 2; i++) {
		CosThetaL = 1 + (i+1) * .02;
		d->add(RooArgSet(CosThetaL));
		CosThetaL = -1 - (i+1) * .02;
		d->add(RooArgSet(CosThetaL));
	}
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
	RooFitResult *f_fitresult = f_bkgPeakL_ext.fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"),ExternalConstraints(gaus_nbkgPeak));
	
	f_fitresult->Print();
//	if (f_fitresult->status() != 0) return;
	
//	Draw CosThetaL
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,Binning(20)); 
//	data->plotOn(framecosl,Binning(22)); 
//	f_bkgPeakL.plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.56,.90,TString::Format("#Psi(#mu^{+}#mu^{-}) X MC: 9.81 fb^{-1}(8TeV)"));
	
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
/*	
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgPeakL_c0.getVal();val[1] = bkgPeakL_c0.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c0_%s",outfile), val);
		val[0] = bkgPeakL_c1.getVal();val[1] = bkgPeakL_c1.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c1_%s",outfile), val);
		val[0] = bkgPeakL_c2.getVal();val[1] = bkgPeakL_c2.getError();
		writeParam(iBin, TString::Format("bkgPeakL_c2_%s",outfile), val);
	//	val[0] = bkgPeakL_c3.getVal();val[1] = bkgPeakL_c3.getError();
	//	writeParam(iBin, TString::Format("bkgPeakL_c3_%s",outfile), val);
	//	val[0] = bkgPeakL_c4.getVal();val[1] = bkgPeakL_c4.getError();
	//	writeParam(iBin, TString::Format("bkgPeakL_c4_%s",outfile), val);
	}	
*/
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

