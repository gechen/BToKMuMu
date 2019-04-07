// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------

void angular2D_prior(int iBin, const char outfile[] = "angular2D_prior", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1.1, 1.1);
	
	// Create combinatorial background distribution
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",0.,-3.,3.);
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-3.,3.);
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-3.,3.);
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
    if (iBin == 2 || iBin == 7 || iBin == 8 || iBin == 9) {
        f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
    }
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
//	RooChebychev f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	// Create peak background distribution
	RooRealVar bkgGauss_mean0("bkgGauss_mean0","cos#theta_{l}",  0.2, -0.2,  0.5);
	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}",  0.2, 0.01,  0.3);
	RooRealVar bkgGauss_mean2("bkgGauss_mean2","cos#theta_{l}",  0.01, -0.2,  0.2);
	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.2, -0.3, -0.1);
	RooRealVar bkgGauss_mean6("bkgGauss_mean6","cos#theta_{l}", -0.6, -0.9, -0.5);
	RooRealVar bkgGauss_mean7("bkgGauss_mean7","cos#theta_{l}", -0.7, -0.9, -0.5);
	RooRealVar bkgGauss_mean8("bkgGauss_mean8","cos#theta_{l}",  0.3,  0.1,  0.35);
	RooRealVar bkgGauss_mean9("bkgGauss_mean9","cos#theta_{l}",  0.2,  0.08,  0.25);
	RooRealVar bkgGauss_mean10("bkgGauss_mean10","cos#theta_{l}",  0.1,  0.0,  0.2);
	
	RooRealVar bkgGauss_sigma0("bkgGauss_sigma0","#sigma_{0}",  .20,  .01,  2);
	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",  .20,  .01,0.4);
	RooRealVar bkgGauss_sigma2("bkgGauss_sigma2","#sigma_{2}",  .12,  .10,0.3);
	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  .20,  .01,0.4);
	RooRealVar bkgGauss_sigma6("bkgGauss_sigma6","#sigma_{6}",  0.6,  .1,   1);
	RooRealVar bkgGauss_sigma7("bkgGauss_sigma7","#sigma_{7}",  .20,  .1,  0.5);
	RooRealVar bkgGauss_sigma8("bkgGauss_sigma8","#sigma_{8}",  .10,  .06, 0.3);
	RooRealVar bkgGauss_sigma9("bkgGauss_sigma9","#sigma_{9}",  .10,  .06,  0.2);
	RooRealVar bkgGauss_sigma10("bkgGauss_sigma10","#sigma_{10}",  .20,  .1,   1);

	
	RooGaussian f_bkgCombLGauss00("f_bkgCombLGauss00","f_bkgCombLGauss00", CosThetaL, bkgGauss_mean0, bkgGauss_sigma0);
	RooGaussian f_bkgCombLGauss11("f_bkgCombLGauss11","f_bkgCombLGauss11", CosThetaL, bkgGauss_mean1, bkgGauss_sigma1);
	RooGaussian f_bkgCombLGauss22("f_bkgCombLGauss22","f_bkgCombLGauss22", CosThetaL, bkgGauss_mean2, bkgGauss_sigma2);
	RooGaussian f_bkgCombLGauss44("f_bkgCombLGauss44","f_bkgCombLGauss44", CosThetaL, bkgGauss_mean4, bkgGauss_sigma4);
	RooGaussian f_bkgCombLGauss66("f_bkgCombLGauss66","f_bkgCombLGauss66", CosThetaL, bkgGauss_mean6, bkgGauss_sigma6);
	RooGaussian f_bkgCombLGauss77("f_bkgCombLGauss77","f_bkgCombLGauss77", CosThetaL, bkgGauss_mean7, bkgGauss_sigma7);
	RooGaussian f_bkgCombLGauss88("f_bkgCombLGauss88","f_bkgCombLGauss88", CosThetaL, bkgGauss_mean8, bkgGauss_sigma8);
	RooGaussian f_bkgCombLGauss99("f_bkgCombLGauss99","f_bkgCombLGauss99", CosThetaL, bkgGauss_mean9, bkgGauss_sigma9);
	RooGaussian f_bkgCombLGauss1010("f_bkgCombLGauss1010","f_bkgCombLGauss1010", CosThetaL, bkgGauss_mean10, bkgGauss_sigma10);
	
	RooRealVar bkg_frac("bkg_frac","bkg_frac",1.,0.,1.);

	RooAddPdf *f_bkgCombL = 0;
	switch (iBin) {
		case 9:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss99, f_bkgCombL_P), bkg_frac);
			break;
		case 4:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss44, f_bkgCombL_P), bkg_frac);
			break;
		case 0:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss00, f_bkgCombL_P), bkg_frac);
			break;
		case 10:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss1010, f_bkgCombL_P), bkg_frac);
			break;
		case 1:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss11, f_bkgCombL_P), bkg_frac);
			break;
		case 2:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss22, f_bkgCombL_P), bkg_frac);
			break;
		case 6:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss66, f_bkgCombL_P), bkg_frac);
			break;
		case 7:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss77, f_bkgCombL_P), bkg_frac);
			break;
		case 8:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss88, f_bkgCombL_P), bkg_frac);
			break;
	}

//  Get data and apply unbinned fit
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.384 || Bmass < 5.174)",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.349 || Bmass < 5.209)",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.379 || Bmass < 5.179)",Q2range[iBin]),0);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && ( (Bmass > 5.349 && Bmass < 5.458) || (Bmass < 5.209))",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && ( (Bmass > 5.379 && Bmass < 5.458) || (Bmass < 5.179))",Q2range[iBin]),0);
//	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.458)",Q2range[iBin]),0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	if ( iBin == 1 || iBin == 0) {
	    for (int i = 0; i < 1; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 8 || iBin == 7 || iBin == 2) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 6 || iBin == 9) {
	    for (int i = 0; i < 3; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 4) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
//	if (iBin == 9) {
//	    for (int i = 0; i < 3; i++) {
//		    CosThetaL = 1 + 1 * .08;
//		    d->add(RooArgSet(CosThetaL));
//		    CosThetaL = -1 - 1 * .08;
//		    d->add(RooArgSet(CosThetaL));
//	    } 
//    }
	if (iBin == 4) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .08;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .08;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 10) {
	    for (int i = 0; i < 20; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Extended(),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Extended(),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
//	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*data,Save(kTRUE),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) return;
    int NNCOMB[11] = {273, 591, 841, 0, 503, 0, 164, 164, 164, 1227, 269};
    TString otoyspath = TString::Format("./RootFiles");
	TFile *fout = new TFile(TString::Format("%s/%s_cosl_bin%d.root",otoyspath.Data(),outfile,iBin),"RECREATE");
////    TH1* prior = d1->createHistogram("CosThetaL",22);
////    TH1* prior_gen = f_bkgCombL->createHistogram("CosThetaL",22);
    TH1* prior = d1->createHistogram("data",CosThetaL,Binning(22,-1.1,1.1));
    TH1* prior_gen = f_bkgCombL->createHistogram("pdf",CosThetaL,Binning(22,-1.1,1.1));
    prior_gen->Scale(NNCOMB[iBin]);
	fout->Write();
	fout->Close();
	
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
//	data->plotOn(framecosl,Binning(25)); 
	data->plotOn(framecosl,Binning(22)); 
	f_bkgCombL->plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
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
	delete t1;
	delete c;
	delete data;
// Prepare datacard
    keepParam = 0;
	if (keepParam){
		double val[3] = {0,0,0};
		val[0] = bkgCombL_c0.getVal();val[1] = bkgCombL_c0.getError();
		writeParam(iBin, "bkgCombL_c0", val);
		val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
		writeParam(iBin, "bkgCombL_c1", val);
		val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
		writeParam(iBin, "bkgCombL_c2", val);
		val[0] = 0; val[1] = 0;
		writeParam(iBin, "bkgCombL_c3", val);
		if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
			val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
			writeParam(iBin, "bkgCombL_c3", val);
	    }
		val[0] = 0; val[1] = 0;
		if ( iBin == 1) {
			val[0] = bkgGauss_mean1.getVal();val[1] = bkgGauss_mean1.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 2 ) {
			val[0] = bkgGauss_mean2.getVal();val[1] = bkgGauss_mean2.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma2.getVal();val[1] = bkgGauss_sigma2.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 10) {
			val[0] = bkgGauss_mean10.getVal();val[1] = bkgGauss_mean10.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma10.getVal();val[1] = bkgGauss_sigma10.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		} else if (iBin == 7 ) {
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		} else if (iBin == 9 ) {
			val[0] = bkgGauss_mean9.getVal();val[1] = bkgGauss_mean9.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma9.getVal();val[1] = bkgGauss_sigma9.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 6 ) {
			val[0] = bkgGauss_mean6.getVal();val[1] = bkgGauss_mean6.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma6.getVal();val[1] = bkgGauss_sigma6.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 4 ) {
			val[0] = bkgGauss_mean4.getVal();val[1] = bkgGauss_mean4.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma4.getVal();val[1] = bkgGauss_sigma4.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 8 ) {
			val[0] = bkgGauss_mean8.getVal();val[1] = bkgGauss_mean8.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma8.getVal();val[1] = bkgGauss_sigma8.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	    } else if (iBin == 0) {
			val[0] = bkgGauss_mean0.getVal();val[1] = bkgGauss_mean0.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma0.getVal();val[1] = bkgGauss_sigma0.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		}
		val[0] = bkg_frac.getVal();val[1] = bkg_frac.getError();
		writeParam(iBin, "bkg_frac", val);
	}	
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_prior_C(int iBin, const char outfile[] = "angular2D_prior", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
//	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1.1, 1.1);
	
	// Create combinatorial background distribution
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",0.,-3.,3.);
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",0.,-3.,3.);
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",0.,-3.,3.);
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
    if (iBin == 2 || iBin == 7 || iBin == 8 || iBin == 9) {
        f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
    }
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
//	RooChebychev f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	// Create peak background distribution
	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}",  0.2,  0.08,  0.25);
	RooRealVar bkgGauss_mean2("bkgGauss_mean2","cos#theta_{l}",  0.2, -0.2,  0.5);
	RooRealVar bkgGauss_mean3("bkgGauss_mean3","cos#theta_{l}",  0.3,  0.1,  0.35);
	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.6, -0.9, -0.5);
	RooRealVar bkgGauss_mean5("bkgGauss_mean5","cos#theta_{l}", -0.2, -0.3, -0.1);
	RooRealVar bkgGauss_mean6("bkgGauss_mean6","cos#theta_{l}",  0.01, -0.2,  0.2);
	RooRealVar bkgGauss_mean7("bkgGauss_mean7","cos#theta_{l}",  0.2, 0.01,  0.3);
	RooRealVar bkgGauss_mean8("bkgGauss_mean8","cos#theta_{l}", -0.7, -0.9, -0.5);
	RooRealVar bkgGauss_mean9("bkgGauss_mean9","cos#theta_{l}",  0.1,  0.0,  0.2);
	
	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",  .10,  .06,  0.2);
	RooRealVar bkgGauss_sigma2("bkgGauss_sigma2","#sigma_{2}",  .20,  .01,  2);
	RooRealVar bkgGauss_sigma3("bkgGauss_sigma3","#sigma_{3}",  .10,  .06, 0.3);
	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  0.6,  .1,   1);
	RooRealVar bkgGauss_sigma5("bkgGauss_sigma5","#sigma_{5}",  .20,  .01,0.4);
	RooRealVar bkgGauss_sigma6("bkgGauss_sigma6","#sigma_{6}",  .12,  .10,0.3);
	RooRealVar bkgGauss_sigma7("bkgGauss_sigma7","#sigma_{7}",  .20,  .01,0.4);
	RooRealVar bkgGauss_sigma8("bkgGauss_sigma8","#sigma_{8}",  .20,  .1,  0.5);
	RooRealVar bkgGauss_sigma9("bkgGauss_sigma9","#sigma_{9}",  .20,  .1,   1);
	
	RooGaussian f_bkgCombLGauss11("f_bkgCombLGauss11","f_bkgCombLGauss11", CosThetaL, bkgGauss_mean1, bkgGauss_sigma1);
	RooGaussian f_bkgCombLGauss22("f_bkgCombLGauss22","f_bkgCombLGauss22", CosThetaL, bkgGauss_mean2, bkgGauss_sigma2);
	RooGaussian f_bkgCombLGauss33("f_bkgCombLGauss33","f_bkgCombLGauss33", CosThetaL, bkgGauss_mean3, bkgGauss_sigma3);
	RooGaussian f_bkgCombLGauss44("f_bkgCombLGauss44","f_bkgCombLGauss44", CosThetaL, bkgGauss_mean4, bkgGauss_sigma4);
	RooGaussian f_bkgCombLGauss55("f_bkgCombLGauss55","f_bkgCombLGauss55", CosThetaL, bkgGauss_mean5, bkgGauss_sigma5);
	RooGaussian f_bkgCombLGauss66("f_bkgCombLGauss66","f_bkgCombLGauss66", CosThetaL, bkgGauss_mean6, bkgGauss_sigma6);
	RooGaussian f_bkgCombLGauss77("f_bkgCombLGauss77","f_bkgCombLGauss77", CosThetaL, bkgGauss_mean7, bkgGauss_sigma7);
	RooGaussian f_bkgCombLGauss88("f_bkgCombLGauss88","f_bkgCombLGauss88", CosThetaL, bkgGauss_mean8, bkgGauss_sigma8);
	RooGaussian f_bkgCombLGauss99("f_bkgCombLGauss99","f_bkgCombLGauss99", CosThetaL, bkgGauss_mean9, bkgGauss_sigma9);
	
	RooGaussian f_bkgCombLGauss62("f_bkgCombLGauss62","f_bkgCombLGauss62", CosThetaL, bkgGauss_mean6, bkgGauss_sigma2);
	RooGaussian f_bkgCombLGauss42("f_bkgCombLGauss42","f_bkgCombLGauss42", CosThetaL, bkgGauss_mean4, bkgGauss_sigma2);
	RooGaussian f_bkgCombLGauss31("f_bkgCombLGauss31","f_bkgCombLGauss31", CosThetaL, bkgGauss_mean3, bkgGauss_sigma1);
	RooRealVar bkg_frac_G("bkg_frac_G","bkg_frac_G",1.,0.,1.);
	RooAddPdf f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", RooArgList(f_bkgCombLGauss77,f_bkgCombLGauss66), bkg_frac_G);
	
	RooRealVar bkg_frac("bkg_frac","bkg_frac",1.,0.,1.);

	RooAddPdf *f_bkgCombL = 0;
	switch (iBin) {
		case 9:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss11, f_bkgCombL_P), bkg_frac);
			break;
		case 4:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss55, f_bkgCombL_P), bkg_frac);
			break;
		case 0:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss22, f_bkgCombL_P), bkg_frac);
			break;
		case 10:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss99, f_bkgCombL_P), bkg_frac);
			break;
		case 1:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss77, f_bkgCombL_P), bkg_frac);
			break;
		case 2:
			f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss66, f_bkgCombL_P), bkg_frac);
			break;
		case 6:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss44, f_bkgCombL_P), bkg_frac);
			break;
		case 7:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss88, f_bkgCombL_P), bkg_frac);
			break;
		case 8:
		    f_bkgCombL = new RooAddPdf("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss33, f_bkgCombL_P), bkg_frac);
			break;
	}

//  Get data and apply unbinned fit
////	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && ( (Bmass > 5.349 && Bmass < 5.458) || (Bmass < 5.209))",Q2range[iBin]),0);
////	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.209)",Q2range[iBin]),0);
////	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.349 && Bmass < 5.458) ",Q2range[iBin]),0);
////	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.349 && Bmass > 5.209) ",Q2range[iBin]),0);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.458) ",Q2range[iBin]),0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	if ( iBin == 1 || iBin == 0) {
	    for (int i = 0; i < 1; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 8 || iBin == 7 || iBin == 2) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .03;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .03;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 6 || iBin == 9) {
	    for (int i = 0; i < 3; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	if (iBin == 4) {
	    for (int i = 0; i < 4; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
//	if (iBin == 9) {
//	    for (int i = 0; i < 3; i++) {
//		    CosThetaL = 1 + 1 * .08;
//		    d->add(RooArgSet(CosThetaL));
//		    CosThetaL = -1 - 1 * .08;
//		    d->add(RooArgSet(CosThetaL));
//	    } 
//    }
////	if (iBin == 4) {
////	    for (int i = 0; i < 4; i++) {
////		    CosThetaL = 1 + 1 * .08;
////		    d->add(RooArgSet(CosThetaL));
////		    CosThetaL = -1 - 1 * .08;
////		    d->add(RooArgSet(CosThetaL));
////	    } 
////    }
	if (iBin == 10) {
	    for (int i = 0; i < 20; i++) {
		    CosThetaL = 1 + 1 * .05;
		    d->add(RooArgSet(CosThetaL));
		    CosThetaL = -1 - 1 * .05;
		    d->add(RooArgSet(CosThetaL));
	    } 
    }
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
	f_fitresult->Print();
////	if (f_fitresult->status() != 0) return;
	
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
//	data->plotOn(framecosl,Binning(25)); 
	data->plotOn(framecosl,Binning(22)); 
	f_bkgCombL->plotOn(framecosl); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
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
	delete t1;
	delete c;
	delete data;
// Prepare datacard
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void angular2D_prior_Check(int iBin, const char outfile[] = "angular2D_prior", bool keepParam = true)
{//{{{
	setTDRStyle();
	// Fit to signal simulation by YsSm+YcCm to determine Sm
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	

//  Get data and apply unbinned fit
////	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && ( (Bmass > 5.349 && Bmass < 5.458) || (Bmass < 5.209))",Q2range[iBin]),0);
	RooDataSet *data_L = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.209)",Q2range[iBin]),0);
	RooDataSet *data_H = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.349) ",Q2range[iBin]),0);
////	RooDataSet *data_H = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.349 && Bmass < 5.458) ",Q2range[iBin]),0);
	RooDataSet *data_S = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.349 && Bmass > 5.209) ",Q2range[iBin]),0);
////	RooDataSet *data_T = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.458) ",Q2range[iBin]),0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
////	RooDataSet *data_L = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.219)",Q2range[iBin]),0);
////	RooDataSet *data_H = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass > 5.339 && Bmass < 5.458) ",Q2range[iBin]),0);
////	RooDataSet *data_S = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && (Bmass < 5.339 && Bmass > 5.219) ",Q2range[iBin]),0);
	
	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	RooPlot* framecosl = CosThetaL.frame(); 
////	data_S->plotOn(framecosl,Binning(20),RooFit::Name("Signal region"),MarkerColor(2),LineColor(2),LineWidth(2)); 
	data_S->plotOn(framecosl,Binning(10),RooFit::Name("Signal"),MarkerStyle(20),MarkerColor(2),LineColor(2),LineWidth(2)); 
	data_L->plotOn(framecosl,Binning(10),RooFit::Name("Lower"),MarkerStyle(22),MarkerColor(3),LineColor(3),LineWidth(2)); 
	data_H->plotOn(framecosl,Binning(10),RooFit::Name("Higher"),MarkerStyle(23),MarkerColor(4),LineColor(4),LineWidth(2)); 
////	data_T->plotOn(framecosl,Binning(10),RooFit::Name("Rest"),MarkerStyle(33),MarkerColor(6),LineColor(6),LineWidth(2)); 
	framecosl->SetTitle("");
	framecosl->SetMinimum(0);
	framecosl->SetAxisRange(-1.,1.,"X"); 
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->Draw();

	TLegend *leg =new TLegend(0.71,0.69,0.89,0.86,NULL);
	leg->AddEntry("Signal"," [5.209, 5.349]"," PEL ");
	leg->AddEntry("Lower"," [5.100, 5.209]"," PEL ");
	leg->AddEntry("Higher"," [5.349, 5.600]"," PEL ");
////	leg->AddEntry("Higher"," [5.349, 5.458]"," PEL ");
////	leg->AddEntry("Rest"," [5.458, 5.600]"," PEL ");
////	leg->AddEntry("Signal"," [5.219, 5.339]"," PEL ");
////	leg->AddEntry("Lower"," [5.100, 5.219]"," PEL ");
////	leg->AddEntry("Higher"," [5.339, 5.458]"," PEL ");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->Draw();
	
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
	delete t1;
	delete c;
	delete data_S;
	delete data_L;
	delete data_H;
// Prepare datacard
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
