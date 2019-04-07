// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------

//void angular2D_bin(int iBin, const char outfile[] = "angular2D")
std::vector<double> angular2D_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular2D")
//std::vector<double> angular_reco_bin(int iBin, const char outfile[] = "angular_reco")
{//{{{
    setTDRStyle();
	// Remark: You must use RooFit!! It's better in unbinned fit.
	//         Extended ML fit is adopted by Mauro, just follow!!
	//         Need some modification for accXrecoEff.
	cout<<endl<<"iBin = "<<iBin<<endl<<endl; 
	// Create parameters and PDFs
////	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
////	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar CosThetaL("CosThetaL", "", -1., 1.);
	RooRealVar Bmass("Bmass","",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
////////////////////////////////////// scan   //////////////////////////////////.........................
////	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
////	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
////////////////////////////////////// refit  //////////////////////////////////.........................
		
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
//	RooRealVar effP7("effP7","effP7",readParam(iBin,"accXrecoEff", 7)); 
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
//	effP10.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
//	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
//	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
//	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
	
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
	// // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.26,5.30);
//	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.269,5.289);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
	sigGauss1_sigma.setConstant(kTRUE);
//	sigGauss1_sigma.setError(readParam(iBin,"sigGauss1_sigma",1));
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
	sigGauss2_sigma.setConstant(kTRUE);
//	sigGauss2_sigma.setError(readParam(iBin,"sigGauss2_sigma",1));
	RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
	sigM_frac.setConstant(kTRUE);
//	sigM_frac.setError(readParam(iBin,"sigM_frac",1));
	// // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
	
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format;
//	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
//	} else {
////		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
//	}	
	if (iBin != 0 && iBin != 1 && iBin != 9) {
//	if (iBin != 0 ) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	//	f_sigA_argset.add(RooArgSet(effP7, effP8, effP9));
	//	f_rec_format = "( effP0*exp(-0.5* (((CosThetaL-effP1)/effP2)**2)) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 ) ";
	//	f_sigA_argset.add(RooArgSet(effP10));
	//	f_rec_format = "( effP0 + effP1*exp(-0.5* (((CosThetaL-effP2)/effP3)**2)) ) * ( effP4+effP5*CosThetaL+effP6*CosThetaL**2+effP7*CosThetaL**3+effP8*CosThetaL**4+effP9*CosThetaL**5+effP10*CosThetaL**6 ) ";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	//	f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
	//	f_sigA_argset.add(RooArgSet(effP7, effP8, effP9));
	//	f_rec_format = "( effP0*exp(-0.5* (((CosThetaL-effP1)/effP2)**2)) ) * ( effP3+effP4*CosThetaL+effP5*CosThetaL**2+effP6*CosThetaL**3+effP7*CosThetaL**4+effP8*CosThetaL**5+effP9*CosThetaL**6 ) ";
	//	f_sigA_argset.add(RooArgSet(effP10));
	//	f_rec_format = "( effP0 + effP1*exp(-0.5* (((CosThetaL-effP2)/effP3)**2)) ) * ( effP4+effP5*CosThetaL+effP6*CosThetaL**2+effP7*CosThetaL**3+effP8*CosThetaL**4+effP9*CosThetaL**5+effP10*CosThetaL**6 ) ";
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
    f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
    bkgCombL_c0.setConstant(kTRUE);
    bkgCombL_c1.setConstant(kTRUE);
    bkgCombL_c2.setConstant(kTRUE);
    if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
        f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
        bkgCombL_c3.setConstant(kTRUE);
    }
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	// Create peak background distribution
	RooRealVar bkgGauss_mean("bkgGauss_mean","cos#theta_{l}", readParam(iBin,"bkgGauss_mean",0));
	bkgGauss_mean.setConstant(kTRUE);
    RooRealVar bkgGauss_sigma("bkgGauss_sigma","#sigma",  readParam(iBin,"bkgGauss_sigma",0));
	bkgGauss_sigma.setConstant(kTRUE);
//	RooRealVar bkgGauss_mean_1("bkgGauss_mean_1","cos#theta_{l}", readParam(iBin,"bkgGauss_mean_1",0));
//	bkgGauss_mean_1.setConstant(kTRUE);
//  RooRealVar bkgGauss_sigma_1("bkgGauss_sigma_1","#sigma",  readParam(iBin,"bkgGauss_sigma_1",0));
//	bkgGauss_sigma_1.setConstant(kTRUE);
    RooRealVar bkg_frac("bkg_frac","bkg_frac",readParam(iBin,"bkg_frac",0));
	bkg_frac.setConstant(kTRUE);
	
	RooGaussian f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", CosThetaL, bkgGauss_mean, bkgGauss_sigma);
//	RooGaussian f_bkgCombLGauss_1("f_bkgCombLGauss_1","f_bkgCombLGauss_1", CosThetaL, bkgGauss_mean_1, bkgGauss_sigma_1);
//	RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombLGauss_1), bkg_frac);
	RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombL_P), bkg_frac);
	
	RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL, f_bkgCombM);
	cout<<">>>>>>>>>>>>>>>> INFO: f_bkgComb prepared. <<<<<<<<<<<<<<<<<<<<<<"<<endl;
//////  bin2,4,10,  9
    const double insig[11] = {120,  300,  750,  0, 350,  0, 200, 250, 200,  750, 2200};
    const double inbkg[11] = {400, 1000, 1800,  0, 700,  0, 200, 250, 200, 2300, 4800};
	RooRealVar nsig("nsig","nsig",insig[iBin],0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",inbkg[iBin],0,8E3);
//////	bin 0,1,6,7,8
////	RooRealVar nsig("nsig","nsig",50,0,4E3);
////	RooRealVar nbkgComb("nbkgComb","nbkgComb",100,0,6E3);
	
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
	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3));	

	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	


//
//  RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2","Simplex"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit2"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),ExternalConstraints(gausConstraints),Minimizer("Minuit"), Strategy(2), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3) {
        cout<<"STATUS = "<<f_fitresult->status()<<";  "<<f_fitresult->covQual()<<endl;
//	    if (f_fitresult->status() != 0) {
	    std::vector<double> output;
	    output.push_back(fh.getVal());
	    output.push_back(fh.getError());
	    output.push_back(afb.getVal());
	    output.push_back(afb.getError());
	    return output;
	}

//  // Fitting procedure in TMinuit
//  double isMigradConverge[2] = {-1,0};
//  double isMinosValid = -1;
//  RooAbsReal *nll = f.createNLL(*data,Extended(kTRUE),Offset(kFALSE),NumCPU(1));// Minos and Save are unknown.
//  RooMinuit minuit(*nll);
//  printf("INFO\t\t: Start MIGRAD loop\n");
//  for(int iLoop = 0; iLoop < 10; iLoop++){
//      isMigradConverge[0] = minuit.migrad();
//      printf("INFO\t\t: MIGRAD return code=%.0f\n",isMigradConverge[0]);
//      if (isMigradConverge[0] == 0) break;
//  }
//  isMigradConverge[1] = minuit.save()->minNll();
//  int gKeepParam = 0;
//  if (gKeepParam) {
//      writeParam(iBin, "migrad", isMigradConverge);
//      double val[4]={0,0,0,0};
////    val[0] = toBoundedFh( fh.getVal() );
//      val[1] = fh.getError();val[2]=fh.getErrorLo();val[3]=fh.getErrorHi();
//      writeParam(iBin, "fh_migrad", val, 4);
////    val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//      val[1] = afb.getError();val[2]=afb.getErrorLo();val[3]=afb.getErrorHi();
//      writeParam(iBin, "afb_migrad",val, 4);
//  }
//  double isHesseValid = minuit.hesse();
//  // Keep HESSE result as preliminary
//  if (gKeepParam) {
//      writeParam(iBin, "hesse", &isHesseValid, 1);
//      minuit.save();
//      double val[4]={0,0,0,0};
////    val[0] = toBoundedFh( fh.getVal() );
//      val[1] = fh.getError();val[2]=fh.getErrorLo();val[3]=fh.getErrorHi();
//      writeParam(iBin, "fh_hesse", val, 4);
////    val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//      val[1] = afb.getError();val[2]=afb.getErrorLo();val[3]=afb.getErrorHi();
//      writeParam(iBin, "afb_hesse",val, 4);
//  }
//  printf("INFO\t\t: Start MINOS loop\n");
//  for(int iLoop = 0; iLoop < 3; iLoop++){
//      isMinosValid = minuit.minos(RooArgSet(afb,fh,nsig));
//      printf("INFO\t\t: MINOS return code=%.0f\n",isMinosValid);
//      if (isMinosValid == 0) break;
//  }
//  if (gKeepParam) {
//      writeParam(iBin, "minos", &isMinosValid, 1);
//  }
//  minuit.save();


	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
////	TPad *tp1 = new TPad("tp1","",0.0,0.10,1.0,1.00);
	TPad *tp1 = new TPad("tp1","",0.0,0.20,1.0,1.00);
	TPad *tp2 = new TPad("tp2","",0.0,0.00,1.0,0.29);
	tp2->SetTopMargin(0.00);
	tp2->SetBottomMargin(0.30);
    tp1->Draw();
	tp2->Draw();
	tp1->cd();
	RooPlot* framemass = Bmass.frame();
	data->plotOn(framemass,RooFit::Name("data"), XErrorSize(0), Binning(20)); 
	f.plotOn(framemass,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullmass = framemass->pullHist();
	//f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
//	f.plotOn(framemass, Components(f_sig),FillStyle(3001),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2), DrawOption("F"), FillStyle(3005),FillColor(4),VLines());
////	f.plotOn(framemass, RooFit::Name("sig"),Components(f_sig), DrawOption("F"), FillStyle(3005),FillColor(4),VLines());
	f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2));
	f.plotOn(framemass,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(2),LineWidth(4),LineStyle(4));
	
	framemass->SetTitle("");
	framemass->GetYaxis()->SetTitle("Events / 0.025 GeV");
	framemass->SetMinimum(0);
	framemass->SetMaximum(framemass->GetMaximum() * 1.2);
	framemass->GetXaxis()->SetLabelOffset(0.85);
	framemass->SetTitleOffset(0.85, "Y");
////	framemass->GetXaxis()->SetOption("U");
	framemass->Draw();
    
////	TLegend *leg =new TLegend(0.17,0.69,0.35,0.86,NULL,"brNDC");
	TLegend *leg =new TLegend(0.67,0.63,0.88,0.88,NULL,"brNDC");
	leg->AddEntry("data"," Data "," PL ");
	leg->AddEntry("pdf"," Total fit"," L ");
	leg->AddEntry("sig"," Signal "," LF ");
	leg->AddEntry("bkgComb"," Background"," L");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
	leg->SetTextFont(22);
	leg->Draw();

	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t1->DrawLatex(.20,.90,TString::Format("Preliminary"));
	t1->DrawLatex(.18,.82,TString::Format("q^{2}"));
	t2->DrawLatex(.22,.82,Q2String[iBin]);
    if (iBin == 10) {
        t2->DrawLatex(.22,.77, "  14.18 - 22.00 GeV^{2}");
////        t2->DrawLatex(.18,.75, TString::Format("Signal yield: %.0f #pm %.0f", nsig.getVal(), nsig.getError() ));
    } else {
////        t2->DrawLatex(.18,.79, TString::Format("Signal yield: %.0f #pm %.0f", nsig.getVal(), nsig.getError() ));
    }
	t2->DrawLatex(.74,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	
	tp2->cd();
	RooPlot *framepullmass = Bmass.frame();
	framepullmass->addPlotable(pullmass,"P");
	framepullmass->GetYaxis()->SetRangeUser(-4.8,4.8);
	framepullmass->GetYaxis()->SetNdivisions(5);
////	framepullmass->GetYaxis()->SetNdivisions(0, 3.3, 0, kTRUE);
	framepullmass->GetXaxis()->SetTitle("#it{m}(K^{+}#mu^{+}#mu^{-}) [GeV]");
	framepullmass->GetYaxis()->SetTitle("Pull");
	framepullmass->SetTitleSize(0.16, "X");
	framepullmass->SetTitleSize(0.14, "Y");
	framepullmass->SetTitleFont(22, "XY");
	framepullmass->SetTitleOffset(0.29, "Y");
	framepullmass->SetTitleOffset(0.85, "X");
	framepullmass->SetTickLength(0.08,"X");
	framepullmass->SetTickLength(0.02,"Y");
	framepullmass->SetLabelFont(22,"XY");
	framepullmass->SetLabelSize(0.13,"X");
	framepullmass->SetLabelSize(0.15,"Y");
	TLine *tl1 =new TLine(5.10, 0., 5.60, 0.);
	tl1->SetLineColor(1);
	tl1->SetLineWidth(1.5);
	framepullmass->Draw();
	tl1->Draw();

	cout<<endl<<"Chi^{2}_mass  = "<<framemass->chiSquare()<<endl<<endl; 
	c->Update();
	c->Print(TString::Format("./plots/%s_bin%d_Index_%d.pdf",outfile,iBin,Index));
	
	// Draw projection to CosThetaL
	tp1->cd();
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,RooFit::Name("data"), XErrorSize(0), Binning(20)); 
	f.plotOn(framecosl,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullcosl = framecosl->pullHist();
	//f.plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
	//f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),FillStyle(3005),FillColor(4),VLines(), DrawOption("F"));
//	f.plotOn(framecosl, Components(f_sig),FillStyle(3001),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2), DrawOption("F"), FillStyle(3005),FillColor(4),VLines());
////	f.plotOn(framecosl, RooFit::Name("sig"),Components(f_sig), DrawOption("F"), FillStyle(3005),FillColor(4),VLines());
	f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(4),LineWidth(2));
	f.plotOn(framecosl,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(2),LineWidth(4),LineStyle(4));

	
	framecosl->SetTitle("");
	framecosl->GetYaxis()->SetTitle("Events / 0.1");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.2);
	framecosl->GetXaxis()->SetLabelOffset(0.85);
	framecosl->SetTitleOffset(0.85, "Y");
	framecosl->Draw();
	leg->Draw();
	t1->DrawLatex(.20,.90,TString::Format("Preliminary"));
	t1->DrawLatex(.18,.82,TString::Format("q^{2}"));
	t2->DrawLatex(.22,.82,Q2String[iBin]);
    if (iBin == 10) {
        t2->DrawLatex(.22,.77, "  14.18 - 22.00 GeV^{2}");
    }
	t2->DrawLatex(.74,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));

	tp2->cd();
	RooPlot *framepullcosl = CosThetaL.frame();
	framepullcosl->addPlotable(pullcosl,"P");
	framepullcosl->GetYaxis()->SetRangeUser(-4.8,4.8);
	framepullcosl->GetYaxis()->SetNdivisions(5);
	framepullcosl->GetXaxis()->SetTitle("cos#it{#theta_{l}}");
	framepullcosl->GetYaxis()->SetTitle("Pull");
	framepullcosl->SetTitleSize(0.16, "X");
	framepullcosl->SetTitleSize(0.14, "Y");
	framepullcosl->SetTitleFont(22, "XY");
	framepullcosl->SetTitleOffset(0.29, "Y");
	framepullcosl->SetTitleOffset(0.85, "X");
	framepullcosl->SetLabelFont(22,"XY");
	framepullcosl->SetTickLength(0.08,"X");
	framepullcosl->SetTickLength(0.02,"Y");
	framepullcosl->SetLabelSize(0.13,"X");
	framepullcosl->SetLabelSize(0.15,"Y");
	TLine *tl2 =new TLine(-1., 0., 1., 0.);
	tl2->SetLineColor(1);
	tl2->SetLineWidth(1.5);
	framepullcosl->Draw();
	tl2->Draw();      
	
	cout<<endl<<"Chi^{2}_cosl  = "<<framecosl->chiSquare()<<endl<<endl; 
	c->Update();
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.pdf",outfile,iBin,Index));

	delete c;
	delete t1;
	delete data;


// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
        // map of signal
        bkgCombM_c.setConstant(kTRUE);
        offset.setConstant(kTRUE);
        fh.setConstant(kTRUE);
        afb.setConstant(kTRUE);
        nsig.setConstant(kTRUE);
        nbkgComb.setConstant(kTRUE);
        sigGauss_mean.setConstant(kTRUE);
        RooWorkspace *wspace2 = new RooWorkspace("wspace","wspace");
        wspace2->import(f);
        wspace2->writeToFile(TString::Format("%s/wspace_pdf_bin%d.root",owspacepath.Data(),iBin),true);
	//	val[0] = fh.getVal();val[1] = fh.getError();val[2] = fh.getError();
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, "fh", val, 3);
	//	val[0] = afb.getVal();val[1] = afb.getError();val[2] = afb.getError();
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, "afb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "FCN", val);
        val[0] = nsig.getVal(); val[1] = nsig.getError();
		writeParam(iBin, "Nsig", val);
        val[0] = nbkgComb.getVal(); val[1] = nbkgComb.getError();
		writeParam(iBin, "NbkgComb", val);
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
		writeOutput(outfile,iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(outfile,iBin, Index, "fh", val);
		val[1]=0; val[2]=0;
//////	bin 0,1,6,7,8,9
////		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrained
//////	bin 2,4,10
		val[0] = afb.getVal();   // Unconstrained
		writeOutput(outfile,iBin, Index, "F_afb", val);
////		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrained
		val[0] = fh.getVal();    // Unconstrained
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

void angular2D( const char outfile[] = "angular2D")
{//{{{
	setTDRStyle();
//	doFit = false;
//	double SMAfb[11]={0.00, 0.00,  0.00, 0.00,  0.00, 0.00,  0.00,  0.00, 0.00, 0.00,  0.00};   // .....
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
////	double conYA[2]   = {-0.05,  -0.05};
////	double conYAerr[2]= { 0.45,   0.45};
////	double conYF[2]   = { 0.70,   0.70};
////	double conYFerr[2]= { 0.80,   0.80};
	double conYA[2]   = { 0.00,   0.00};
	double conYAerr[2]= { 1.00,   1.00};
	double conYF[2]   = { 1.40,   1.40};
	double conYFerr[2]= { 1.60,   1.60};

//	double SMFh[9]    = { 0.10,  0.07,    0.04,  0.00,   0.02,   0.02,   0.03, 0.05, 0.02}; // LHCb
//	double SMFherr[7] = { 0.008, 0.001, 0.0005,  0.00, 0.0005, 0.0005, 0.0015};
	double CBFh[7]     = { 0.04710, 0.02370, 0.007050,  0.006920, 0.008050, 0.007400, 0.025400}; // Christoph Bobeth
	double CBFhuerr[7] = { 0.00028, 0.00029, 0.000370,  0.000410, 0.000610, 0.000470, 0.002900};
	double CBFhderr[7] = { 0.00023, 0.00024, 0.000420,  0.000450, 0.000610, 0.000510, 0.002800};
	double CB[7]   ={1.50, 3.15, 15.09, 17.0, 20.0, 18.5, 3.5};
	double CBerr[7]={0.50, 1.15,  0.91,  1.0,  2.0,  3.5, 2.5};
	double IFAEFh[7]    = { 0.046, 0.023, 0.012,  0.007, 0.007, 0.009, 0.025}; // Matias from IFAE
	double IFAEFherr[7] = { 0.005, 0.003, 0.001,  0.001, 0.001, 0.001, 0.003};
	double IFAE[7]   ={1.50, 3.15, 6.49, 15.09, 17.0, 20.0, 3.5};
	double IFAEerr[7]={0.50, 1.15, 2.19,  0.91,  1.0,  2.0, 2.5};
	double TUMFh[8]    = { 0.0445, 0.0216, 0.0110,  0.0068, 0.0060, 0.0060, 0.0075, 0.0196}; // flavio from TUM
	double TUMFherr[8] = { 0.00022, 0.00020, 0.00020,  0.00017, 0.00015, 0.00017, 0.00024, 0.00025};
	double TUM[8]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5};
	double TUMerr[8]={0.50, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5};

	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.50, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
//	double yfh[9], yerrfh[9], yafb[9], yerrafb[9];
	double yfh[9], yuerrfh[9],yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9];
    // approval old
////	double sysfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
////	double sysafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
    // approval new
	double sysafb[9] = { 0.0523, 0.0747, 0.0186, 0.0538, 0.0155, 0.0252, 0.0249, 0.0331, 0.0296 };
	double sysfh[9]  = { 0.3889, 0.1379, 0.0442, 0.0616, 0.0715, 0.0660, 0.0857, 0.0923, 0.0569 };

	double dataUfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
	double dataDfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
	double dataUafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
	double dataDafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
	for (int i = 0; i < 9; i++) {
	//	yfh[i]  = 0; yerrfh[i]  = 0; 
		yfh[i]  = 0; yuerrfh[i] = 0; yderrfh[i] = 0;
		yafb[i] = 0; yuerrafb[i]= 0; yderrafb[i]= 0;
	}
// Checkout input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
	//	if (ibin != 0 && ibin != 1 && ibin != 6 && ibin != 9) continue;
	//	if (ibin == 1) continue;
		yafb[i]      = readParam(ibin,"afb",0);
		yuerrafb[i]  = fabs(readParam(ibin,"FCErrAfb",1));
		yderrafb[i]  = fabs(readParam(ibin,"FCErrAfb",0));
		yfh[i]       = readParam(ibin,"fh",0);
	//	yerrfh[i]   = fabs(readParam(ibin,"fh",1));
		yuerrfh[i]   = fabs(readParam(ibin,"FCErrFh",1));
		yderrfh[i]   = fabs(readParam(ibin,"FCErrFh",0));
        dataUfh[i] = sysfh[i] + yuerrfh[i];
        dataDfh[i] = sysfh[i] + yderrfh[i];
        if (dataDfh[i] > yfh[i]) dataDfh[i] = yfh[i];
        dataUafb[i] = sysafb[i] + yuerrafb[i];
        dataDafb[i] = sysafb[i] + yderrafb[i];

		if (yderrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
	//	else { yderrfh[i] = yuerrfh[i]; }
		printf("Afb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
	//	printf("Fh [%d]=%6.4f +- %6.4f\n",ibin,yfh[i], yerrfh[i]);
		printf("Fh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i], yuerrfh[i], yderrfh[i]);
	}
//	Draw
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitleOffset(1.1,"XY");
	frame->SetTitle("");
	frame->Draw("2");	
	
	frame->SetXTitle("#it{q}^{2} [GeV^{2}]");
	frame->SetYTitle("#it{F_{H}}");
	//frame->SetAxisRange(-0.005,0.3,"Y");
////	frame->SetAxisRange(-0.10,1.40,"Y");
	frame->SetAxisRange(-0.20,3.00,"Y");
//	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yerrfh,yerrfh);
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
////	d_fh->SetMarkerColor(1);
////	d_fh->SetMarkerStyle(20);
////	d_fh->Draw("P");
////	d_fh->SetFillColor(4);
////	d_fh->SetFillStyle(3144);
////	d_fh->Draw("2");
    d_fh->SetLineWidth(3);
	d_fh->Draw("p[]");
////	TGraphAsymmErrors *sys_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,sysfh,sysfh);
////	sys_fh->SetFillColor(4);
////	sys_fh->SetFillStyle(3244);
////	sys_fh->Draw("2");
	TGraphAsymmErrors *sys_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,dataDfh,dataUfh);
	sys_fh->SetMarkerSize(1.8);
    sys_fh->SetLineWidth(3);
    sys_fh->SetFillColor(1);
    sys_fh->SetFillStyle(1001);
    sys_fh->SetMarkerColor(1);
    sys_fh->SetMarkerStyle(20);
	sys_fh->Draw("Pz");
	TGraphAsymmErrors *IFAE_fh  = new TGraphAsymmErrors(6,IFAE,IFAEFh,IFAEerr,IFAEerr,IFAEFherr,IFAEFherr);
	IFAE_fh->SetFillColor(2);
	IFAE_fh->Draw("2");
	TGraphAsymmErrors *TUM_fh  = new TGraphAsymmErrors(7,TUM,TUMFh,TUMerr,TUMerr,TUMFherr,TUMFherr);
	TUM_fh->SetFillColor(3);
////	TUM_fh->Draw("2");
//	TGraphAsymmErrors *s_fh  = new TGraphAsymmErrors(7,x,SMFh,xerr,xerr,SMFherr,SMFherr);
	TGraphAsymmErrors *s_fh  = new TGraphAsymmErrors(5,CB,CBFh,CBerr,CBerr,CBFhderr,CBFhuerr);
	s_fh->SetFillColor(6);
////	s_fh->Draw("2");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(17);
	c_fh->SetFillStyle(3001);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();

////	TLegend *leg =new TLegend(0.65,0.56,0.90,0.86);
////	TLegend *leg =new TLegend(0.65,0.74,0.90,0.86);
	TLegend *leg =new TLegend(0.68,0.70,0.90,0.86);
	leg->AddEntry(sys_fh,"    Data","lep");
////	leg->AddEntry(d_fh,"stat. uncertainty","f");
////	leg->AddEntry(sys_fh,"total uncertainty","e");
////	leg->AddEntry(sys_fh," CMS data","lep");
////	leg->AddEntry(d_fh," F-C uncertainty","f");
////	leg->AddEntry(sys_fh," Systematic uncertainty ","f");
////	leg->AddEntry(s_fh," #LT EOS #GT","f");
	leg->AddEntry(IFAE_fh," #LT DHMV #GT","f");
////	leg->AddEntry(TUM_fh," #LT FLAVIO #GT","f");
////	leg->AddEntry(s_fh," #LT SM-EOS #GT","f");
////	leg->AddEntry(IFAE_fh," #LT SM-DHMV #GT","f");
////	leg->AddEntry(TUM_fh," #LT SM-FLAVIO #GT","f");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->Draw();
	
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t1->DrawLatex(.23,.90,TString::Format("Preliminary"));
	t2->DrawLatex(.68,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	//c->Print(TString::Format("./plots/%s_fh_SM.pdf",outfile));
	c->Clear();
	
	frame->SetXTitle("#it{q}^{2} [GeV^{2}]");
	frame->SetYTitle("#it{A_{FB}}");
//	frame->SetAxisRange(-0.2,0.2,"Y");
////	frame->SetAxisRange(-0.4,0.4,"Y");
	frame->SetAxisRange(-1.0,1.0,"Y");
	frame->Draw();
	TGraphAsymmErrors *d_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
////	d_afb->SetMarkerColor(1);
////	d_afb->SetMarkerStyle(20);
////	d_afb->Draw("P");
//	d_afb->SetLineWidth(0.005);	
////	d_afb->SetFillColor(4);
////	d_afb->SetFillStyle(3144);
////	d_afb->Draw("2");
    d_afb->SetLineWidth(3);
	d_afb->Draw("p[]");
////	TGraphAsymmErrors *sys_afb  = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,sysafb,sysafb);
////	sys_afb->SetFillColor(4);
////	sys_afb->SetFillStyle(3244);
////	sys_afb->Draw("2");
	TGraphAsymmErrors *sys_afb  = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,dataDafb,dataUafb);
	sys_afb->SetMarkerSize(1.8);
    sys_afb->SetLineWidth(3);
    sys_afb->SetFillColor(1);
    sys_afb->SetFillStyle(1001);
    sys_afb->SetMarkerColor(1);
    sys_afb->SetMarkerStyle(20);
	sys_afb->Draw("Pz");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(17);
	c_afb->SetFillStyle(3001);
	c_afb->Draw("2");
   
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(1);
	tl1->Draw();
////	TLegend *leg_1 =new TLegend(0.65,0.74,0.90,0.86);
	TLegend *leg_1 =new TLegend(0.68,0.78,0.90,0.86);
	leg_1->AddEntry(sys_afb,"    Data","lep");
////	leg_1->AddEntry(d_afb,"stat. uncertainty","f");
////	leg_1->AddEntry(sys_afb,"total uncertainty","e");
////	leg_1->AddEntry(sys_afb," CMS data","lep");
////	leg_1->AddEntry(d_afb," F-C uncertainty","f");
////	leg_1->AddEntry(d_afb," Data w/ F&C uncertainty","lep");
////	leg_1->AddEntry(sys_afb," Systematic uncertainty ","f");
	leg_1->SetLineColor(0);
	leg_1->SetFillColor(0);
	leg_1->SetBorderSize(0);
	leg_1->SetTextSize(0.04);
	leg_1->Draw();
	
	t1->DrawLatex(.23,.90,TString::Format("Preliminary"));
	t2->DrawLatex(.68,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}


void Average( const char outfile[] = "Average")
{//{{{
	setTDRStyle();
	float genRfh[9], genRafb[9], genRfherr[9], genRafberr[9], genRafb_w[9], genRfh_w[9];
	float recofh[9], recoafb[9], recofherr[9], recoafberr[9], recoafb_w[9], recofh_w[9];
	float Sum_genRafb_w = 0, Sum_genRfh_w = 0, Sum_genRafb_wx = 0, Sum_genRfh_wx = 0;
	float Sum_recoafb_w = 0, Sum_recofh_w = 0, Sum_recoafb_wx = 0, Sum_recofh_wx = 0;
	float A_genRfh, A_genRfherr, A_genRafb, A_genRafberr;
	float A_recofh, A_recofherr, A_recoafb, A_recoafberr;
	
	for(int i = 0, ibin = 0; i < 7 && ibin < 9; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
	//	cout<<"iBin = "<<ibin<<endl;
		recoafb[i]      = readParam(ibin,"recoafb",0);
		recoafberr[i]   = fabs(readParam(ibin,"recoafb",1));
		recofh[i]       = readParam(ibin,"recofh",0);
		recofherr[i]    = fabs(readParam(ibin,"recofh",1));
		
		genRafb[i]      = readParam(ibin,"genafb_R",0);
		genRafberr[i]   = fabs(readParam(ibin,"genafb_R",1));
		genRfh[i]       = readParam(ibin,"genfh_R",0);
		genRfherr[i]    = fabs(readParam(ibin,"genfh_R",1));
		
		recoafb_w[i] = 1. / ( recoafberr[i] * recoafberr[i] );
		recofh_w[i]  = 1. / ( recofherr[i]  * recofherr[i] );
		Sum_recoafb_w  = Sum_recoafb_w + recoafb_w[i];
		Sum_recofh_w   = Sum_recofh_w  + recofh_w[i];
		Sum_recoafb_wx = Sum_recoafb_wx + recoafb_w[i] * recoafb[i];
		Sum_recofh_wx  = Sum_recofh_wx  + recofh_w[i] * recofh[i];
		genRafb_w[i] = 1. / ( genRafberr[i] * genRafberr[i] );
		genRfh_w[i]  = 1. / ( genRfherr[i]  * genRfherr[i] );
		Sum_genRafb_w = Sum_genRafb_w + genRafb_w[i];
		Sum_genRfh_w  = Sum_genRfh_w  + genRfh_w[i];
		Sum_genRafb_wx = Sum_genRafb_wx + genRafb_w[i] * genRafb[i];
		Sum_genRfh_wx  = Sum_genRfh_wx  + genRfh_w[i] * genRfh[i];
	}
	A_genRfh     = Sum_genRfh_wx  / Sum_genRfh_w; 
	A_genRafb    = Sum_genRafb_wx / Sum_genRafb_w; 
	A_genRfherr  = 1. / sqrt(Sum_genRfh_w);
	A_genRafberr = 1. / sqrt(Sum_genRafb_w);
	printf("A_genAfb_R =%6.6f +- %6.6f\n",A_genRafb,A_genRafberr);
	printf("A_genFh_R  =%6.6f +- %6.6f\n",A_genRfh,A_genRfherr);
	
	A_recofh     = Sum_recofh_wx  / Sum_recofh_w; 
	A_recoafb    = Sum_recoafb_wx / Sum_recoafb_w; 
	A_recofherr  = 1. / sqrt(Sum_recofh_w);
	A_recoafberr = 1. / sqrt(Sum_recoafb_w);
	printf("A_recoAfb  =%6.6f +- %6.6f\n",A_recoafb,A_recoafberr);
	printf("A_recoFh   =%6.6f +- %6.6f\n",A_recofh,A_recofherr);

// write output
	double val[2]={0,0};
	val[0] = A_genRafb; val[1] = A_genRafberr;
	writeParam(999, "A_genRafb", val);
	val[0] = A_genRfh; val[1] = A_genRfherr;
	writeParam(999, "A_genRfh", val);
	val[0] = A_recoafb; val[1] = A_recoafberr;
	writeParam(999, "A_recoafb", val);
	val[0] = A_recofh; val[1] = A_recofherr;
	writeParam(999, "A_recofh", val);
	 
}
void BranchFraction( int iBin, const char outfile[] = "BranchFraction")
{//{{{
	setTDRStyle();
	double BRsig[10], BRsigErr[10];
	double q2Width[9]   = {  1.00,   2.30,   4.38,    2.77,    1.82,   2.00,   4.00,   5.00,    21.0};
	double Nsig[9]      = { 169.3,  330.9,  785.1,   364.7,   214.7,  261.9,  225.9,  778.3,  2286.1};
	double NsigErr[9]   = {  21.9,   31.7,   41.8,    28.6,    19.3,   21.0,   20.4,   46.5,    72.5};
	double nsigMC[9]    = { 52492, 124757, 235372,  126589,   77442,  79723,  75543, 268733,  771916};
	double nsigMCErr[9] = {   229,    353,    485,     356,     278,    282,    275,    518,     879};
	double NsigMC[9]    = { 832529, 1938823, 3941409, 2703747, 1747022, 1748609, 2033121, 4236991, 18441527};  // 1. full
//	double NsigMC[9]    = { 832333, 1938388, 3940533, 2703291, 1746758, 1748383, 2032863, 4236039, 14942549};  // 2. eta, pt
//	double NsigMC[9]    = { 818990, 1915353, 3909301, 2688241, 1738252, 1740501, 2023748, 4186452, 14834386};  // 3. eta, pt, Beta < 2.5

	double Njpsi    = 788213;	 double NjpsiErr   = 841;  // 1. loose
//	double Njpsi    = 406422;	 double NjpsiErr   = 638;  // 1. tight
	double NjpsiMC  = 9707803;  double NjpsiMCErr = 0;    // 1. full
//	double NjpsiMC  = 9712053;  double NjpsiMCErr = 0;    // 1. full no
//	double NjpsiMC  = 9448460;  double NjpsiMCErr = 0;    // 2. eta, pt
	double njpsiMC  = 1076117;  double njpsiMCErr = 982; // 1. loose
//	double njpsiMC  = 590551;    double njpsiMCErr = 199; // 1. tight
	
    double Npsip    = 68830;	 double NpsipErr   = 179; // 1. loose 
//   double Npsip    = 33411;	 double NpsipErr   = 183; // 1. tight
	double NpsipMC  = 9850824;  double NpsipMCErr = 0;    // 1. full
//	double NpsipMC  = 9449549;  double NpsipMCErr = 0;    // 2. eta, pt
	double npsipMC  = 1052067;  double npsipMCErr = 1026; // 1. loose
//	double npsipMC  = 543415;   double npsipMCErr = 737;  // 1. tight

	double BRjpsiK  = 1.026e-3;	double BRjpsiKErr = 0.031e-3;
//	double BRjpsiK  = 1.016e-3;	double BRjpsiKErr = 0.034e-3;
	double BRJpsi   = 5.961e-2;   double BRJpsiErr  = 0.033e-2;
//	double BRJpsi   = 5.93e-2;    double BRJpsiErr  = 0.06e-2;
	double BRpsip   = 7.900e-3;   double BRpsipErr  = 0.9e-3;

	double lumiData = 20.5;

// Acceptance
	double gensig[9]   = {  9637,  22367,  46309,  32063,  20132,  20667,  23870,  49087, 175045 };  // 1. eta, pt, Beta  
	double Gensig[9]   = {290114, 643113,1185047, 684319, 376578, 344678, 349424,1392595,3873273 };  // 1. eta
//	double gensig[9]   = {  9805,  22645,  46697,  32246,  20218,  20775,  23979,  49703, 176365 };  // 4. eta, pt  
//	double Gensig[9]   = {479958,1093224,2066139,1212858, 671508, 614788, 624217,2368482,6762692 };  // 4. full
	double genjpsi  = 33125;     double Genjpsi    = 763031; // 1. eta, pt, Beta   
	double genpsip  = 15188;     double Genpsip    = 301155; // 1. eta, pt, Beta
//	double genjpsi  = 15979;     double Genjpsi    = 648914; // 4. signalMC full   
//	double genpsip  = 15248;     double Genpsip    = 535431; // 4. signalMC full
	double BRpsipjpsi    = 0; double BRkmumujpsi    = 0;
	double BRpsipjpsiErr = 0; double BRkmumujpsiErr = 0;
	double CSjpsi   = 0;      double CSjpsiErr    = 0;
	
	double d_Nsig2[9],  d_Esig2[9];
	double d_gjpsi2, d_gpsip2;
    double d_gsig2[9];
	double d_BjpsiK2, d_BJpsi2, d_Ejpsi2, d_Njpsi2, d_Epsip2, d_Npsip2, d_Bpsip2;//, d_BpsipK2;
	for(int i = 0; i<9; i++) d_gsig2[i]   = (Gensig[i] - gensig[i]) / ( Gensig[i] * Gensig[i] * gensig[i]);	
    d_gjpsi2  = (Genjpsi - genjpsi) / ( Genjpsi * Genjpsi * genjpsi);
	d_gpsip2  = (Genpsip - genpsip) / ( Genpsip * Genpsip * genpsip);
	d_Ejpsi2  = (NjpsiMC - njpsiMC) / ( NjpsiMC * NjpsiMC * njpsiMC );  
	d_Epsip2  = (NpsipMC - npsipMC) / ( NpsipMC * NpsipMC * npsipMC );  
	d_Njpsi2  = NjpsiErr   * NjpsiErr   / ( Njpsi   * Njpsi );
	d_Npsip2  = NpsipErr   * NpsipErr   / ( Npsip   * Npsip );
	d_BjpsiK2 = BRjpsiKErr * BRjpsiKErr / ( BRjpsiK * BRjpsiK );
//	d_BpsipK2 = BRpsipKErr * BRpsipKErr / ( BRpsipK * BRpsipK );
	d_BJpsi2  = BRJpsiErr  * BRJpsiErr  / ( BRJpsi  * BRJpsi );
	d_Bpsip2  = BRpsipErr  * BRpsipErr  / ( BRpsip  * BRpsip );
	

	BRsig[9] = 0.;
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
//	for(int i = 0; i < 9; i++){
		cout<<"iBin = "<<ibin<<endl;
		if (i != 8) { 
			BRsig[i] = 1. * pow(10,7)  * Nsig[i] * ( genjpsi/Genjpsi ) * (njpsiMC / NjpsiMC) * BRjpsiK * BRJpsi / ( q2Width[i] * Njpsi * ( gensig[i] / Gensig[i]) * nsigMC[i] / NsigMC[i] );
		//	BRsig[i] = 1. * pow(10,7)  * Nsig[i] * (njpsiMC / NjpsiMC) * BRjpsiK * BRJpsi / ( q2Width[i] * Njpsi * nsigMC[i] / NsigMC[i] );
		//	BRsig[i] = 1. * Nsig[i] * (njpsiMC / NjpsiMC) * BRjpsiK * BRJpsi / (  Njpsi * nsigMC[i] / NsigMC[i] );
		} else { 
			BRsig[i] = 1. * pow(10,7)  * Nsig[i] * ( genjpsi/Genjpsi ) * (njpsiMC / NjpsiMC) * BRjpsiK * BRJpsi / ( Njpsi * ( gensig[8] / Gensig[8]) * nsigMC[i] / NsigMC[i] );
		//	BRsig[i] = 1. * pow(10,7)  * Nsig[i] * (njpsiMC / NjpsiMC) * BRjpsiK * BRJpsi / ( Njpsi * nsigMC[i] / NsigMC[i] );
		}
		d_Nsig2[i]  = NsigErr[i] * NsigErr[i] / (Nsig[i] * Nsig[i]);
		d_Esig2[i]  = (NsigMC[i] - nsigMC[i]) /( NsigMC[i] * NsigMC[i] * nsigMC[i] );
		
		BRpsipjpsi  = Npsip * ( genjpsi/Genjpsi ) * ( njpsiMC/NjpsiMC ) * BRJpsi / (( genpsip / Genpsip) * ( npsipMC/NpsipMC ) * Njpsi * BRpsip ); 
		BRpsipjpsiErr = BRpsipjpsi * sqrt( d_Npsip2 + d_gpsip2 + d_Epsip2 + d_gjpsi2 + d_Ejpsi2 + d_Njpsi2 + d_BjpsiK2 );
	//	BRpsipjpsi  = Npsip * ( njpsiMC/NjpsiMC ) * BRJpsi / ( ( npsipMC/NpsipMC ) * Njpsi * BRpsip ); 
	//	BRpsipjpsiErr = BRpsipjpsi * sqrt( d_Npsip2 + d_Epsip2 + d_Ejpsi2 + d_Njpsi2 + d_BjpsiK2 );
	//	BRpsipjpsiErr = BRpsipjpsi * sqrt( d_Npsip2 + d_gpsip2 + d_Epsip2 + d_gjpsi2 + d_Ejpsi2 + d_Njpsi2 + d_BjpsiK2 + d_Bpsip2);
			
		BRkmumujpsi = Nsig[8] * ( genjpsi/Genjpsi ) * ( njpsiMC/NjpsiMC ) * BRJpsi / (( gensig[8] / Gensig[8]) * ( nsigMC[8]/NsigMC[8] ) * Njpsi * 1. );
		BRkmumujpsiErr = BRkmumujpsi * sqrt( d_Nsig2[8] + d_gsig2[8] + d_Esig2[8] + d_gjpsi2 + d_Ejpsi2 + d_Njpsi2 + d_BjpsiK2 );
	//	BRkmumujpsi = Nsig[8] * ( njpsiMC/NjpsiMC ) * BRJpsi / ( ( nsigMC[8]/NsigMC[8] ) * Njpsi * 1. );
	//	BRkmumujpsiErr = BRkmumujpsi * sqrt( d_Nsig2[8]  + d_Esig2[8]  + d_Ejpsi2 + d_Njpsi2 + d_BjpsiK2 );

		CSjpsi = 1e-9 * Njpsi / ( 2. * ( genjpsi/Genjpsi ) * ( njpsiMC/NjpsiMC ) * lumiData * BRjpsiK * BRJpsi );
		CSjpsiErr = CSjpsi * sqrt( d_Njpsi2 + d_gjpsi2 + d_Ejpsi2 + d_BjpsiK2 + d_BJpsi2 );
	//	CSjpsi = 1e-9 * Njpsi / ( ( njpsiMC/NjpsiMC ) * lumiData * BRjpsiK * BRJpsi );
	//	CSjpsiErr = CSjpsi * sqrt( d_Njpsi2 + d_Ejpsi2 + d_BjpsiK2 + d_BJpsi2 );

		BRsigErr[i] = BRsig[i] * sqrt( d_Nsig2[i] + d_gjpsi2 + d_Ejpsi2 + d_Njpsi2 + d_gsig2[i] + d_Esig2[i] + d_BjpsiK2 + d_BJpsi2 );
	//	BRsigErr[i] = BRsig[i] * sqrt( d_Nsig2[i] + d_Ejpsi2 + d_Njpsi2 + d_Esig2[i] + d_BjpsiK2 + d_BJpsi2 );
		cout<<"iBR = "<<BRsig[i]<<" +- "<<BRsigErr[i]<<endl;
	//	cout<<sqrt( d_Nsig2[i] + d_Ejpsi2 - d_Njpsi2 - d_Esig2[i] + d_BjpsiK2 + d_BJpsi2 )<<endl;
	//	cout<<d_Nsig2[i] + d_Ejpsi2 - d_Njpsi2 - d_Esig2[i] + d_BjpsiK2 + d_BJpsi2<<endl;
	//	cout<<d_Nsig2[i]<<endl;
	//	cout<<d_Ejpsi2<<endl;
	//	cout<<d_Njpsi2<<endl;
	//	cout<<d_Esig2[i]<<endl;
	//	cout<<d_BjpsiK2<<endl;
	//	cout<<d_BJpsi2<<endl;
	}
	cout<<endl<<"psipK / jpsiK = "<<BRpsipjpsi<<" +- "<<BRpsipjpsiErr<<endl;
	cout<<"mumuK / jpsiK = ("<<BRkmumujpsi*1000<<" +- "<<BRkmumujpsiErr*1000<<") * E-3"<<endl;
	cout<<endl<<"JpsiK CR  = ("<<CSjpsi<<" +- "<<CSjpsiErr<<") * ub"<<endl<<endl;
	
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYF[2]   = { 0.35,   0.35};
	double conYFerr[2]= { 0.35,   0.35};
	double SMBR[6]    = { 0.381, 0.381, 0.371, 0.228, 0.193, 0.095};
	double SMBRerr[6] = { 0.086, 0.092, 0.097, 0.075, 0.073, 0.032};
	double LHCbBR[7]    = { 0.285, 0.249, 0.229, 0.204, 0.207, 0.177, 0.078};
	double LHCbBRerr[7] = { 0.027, 0.023, 0.016, 0.018, 0.020, 0.018, 0.010};
	double LHCb3BR[17]    = { 0.332, 0.233, 0.282, 0.254, 0.221, 0.231, 0.245, 0.231, 0.177, 0.193, 0.161, 0.164, 0.206, 0.137, 0.074, 0.059, 0.043 }; 
	double LHCb3BRerr[17] = { 0.018, 0.015, 0.016, 0.015, 0.014, 0.014, 0.014, 0.014, 0.013, 0.012, 0.010, 0.010, 0.011, 0.010, 0.008, 0.007, 0.007 }; 
	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={ 0.5, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double SMx[6]   ={ 1.025, 3.15, 6.49, 15.09, 17.0, 20.0};
	double SMxerr[6]={ 0.975, 1.15, 2.19,  0.91,  1.0,  2.0};
	double LHCbx[9]   ={1.025, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double LHCbxerr[9]={0.975, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double LHCb3x[17]   ={ 0.54, 1.55, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50, 11.40, 12.15, 15.50, 16.50, 17.50, 18.50, 19.50, 20.50, 21.50 };
	double LHCb3xerr[17]={ 0.44, 0.45, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50,  0.40,  0.35,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50 };
	
    double CBBR[7]     = { 0.349, 0.816, 0.399,  0.358, 0.392, 0.984, 1.760}; // Christoph Bobeth
	double CBBRuerr[7] = { 0.071, 0.130, 0.040,  0.045, 0.033, 0.078, 0.310};
	double CBBRderr[7] = { 0.052, 0.130, 0.036,  0.025, 0.040, 0.100, 0.260};
	double CB[7]   ={1.50, 3.15, 15.09, 17.0, 20.0, 18.5, 3.5};
	double CBerr[7]={0.50, 1.15,  0.91,  1.0,  2.0,  3.5, 2.5};
	double IFAEBR[7]    = { 0.353, 0.801, 1.490,  0.417, 0.372, 0.384, 1.740}; // Matias from IFAE
	double IFAEBRerr[7] = { 0.103, 0.241, 0.507,  0.062, 0.052, 0.052, 0.528};
	double IFAE[7]   ={1.50, 3.15, 6.49, 15.09, 17.0, 20.0, 3.5};
	double IFAEerr[7]={0.50, 1.15, 2.19,  0.91,  1.0,  2.0, 2.5};
	double TUMBR[8]     = { 0.369, 0.366, 0.357, 0.339, 0.258, 0.209, 0.108, 0.}; // flavio from TUM
	double TUMBRerr[8]  = { 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.};
	double TUM[8]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5};
	double TUMerr[8]={0.50, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5};
//	Draw
	TCanvas *c = new TCanvas();
	
	TGraphAsymmErrors *s_br  = new TGraphAsymmErrors(6,SMx,SMBR,SMxerr,SMxerr,SMBRerr,SMBRerr);
	s_br->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
	s_br->GetYaxis()->SetTitle("dB/dq^{2} [10^{-7}/GeV^{2}]");
	s_br->GetXaxis()->SetLimits(0.,22);
	s_br->SetMinimum(0.);
    s_br->SetMaximum(0.7);
	s_br->SetFillColor(46);
	s_br->Draw("A2");
	TGraphAsymmErrors *IFAE_br  = new TGraphAsymmErrors(6,IFAE,IFAEBR,IFAEerr,IFAEerr,IFAEBRerr,IFAEBRerr);
	IFAE_br->SetFillColor(3);
////	IFAE_br->Draw("2");
	TGraphAsymmErrors *cb_br  = new TGraphAsymmErrors(5,CB,CBBR,CBerr,CBerr,CBBRderr,CBBRuerr);
	cb_br->SetFillColor(6);
////	cb_br->Draw("2");
	TGraphAsymmErrors *l_br  = new TGraphAsymmErrors(7,LHCbx,LHCbBR,LHCbxerr,LHCbxerr,LHCbBRerr,LHCbBRerr);
	l_br->SetFillColor(4);
	l_br->SetFillStyle(3002);
	l_br->Draw("2");
	TGraphAsymmErrors *l_br3  = new TGraphAsymmErrors(17,LHCb3x,LHCb3BR,LHCb3xerr,LHCb3xerr,LHCb3BRerr,LHCb3BRerr);
	l_br3->SetFillColor(3);
	l_br3->SetFillStyle(3008);
	l_br3->Draw("2");
	TGraphAsymmErrors *TUM_br  = new TGraphAsymmErrors(7,TUM,TUMBR,TUMerr,TUMerr,TUMBRerr,TUMBRerr);
	TUM_br->SetFillColor(2);
////	TUM_br->Draw("2");
	TGraphAsymmErrors *d_br  = new TGraphAsymmErrors(7,x,BRsig,xerr,xerr,BRsigErr,BRsigErr);
	d_br->SetMarkerColor(1);
	d_br->SetMarkerStyle(20);
	d_br->Draw("P");
	TGraphAsymmErrors *c_br  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_br->SetFillColor(1);
	c_br->SetFillStyle(3003);
	c_br->Draw("2");

//	c->BuildLegend();
	TLegend *leg =new TLegend(0.66,0.69,0.90,0.86);
	leg->AddEntry(d_br," CMS 8TeV ","lep");
	leg->AddEntry(l_br," LHCb 7TeV ","f");
	leg->AddEntry(l_br3," LHCb 7TeV+8TeV ","f");
	leg->AddEntry(s_br," SM prediction ","f");
////	leg->AddEntry(cb_br," Bobeth et al. Theorists ","f");
////	leg->AddEntry(IFAE_br," Matias et al. Theorists ","f");
////	leg->AddEntry(TUM_br," Straub et al. Theorists ","f");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->Draw();
	
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t2->DrawLatex(.68,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.15,.90,TString::Format("CMS"));
	c->Print(TString::Format("./plots/%s_br.pdf",outfile));
	c->Clear();

}

double tosigPDF(double afb, double fh, double cosL){
    return (fh/2. + 3./4.*(1.-fh)*cosL) - 1./4.*(1-fh)*cosL*cosL*cosL + 0.5*afb*cosL*cosL ;
}
void read_pdf()
{
    ofstream myfile;
    myfile.open (TString::Format("./RootFiles/Para_bins.txt"));
    for (int iBin = 0; iBin <= 10; iBin++) {
        if (iBin == 3 || iBin == 5) continue;
        double fh = readParam(iBin,"fh", 0);
        double afb = readParam(iBin,"afb", 0);
        TH1 *sig__CosThetaL = 0;
    ////    TH1 *comb__CosThetaL = 0;
    ////    TH1 *pdf__CosThetaL = 0;
        TH1F *h2_eff_fine = 0;
    	TFile *fout = new TFile(TString::Format("./RootFiles/angular2D_bin%d.root",iBin));
    	TFile *feff = new TFile(TString::Format("./RootFiles/Efficiency_%d.root",iBin));
        h2_eff_fine = (TH1F*)feff->Get("h2_eff_fine");
        sig__CosThetaL = (TH1*)fout->Get("sig__CosThetaL");
    ////    pdf__CosThetaL = (TH1*)fout->Get("pdf__CosThetaL");
    ////    myfile.open (TString::Format("./RootFiles/Para_bin%d.txt",iBin));
        myfile << TString::Format("bin%d\n",iBin);
        myfile << TString::Format("AFB=%f\tFH=%f\n",afb,fh);
        myfile << "AngularFormat:\t[F_H/2 + 3/4*(1-F_H)*(1-CosThetaL^2) + AFB*CosThetaL]\n";
////        myfile << "SignalPDF:\t[ (F_H/2 + 3/4*(1-F_H)*CosThetaL) -(1./4.*(1-F_H)*CosThetaL^3) + 1./2.*AFB*CosThetaL^2 ]\n";
////        myfile << "Index\tCosThetaL\tSignalYield\tSignalPDF\tEfficiency\tSignalPDF*Efficiency\n";
        myfile << "Index\tCosThetaL\tSignalYield\tEfficiency\n";
        for (int i=1; i<=20; i++) {
            double sig_i = sig__CosThetaL->GetBinContent(i);
            double eff_i = h2_eff_fine->GetBinContent(i);
            double cosL = i*0.1-1.05;
////            double sig_pdf = tosigPDF(afb, fh, cosL+0.05) - tosigPDF(afb, fh, cosL-0.05); 
////            double pdfXeff = sig_pdf*eff_i;
////            myfile << TString::Format("%d\t%f\t%f\t%f\t%f\t%f\n",i,cosL,sig_i,sig_pdf,eff_i,pdfXeff);
            myfile << TString::Format("%d\t%f\t%f\t%f\n",i,cosL,sig_i,eff_i);
        }
        myfile << "\n\n";
    }
    myfile.close();
}


void angular2D_LHCb( const char outfile[] = "angular2D")
{//{{{
	setTDRStyle();
//	doFit = false;
//	double SMAfb[11]={0.00, 0.00,  0.00, 0.00,  0.00, 0.00,  0.00,  0.00, 0.00, 0.00,  0.00};   // .....
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = {-0.05,  -0.05};
	double conYAerr[2]= { 0.45,   0.45};
	double conYF[2]   = { 0.70,   0.70};
	double conYFerr[2]= { 0.80,   0.80};
////	double conYA[2]   = { 0.00,   0.00};
////	double conYAerr[2]= { 1.00,   1.00};
////	double conYF[2]   = { 1.40,   1.40};
////	double conYFerr[2]= { 1.60,   1.60};


	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.50, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9],yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9];
    // approval old
////	double sysfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
////	double sysafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
    // approval new
	double sysafb[9] = { 0.0523, 0.0747, 0.0186, 0.0538, 0.0155, 0.0252, 0.0249, 0.0331, 0.0296 };
	double sysfh[9]  = { 0.3889, 0.1379, 0.0442, 0.0616, 0.0715, 0.0660, 0.0857, 0.0923, 0.0569 };

	double dataUfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
	double dataDfh[9]  = { 0.0752, 0.0871, 0.0435, 0.0292, 0.0313, 0.0199, 0.0391, 0.0731, 0.0276 };
	double dataUafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
	double dataDafb[9] = { 0.0379, 0.0198, 0.0185, 0.0118, 0.0093, 0.0096, 0.0176, 0.0179, 0.0117 };
	for (int i = 0; i < 9; i++) {
	//	yfh[i]  = 0; yerrfh[i]  = 0; 
		yfh[i]  = 0; yuerrfh[i] = 0; yderrfh[i] = 0;
		yafb[i] = 0; yuerrafb[i]= 0; yderrafb[i]= 0;
	}
// Checkout input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]      = readParam(ibin,"afb",0);
		yuerrafb[i]  = fabs(readParam(ibin,"FCErrAfb",1));
		yderrafb[i]  = fabs(readParam(ibin,"FCErrAfb",0));
		yfh[i]       = readParam(ibin,"fh",0);
		yuerrfh[i]   = fabs(readParam(ibin,"FCErrFh",1));
		yderrfh[i]   = fabs(readParam(ibin,"FCErrFh",0));
        dataUfh[i] = sysfh[i] + yuerrfh[i];
        dataDfh[i] = sysfh[i] + yderrfh[i];
        if (dataDfh[i] > yfh[i]) dataDfh[i] = yfh[i];
        dataUafb[i] = sysafb[i] + yuerrafb[i];
        dataDafb[i] = sysafb[i] + yderrafb[i];

		if (yderrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		printf("Afb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("Fh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i], yuerrfh[i], yderrfh[i]);
	}

    ////////////////////////// LHCb 1/fb////////////////////////////////////
    const Int_t n2 = 7;
    //Afb
    Float_t x2[n2]  = {1.025, 3.15, 6.49, 11.475, 15.09, 17.0, 20.0};
    Float_t ex2l[n2] = {.975, 1.15, 2.19, 1.385, 0.91, 1.0, 2.0};
    Float_t ex2h[n2] = {.975, 1.15, 2.19, 1.385, 0.91, 1.0, 2.0};
    Float_t y2[n2]  =  {0.00, 0.07, -0.02, -0.03, -0.01, -0.09, 0.02};
    //stat uncertainty
    Float_t ey2h[n2] = {0.06, 0.08, 0.03, 0.07,  0.12,  0.07,  0.11};
    Float_t ey2l[n2] = {0.05, 0.05, 0.05 ,0.07,  0.06,  0.09,  0.11};
    // total uncertainty
    Float_t ey02h[n2] = {0.06708, 0.08246, 0.04243, 0.07071,  0.12042,  0.07280,  0.11045};
    Float_t ey02l[n2] = {0.05099, 0.05099, 0.05831 ,0.07071,  0.06083,  0.09055,  0.11045};
    //Fh
    Float_t x2_F[n2]  = {1.025, 3.15, 6.49, 11.475, 15.09, 17.0, 20.0};
    Float_t ex2l_F[n2] = {.975, 1.15, 2.19, 1.385, 0.91, 1.0, 2.0};
    Float_t ex2h_F[n2] = {.975, 1.15, 2.19, 1.385, 0.91, 1.0, 2.0};
    Float_t y2_F[n2]  =  {0.00, 0.14, 0.04, 0.11, 0.08, 0.18, 0.14};
    //stat uncertainty
    Float_t ey2l_F[n2] = {0.00, 0.10, 0.04 ,0.08,  0.08,  0.14,  0.14};
    Float_t ey2h_F[n2] = {0.12, 0.16, 0.10, 0.20,  0.28,  0.22,  0.31};
    // total uncertainty
    Float_t ey02l_F[n2] = {0.0000, 0.1020, 0.0566, 0.0806, 0.0806, 0.1456, 0.1414};
    Float_t ey02h_F[n2] = {0.1236, 0.1649, 0.1166, 0.2010, 0.2807, 0.2202, 0.3102};

    //////////////////////////LHCb 3/fb //////////////////////////////////
    const Int_t n3 = 17;
    //Afb
    Float_t x3[n3] = {0.54, 1.55, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50, 11.375, 12.125, 15.50, 16.50, 17.50, 18.50, 19.50, 20.50, 21.50};
    Float_t ex3l[n3] = {0.44, 0.45, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.375, 0.375, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
    Float_t ex3h[n3] = {0.44, 0.45, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.375, 0.375, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
    Float_t y3[n3] = {-0.05, 0.05, 0.06, -0.005, 0.02, 0.00, 0.015, 0.00, 0.075, 0.025, -0.065, -0.025, -0.03, 0.01, 0.015, 0.03, 0.11 };
    // Afb range
    Float_t ey3l_0[n3] = {-0.09, 0.00, 0.01, -0.02, -0.01, -0.04, -0.01, -0.02, 0.03, 0.00, -0.10, -0.05, -0.06, -0.03, -0.02, -0.01, 0.03};
    Float_t ey3h_0[n3] = {-0.01, 0.10, 0.11,  0.01,  0.05,  0.04,  0.04,  0.02, 0.12, 0.05, -0.03,  0.00,  0.00,  0.05,  0.05,  0.07, 0.19};
    //stat uncertainty
    Float_t ey3l[n3] = {-0.09, 0.00, 0.01, -0.02, -0.01, -0.04, -0.01, -0.02, 0.03, 0.00, -0.10, -0.05, -0.06, -0.03, -0.02, -0.01, 0.03};
    Float_t ey3h[n3] = {-0.01, 0.10, 0.11,  0.01,  0.05,  0.04,  0.04,  0.02, 0.12, 0.05, -0.03,  0.00,  0.00,  0.05,  0.05,  0.07, 0.19};
    // syst uncertainty
    Float_t ey03_s[n3] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02};
    // total uncertainty
    Float_t ey03l[n3] = {0.10, 0.01, 0.01, 0.015, 0.03, 0.04, 0.025, 0.02, 0.045, 0.025, 0.035, 0.025, 0.03, 0.04, 0.035, 0.04, 0.08};
    Float_t ey03h[n3] = {0.02, 0.11, 0.05, 0.015, 0.03, 0.04, 0.025, 0.02, 0.045, 0.025, 0.035, 0.025, 0.03, 0.04, 0.035, 0.04, 0.08};
    //Fh
    Float_t x3_F[n3] = {0.54, 1.55, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50, 11.375, 12.125, 15.50, 16.50, 17.50, 18.50, 19.50, 20.50, 21.50};
    Float_t ex3l_F[n3] = {0.44, 0.45, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.375, 0.375, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
    Float_t ex3h_F[n3] = {0.44, 0.45, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.375, 0.375, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
    Float_t y3_F[n3] = {0.105, 0.105, 0.175, 0.02, 0.045, 0.07, 0.04, 0.015, 0.145, 0.05, 0.13, 0.06, 0.085, 0.14, 0.05, 0.07, 0.225};
    // Fh range
    Float_t ey3l_0_F[n3] = {0.01, 0.00, 0.05, 0.00, 0.00, 0.00,0.00, 0.00, 0.06, 0.00, 0.06, 0.00, 0.01, 0.05, 0.00, 0.00, 0.04};
    Float_t ey3h_0_F[n3] = {0.20, 0.20, 0.30, 0.04, 0.09, 0.14,0.08, 0.03, 0.23, 0.10, 0.20, 0.12, 0.16, 0.23, 0.10, 0.14, 0.41};
    //stat uncertainty
    Float_t ey3l_F[n3] = {0.095, 0.105, 0.125, 0.02, 0.045, 0.07,0.04, 0.015, 0.085, 0.05, 0.07, 0.06, 0.075, 0.09, 0.05, 0.07, 0.185};
    Float_t ey3h_F[n3] = {0.095, 0.105, 0.125, 0.02, 0.045, 0.07,0.04, 0.015, 0.085, 0.05, 0.07, 0.06, 0.075, 0.09, 0.05, 0.07, 0.185};
    // syst uncertainty
    Float_t ey03_s_F[n3] = {0.03, 0.03, 0.03, 0.02, 0.03, 0.02, 0.02, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.04, 0.04, 0.05};
    // total uncertainty
    Float_t ey03l_F[n3] = {0.095, 0.105, 0.125, 0.02, 0.045, 0.07,0.04, 0.015, 0.085, 0.05, 0.07, 0.06, 0.075, 0.09, 0.05, 0.07, 0.185};
    Float_t ey03h_F[n3] = {0.095, 0.105, 0.125, 0.02, 0.045, 0.07,0.04, 0.015, 0.085, 0.05, 0.07, 0.06, 0.075, 0.09, 0.05, 0.07, 0.185};

	for (int iibin = 0; iibin < 17; iibin++) {
		ifstream inputfile(TString::Format("./LHCb_data/%s",LHCb_files[iibin]));
		string s;
		string ss1, ss2;
		int number = 0;
		int number_max = 0;
		float afb_max = -1, afb_min = 1, fh_max = -1, fh_min = 3;
		float a[4000], b[4000], c[4000];
		while (getline( inputfile, s ))
		{
			if (ss2 != "AFB") {
				inputfile >> ss1 >> ss2;
				continue;
			} else {
				inputfile >> a[number] >> b[number] >> c[number];
				number++;
			}
		}
		float pvalue = 0;
		for (int i = 0; i<number; i++) {
			if (pvalue < c[i]) {
				pvalue = c[i];
				number_max = i;
			}
			if (c[i] > 0.6) {
				if (afb_max < a[i]) afb_max = a[i];
				if (afb_min > a[i]) afb_min = a[i];
				if (fh_max < b[i]) fh_max = b[i];
				if (fh_min > b[i]) fh_min = b[i];
			}
		}
////		cout << LHCb_files[iibin] << " " << number_max << ": " << a[number_max] << " " << b[number_max] << " " << c[number_max] << endl;
////		cout << "afb_min = " << afb_min << "; afb_max = " << afb_max << endl;
////		cout << "fh_min  = " << fh_min << "; fh_max  = " << fh_max << endl;
        y3[iibin] = a[number_max];
        ey3l[iibin] = y3[iibin] - ey3l_0[iibin]; if (ey3l[iibin] < 0) cout<<iibin<<" "<<ey3l[iibin]<<endl;
        ey3h[iibin] = ey3h_0[iibin] - y3[iibin]; if (ey3h[iibin] < 0) cout<<iibin<<" "<<ey3h[iibin]<<endl;
        ey03l[iibin] = ey3l[iibin] - ey03_s[iibin];
        ey03h[iibin] = ey3h[iibin] + ey03_s[iibin];
        y3_F[iibin] = b[number_max];
        ey3l_F[iibin] = y3_F[iibin] - ey3l_0_F[iibin]; if (ey3l_F[iibin] < 0) cout<<iibin<<" "<<ey3l_F[iibin]<<endl;
        ey3h_F[iibin] = ey3h_0_F[iibin] - y3_F[iibin]; if (ey3h_F[iibin] < 0) cout<<iibin<<" "<<ey3h_F[iibin]<<endl;
        ey03l_F[iibin] = ey3l_F[iibin] - ey03_s_F[iibin];
        ey03h_F[iibin] = ey3h_F[iibin] + ey03_s_F[iibin];

	}


//	Draw
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitleOffset(1.1,"XY");
	frame->SetTitle("");
	
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////   FH    ////////////////////////////////////////////////////
	frame->SetXTitle("#it{q}^{2} [GeV^{2}]");
	frame->SetYTitle("#it{F_{H}}");
	frame->SetAxisRange(-0.10,1.40,"Y");
////	frame->SetAxisRange(-0.20,3.00,"Y");
	frame->Draw("2");	

////////////////////////////LHCb 1/fb ////////////////////////////////////// 
    TGraphAsymmErrors *grae2_F = new TGraphAsymmErrors(n2,x2_F,y2_F,ex2l_F,ex2h_F,ey2l_F,ey2h_F);
    grae2_F->SetFillColor(2);
    grae2_F->SetMarkerColor(kBlack);
    grae2_F->SetMarkerStyle(24);
    grae2_F->SetMarkerSize(1.);
    // stat uncertainty 
    grae2_F->Draw("p");
    TGraphAsymmErrors *grae02_F  = new TGraphAsymmErrors(n2,x2_F,y2_F,ex2l_F,ex2h_F,ey02l_F,ey02h_F);
    grae02_F->SetFillColor(2);
    grae02_F->SetMarkerColor(kBlack);
    grae02_F->SetMarkerStyle(24);
    grae02_F->SetMarkerSize(1.);
    // total uncertainty  
    grae02_F->Draw("pz");

////////////////////////////LHCb 3/fb ////////////////////////////////////// 
    TGraphAsymmErrors *grae3_F = new TGraphAsymmErrors(n3,x3_F,y3_F,ex3l_F,ex3h_F,ey3l_F,ey3h_F);
    grae3_F->SetMarkerColor(kGray+3);
    grae3_F->SetMarkerStyle(29); // 29
    grae3_F->SetMarkerSize(1.5); // 1.5
    // stat uncertainty    
    grae3_F->Draw("p");
    TGraphAsymmErrors *grae03_F  = new TGraphAsymmErrors(n3,x3_F,y3_F,ex3l_F,ex3h_F,ey03l_F,ey03h_F);
    grae03_F->SetMarkerColor(kGray+3);
    grae03_F->SetMarkerStyle(29); // 29
    grae03_F->SetMarkerSize(1.5); // 1.5
    // total uncertainty   
    grae03_F->Draw("pz"); 

////////////////////////////CMS 20.5/fb ////////////////////////////////////// 
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
    d_fh->SetLineWidth(2);
    d_fh->SetFillColor(1);
    d_fh->SetFillStyle(1001);
    d_fh->SetMarkerColor(1);
    d_fh->SetMarkerStyle(20);
    d_fh->SetMarkerSize(1.5);
    // stat uncertainty    
	d_fh->Draw("p[]");
	TGraphAsymmErrors *sys_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,dataDfh,dataUfh);
    sys_fh->SetLineWidth(2);
    sys_fh->SetFillColor(1);
    sys_fh->SetFillStyle(1001);
    sys_fh->SetMarkerColor(1);
    sys_fh->SetMarkerStyle(20);
	sys_fh->SetMarkerSize(1.5);
    // total uncertainty   
	sys_fh->Draw("pz");

	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();

	TLegend *leg =new TLegend(0.69,0.70,0.90,0.86);
	leg->AddEntry(sys_fh," CMS ","lep");
    leg->AddEntry(grae2_F,"LHCb 1/fb","ple");
    leg->AddEntry(grae3_F,"LHCb 3/fb","ple");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->Draw();
	
	TLatex *tt = new TLatex();
  	tt->SetNDC();
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	TLatex *t2 = new TLatex();
	t2->SetNDC();
	t2->SetTextFont(42);
	t1->DrawLatex(.23,.90,TString::Format("Preliminary"));
	t2->DrawLatex(.68,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	c->Print(TString::Format("./plots/%s_LHCb_fh.pdf",outfile));
	c->Clear();

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////   AFB   ////////////////////////////////////////////////////
	frame->SetXTitle("#it{q}^{2} [GeV^{2}]");
	frame->SetYTitle("#it{A_{FB}}");
	frame->SetAxisRange(-0.4,0.4,"Y");
////	frame->SetAxisRange(-1.0,1.0,"Y");
	frame->Draw();
////////////////////////////LHCb 1/fb ////////////////////////////////////// 
    TGraphAsymmErrors *grae2 = new TGraphAsymmErrors(n2,x2,y2,ex2l,ex2h,ey2l,ey2h);
    grae2->SetFillColor(2);
    grae2->SetMarkerColor(kBlack);
    grae2->SetMarkerStyle(24);
    grae2->SetMarkerSize(1.);
    // stat uncertainty 
    grae2->Draw("p");
    TGraphAsymmErrors *grae02  = new TGraphAsymmErrors(n2,x2,y2,ex2l,ex2h,ey02l,ey02h);
    grae02->SetFillColor(2);
    grae02->SetMarkerColor(kBlack);
    grae02->SetMarkerStyle(24);
    grae02->SetMarkerSize(1.);
    // total uncertainty  
    grae02->Draw("pz");

////////////////////////////LHCb 3/fb ////////////////////////////////////// 
    TGraphAsymmErrors *grae3 = new TGraphAsymmErrors(n3,x3,y3,ex3l,ex3h,ey3l,ey3h);
    grae3->SetFillColor(3);
    grae3->SetMarkerColor(kGray+3);
    grae3->SetMarkerStyle(29); // 29
    grae3->SetMarkerSize(1.5); // 1.5
    // stat uncertainty    
    grae3->Draw("p");
    TGraphAsymmErrors *grae03  = new TGraphAsymmErrors(n3,x3,y3,ex3l,ex3h,ey03l,ey03h);
    grae03->SetFillColor(3);
    grae03->SetMarkerColor(kGray+3);
    grae03->SetMarkerStyle(29); // 29
    grae03->SetMarkerSize(1.5); // 1.5
    // total uncertainty   
    grae03->Draw("pz"); 

////////////////////////////CMS 20.5/fb ////////////////////////////////////// 
	TGraphAsymmErrors *d_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
    d_afb->SetLineWidth(2);
    d_afb->SetFillColor(1);
    d_afb->SetFillStyle(1001);
    d_afb->SetMarkerColor(1);
    d_afb->SetMarkerStyle(20);
    d_afb->SetMarkerSize(1.5);
    // stat uncertainty    
	d_afb->Draw("p[]");
	TGraphAsymmErrors *sys_afb  = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,dataDafb,dataUafb);
    sys_afb->SetLineWidth(2);
    sys_afb->SetFillColor(1);
    sys_afb->SetFillStyle(1001);
    sys_afb->SetMarkerColor(1);
    sys_afb->SetMarkerStyle(20);
	sys_afb->SetMarkerSize(1.5);
    // total uncertainty   
	sys_afb->Draw("pz");

	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2");
   
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	TLegend *leg_1 =new TLegend(0.69,0.70,0.90,0.86);
	leg_1->AddEntry(sys_afb," CMS ","lep");
    leg_1->AddEntry(grae2,"LHCb 1/fb","ple");
    leg_1->AddEntry(grae3,"LHCb 3/fb","ple");
	leg_1->SetLineColor(0);
	leg_1->SetFillColor(0);
	leg_1->SetBorderSize(0);
	leg_1->SetTextSize(0.04);
	leg_1->Draw();
	
	t1->DrawLatex(.23,.90,TString::Format("Preliminary"));
	t2->DrawLatex(.68,.90,TString::Format("20.5 fb^{-1} (8 TeV)"));
  	tt->DrawLatex(.13,.90,TString::Format("CMS"));
	c->Print(TString::Format("./plots/%s_LHCb_afb.pdf",outfile));
}//}}}
