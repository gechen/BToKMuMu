// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
// vim: tw=60 ts=2: 
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
void writeTXT(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
	struct stat fiBuff;
	FILE *fi = 0;
	if (stat(TString::Format("%s/TXTPara%d.txt",odatacardpath.Data(),iBin),&fiBuff) == 0){
		rename(TString::Format("%s/TXTPara%d.txt",odatacardpath.Data(),iBin),TString::Format("%s/TXTPara%d.txt.temp",odatacardpath.Data(),iBin));
		fi = fopen(TString::Format("%s/TXTPara%d.txt.temp",odatacardpath.Data(),iBin),"r");
	}else{
		fi = fopen(TString::Format("%s/TXTPara%d.txt.temp",odatacardpath.Data(),iBin),"w");
	}
	
	bool parExist = false;
	char lineBuff[512];
	char *valBuff = 0;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("%s/TXTPara%d.txt",odatacardpath.Data(),iBin),"w");
	while(fgets(lineBuff,512,fi) != NULL ){
		valBuff = strtok(lineBuff," ");
		if ( strcmp(valBuff,parName) == 0 ){
			fprintf(fp,"%15s",parName);
			int iVal = 0;
			while(iVal < nVal){
				fprintf(fp," %20.15f",val[iVal]);
				iVal++;
			}
			fprintf(fp,"\n");
			parExist = true;
		}else{
			fprintf(fp,"%15s",lineBuff);
			valBuff = strtok(NULL," ");
			while( valBuff != NULL ){
				fprintf(fp," %15s",valBuff);
				valBuff = strtok(NULL," ");
			}
		}
		memset(lineBuff,' ',512*sizeof(char));
	}
	if (parExist == false){
		fprintf(fp,"%15s",parName);
		int iVal = 0;
		while(iVal < nVal){
			fprintf(fp," %20.15f",val[iVal]);
			iVal++;
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	fclose(fi);
	remove(TString::Format("%s/TXTPara%d.txt.temp",odatacardpath.Data(),iBin));
}//}}}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_NLL_bin(int iBin, float Iafb, float Ifh, int Index, const char outfile[] = "angular2D_NLL")
{//{{{
	setTDRStyle();
  idatacardpath="./fitParameters";
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh);
//  fh.setConstant(kTRUE);
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh);
	fh.setConstant(kTRUE);
	RooRealVar afb("afb", "A_{FB}", Iafb);
	afb.setConstant(kTRUE);
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
//	effP0.setConstant(kTRUE);
//	effP1.setConstant(kTRUE);
//	effP2.setConstant(kTRUE);
//	effP3.setConstant(kTRUE);
//	effP4.setConstant(kTRUE);
//	effP5.setConstant(kTRUE);
//	effP6.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
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
//	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
//	} else {
//		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
//	}	
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
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
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
//				bkgCombL_c0.setConstant(kTRUE);
//				bkgCombL_c1.setConstant(kTRUE);
//				bkgCombL_c2.setConstant(kTRUE);
			//	if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
				if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
//				bkgCombL_c3.setConstant(kTRUE);
				}
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

    const double insig[11] = {150,  300,  750,  0, 350,  0, 200, 250, 200,  750, 2200};
    const double inbkg[11] = {500, 1100, 1800,  0, 700,  0, 200, 250, 250, 2300, 4800};
	RooRealVar nsig("nsig","nsig",insig[iBin],0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",inbkg[iBin],0,8E3);
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
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"),  Strategy(2), Warnings(1), PrintEvalErrors(3));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeEffP(iBin, TString::Format("afb_Fix_%d",Index),val, 3);
		val[1]=0; val[2]=0;
		val[0] = fh.getVal();
		writeEffP(iBin, TString::Format("fh_Fix_%d",Index), val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeEffP(iBin, TString::Format("FCN_%d",Index), val);
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); 
		writeOutput(outfile,iBin, Index, "Fix_afb", val);
		val[1]=0; val[2]=0;
		val[0] = Ifh;  val[1] = fh.getVal();
		writeOutput(outfile,iBin, Index, "Fix_fh", val);
		val[1]=0; val[2]=0;
//		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrained
		val[0] = afb.getVal();   // Unconstrained
		writeOutput(outfile,iBin, Index, "F_Fix_afb", val);
//		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrained
		val[0] = fh.getVal();    // Unconstrained
		writeOutput(outfile,iBin, Index, "F_Fix_fh", val); 
		val[0] = f_fitresult->minNll();
		writeOutput(outfile,iBin, Index, "FCN", val);
		val[0] = Index;
		writeOutput(outfile,iBin, Index, "Index", val);
	}
	
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_Profiled_Afb_bin(int iBin, float Iafb, float Ifh, int Index, int idex, const char outfile[] = "angular2D_Profiled_Afb")
{//{{{
	setTDRStyle();
  idatacardpath="./fitParameters";
	cout<<endl<<"iBin = "<<iBin<<"   idex = "<<idex<<endl<<endl; 
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh);
//  fh.setConstant(kTRUE);
//	RooRealVar afb("afb", "A_{FB}", Iafb, -1000, 1000);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh);
	fh.setConstant(kTRUE);
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
//	effP0.setConstant(kTRUE);
//	effP1.setConstant(kTRUE);
//	effP2.setConstant(kTRUE);
//	effP3.setConstant(kTRUE);
//	effP4.setConstant(kTRUE);
//	effP5.setConstant(kTRUE);
//	effP6.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
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
//	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
//	} else {
//		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
//	}	
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
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
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
//				bkgCombL_c0.setConstant(kTRUE);
//				bkgCombL_c1.setConstant(kTRUE);
//				bkgCombL_c2.setConstant(kTRUE);
			//	if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
				if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
//				bkgCombL_c3.setConstant(kTRUE);
				}
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

    const double insig[11] = {150,  300,  750,  0, 350,  0, 200, 250, 200,  750, 2200};
    const double inbkg[11] = {500, 1100, 1800,  0, 700,  0, 200, 250, 250, 2300, 4800};
	RooRealVar nsig("nsig","nsig",insig[iBin],0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",inbkg[iBin],0,8E3);
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
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb)), Strategy(2), Warnings(1), PrintEvalErrors(3));	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3 || f_fitresult->minNll() < -1.E7) {
		if (Index == -1) {
			double val_1[4]={0,0,0,0};
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("afb_Pro_%d",idex),val_1);
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("fh_Fix_%d",idex), val_1);
			val_1[1]= 0;
			writeEffP(iBin, TString::Format("FCN_afb_%d",idex), val_1);
		}
//	if (f_fitresult->status() != 0) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeEffP(iBin, TString::Format("afb_Pro_%d",idex),val, 3);
		val[1]=0; val[2]=0;
		val[0] = fh.getVal();
		writeEffP(iBin, TString::Format("fh_Fix_%d",idex), val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeEffP(iBin, TString::Format("FCN_afb_%d",idex), val);
	} else if (idex == -1) {
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
		writeParam(iBin, TString::Format("afb_Pro_%d",Index-990),val, 3);
		val[1]=0; val[2]=0;
		val[0] = fh.getVal();
		writeParam(iBin, TString::Format("fh_Fix_%d",Index-990), val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, TString::Format("FCN_afb_%d",Index-990), val);
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeEffPOutput(outfile,iBin, idex, Index, "Pro_afb", val);
		val[1]=0; val[2]=0;
		val[0] = Ifh;  val[1] = fh.getVal();
		writeEffPOutput(outfile,iBin, idex, Index, "Fix_fh", val);
		val[1]=0; val[2]=0;
//		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrained
		val[0] = afb.getVal();   // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_Pro_afb", val);
//		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrained
		val[0] = fh.getVal();    // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_Fix_fh", val); 
		val[0] = f_fitresult->minNll();
		writeEffPOutput(outfile,iBin, idex, Index, "FCN_afb", val);
		val[0] = idex;
		writeEffPOutput(outfile,iBin, idex, Index, "idex_afb", val);
	}
	
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_Profiled_Fh_bin(int iBin, float Iafb, float Ifh, int Index, int idex, const char outfile[] = "angular2D_Profiled_Fh")
{//{{{
	setTDRStyle();
  idatacardpath="./fitParameters";
	cout<<endl<<"iBin = "<<iBin<<"   idex = "<<idex<<endl<<endl; 
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	// // Angular parameters
////////////////////////////////////// scan   //////////////////////////////////.........................
//	RooRealVar fh("fh", "F_{H}", Ifh, -1000, 1000 );
//	RooRealVar afb("afb", "A_{FB}", Iafb);
//  afb.setConstant(kTRUE);
////////////////////////////////////// scan   //////////////////////////////////.........................
////////////////////////////////////// refit  //////////////////////////////////.........................
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb);
	afb.setConstant(kTRUE);
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
//	effP0.setConstant(kTRUE);
//	effP1.setConstant(kTRUE);
//	effP2.setConstant(kTRUE);
//	effP3.setConstant(kTRUE);
//	effP4.setConstant(kTRUE);
//	effP5.setConstant(kTRUE);
//	effP6.setConstant(kTRUE);
	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
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
//	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
//	} else {
//		//f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
//		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + afb*CosThetaL )";
//	}	
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
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
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
//				bkgCombL_c0.setConstant(kTRUE);
//				bkgCombL_c1.setConstant(kTRUE);
//				bkgCombL_c2.setConstant(kTRUE);
			//	if (iBin == 1 || iBin == 2 || iBin ==9 || iBin == 10) {
				if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
//				bkgCombL_c3.setConstant(kTRUE);
				}
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

    const double insig[11] = {150,  300,  750,  0, 350,  0, 200, 250, 200,  750, 2200};
    const double inbkg[11] = {500, 1100, 1800,  0, 700,  0, 200, 250, 250, 2300, 4800};
	RooRealVar nsig("nsig","nsig",insig[iBin],0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",inbkg[iBin],0,8E3);
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
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(fh)), Strategy(2), Warnings(1), PrintEvalErrors(3));	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Strategy(2), Warnings(1), PrintEvalErrors(3));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3 || f_fitresult->minNll() < -1.E7) {
		if (Index == -1) {
			double val_1[4]={0,0,0,0};
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("afb_Fix_%d",idex),val_1);
			val_1[0] = -999;
			writeEffP(iBin, TString::Format("fh_Pro_%d",idex), val_1);
			val_1[1]= 0;
			writeEffP(iBin, TString::Format("FCN_fh_%d",idex), val_1);
		}
//	if (f_fitresult->status() != 0) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}
	delete data;
// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = afb.getVal();
		writeEffP(iBin, TString::Format("afb_Fix_%d",idex),val, 3);
		val[1]=0; val[2]=0;
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeEffP(iBin, TString::Format("fh_Pro_%d",idex), val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeEffP(iBin, TString::Format("FCN_fh_%d",idex), val);
	} else if (idex == -1) {
		val[0] = afb.getVal();
		writeParam(iBin, TString::Format("afb_Fix_%d",Index-990),val, 3);
		val[1]=0; val[2]=0;
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
		writeParam(iBin, TString::Format("fh_Pro_%d",Index-990), val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, TString::Format("FCN_fh_%d",Index-990), val);
	} else {
		val[0] = Iafb; val[1] = afb.getVal();
		writeEffPOutput(outfile,iBin, idex, Index, "Fix_afb", val);
		val[1]=0; val[2]=0;
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeEffPOutput(outfile,iBin, idex, Index, "Pro_fh", val);
		val[1]=0; val[2]=0;
	//	val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrained
		val[0] = afb.getVal();   // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_Fix_afb", val);
//		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrained
		val[0] = fh.getVal();   // Unconstrained
		writeEffPOutput(outfile,iBin, idex, Index, "F_Pro_fh", val); 
		val[0] = f_fitresult->minNll();
		writeEffPOutput(outfile,iBin, idex, Index, "FCN_fh", val);
		val[0] = idex;
		writeEffPOutput(outfile,iBin, idex, Index, "idex_fh", val);
	}
	
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotFCN_Profiled_Afb( int iBin, const char outfile[] = "FCN")
{
    setTDRStyle();
    TGraph *gr_afb = new TGraph();
    int n = -1;
    for (int idex = 0; idex < 100; idex++) {	
        double FCN = 999, fcn = 999;
        int Index = 0, index = 0, NIndex = 1000;
        double Iafb, Ifh, afb, fh;
        do {
            index+=1;
            while ( ! (fcn = readEffPOutput(outfile, iBin, idex, index, "FCN_afb", 0)) && index < NIndex) index+=1;
            afb = readEffPOutput(outfile, iBin, idex, index, "F_Pro_afb", 0);  
            fh  = readEffPOutput(outfile, iBin, idex, index, "F_Fix_fh", 0);
//	          if (afb < 0) cout<<idex<<" "<<afb<<endl; 
//	          if (fcn < 0 && (afberr > 0.001) && (fherr > 0.001) ) {
            if ( FCN > fcn ) { FCN = fcn; Index = index; }
//	          if ( FCN > fcn && (fabs(afb) < (fh / 2.))) { FCN = fcn; Index = index; }
        } while ( index < NIndex);
        Iafb   = readEffPOutput(outfile, iBin, idex, Index, "Pro_afb", 0);
        Ifh   = readEffPOutput(outfile, iBin, idex, Index, "Fix_fh", 0);
        afb    = readEffPOutput(outfile, iBin, idex, Index, "F_Pro_afb", 0);  
        fh    = readEffPOutput(outfile, iBin, idex, Index, "F_Fix_fh", 0);  
//        cout<<"rm -rf OutputValues/angular2D_Profiled_Afb/bin2/OutputValues2_EffP"<<idex<<"_Index"<<Index<<".txt"<<endl;
        cout<<"fh_index_"<<idex<<endl;
        cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
        cout<<"afb = "<<afb<<endl;
        cout<<"fh  = "<<fh<<endl;
        double val[3]={0,0,0};
        val[0] = Iafb; val[1] = readEffPOutput(outfile, iBin, idex, Index, "Pro_afb", 1);
        writeEffP(iBin, TString::Format("Iafb_%s_Pro_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = Ifh; val[1] = readEffPOutput(outfile, iBin, idex, Index, "Fix_fh", 1);
        writeEffP(iBin, TString::Format("Ifh_%s_Fix_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = afb;
        writeEffP(iBin, TString::Format("afb_%s_Pro_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = fh;
        writeEffP(iBin, TString::Format("fh_%s_Fix_%d",outfile,idex),val);
        if(Index == 1000) continue;  
//        if(iBin == 4 || iBin == 7 || iBin == 8 || iBin == 0 || iBin == 10 ) { if ( afb < 0 && fh < 0.1) continue; } 
        n+=1;
        gr_afb->SetPoint(n, fh, afb);
    }
    const double upAfb[11] = { 0.12, 0.02, 0.15, 0.00, 0.03, 0.00, 0.02, 0.06, 0.10, 0.10, 0.03};
    const double dnAfb[11] = {-0.01,-0.08,-0.15, 0.00,-0.03, 0.00,-0.01,-0.02,-0.10,-0.30,-0.01};
    const double upFh[11]  = { 0.60, 1.00, 0.30, 0.00, 0.10, 0.00, 0.15, 0.30, 0.30, 0.60, 0.30};
    TCanvas *c = new TCanvas("c","c",800,600);
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    TLine *line = new TLine();
//    c->SetTitle("FCN distribution of A_{FB}");
    gr_afb->GetYaxis()->SetTitle("Profiled A_{FB}");
    gr_afb->GetXaxis()->SetTitle("True F_{H}");
    gr_afb->GetXaxis()->SetLimits(0.,upFh[iBin]);
    gr_afb->GetYaxis()->SetRangeUser(dnAfb[iBin], upAfb[iBin]);
//    gr_afb->GetXaxis()->SetLimits(0.,0.30);
    gr_afb->Draw("AP");
    line->SetLineStyle(1);
    line->SetLineColor(2);
//    line->DrawLine(0,0,0.40,0.20);
    line->DrawLine(0,0, 2*upAfb[iBin],upAfb[iBin]);
    line->DrawLine(0,0,-2*dnAfb[iBin],dnAfb[iBin]);
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->SetTextFont(12);
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
//    TArrow *ar2 = new TArrow(0.30,0.09,0.50,0.09,0.05,"|>");
    TArrow *ar2 = new TArrow(upFh[iBin]*0.70,upAfb[iBin]-0.01,upFh[iBin]*0.9,upAfb[iBin]-0.01,0.05,"|>");
    ar2->SetAngle(60);
    ar2->SetLineWidth(8);
    ar2->SetLineColor(3);
    ar2->SetFillColor(3);
    ar2->Draw();
    const double upg1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    const double upg2[11] = { 0.600, 0.600, 0.300, 0.00, 0.100, 0.00, 0.150, 0.200, 0.250, 0.600, 0.250};
    TF1 *g1    = new TF1("g1","[0]+[1]*x", 0.00, upg1[iBin]);
    TF1 *g2    = new TF1("g2","[0]+[1]*x", upg1[iBin], upg2[iBin]);
    g1->SetLineColor(4);
    g2->SetLineColor(4);
    gr_afb->Fit(g1,"RS","",0.00,upg1[iBin]);
    gr_afb->Fit(g2,"RS","",upg1[iBin],upg2[iBin]);
    g1->Draw("SAME");
    g2->Draw("SAME");
    Double_t par[4];
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[2]);
    writeParam(iBin,"Pro_Afb_shape", par, 4);
    c->Print(TString::Format("./plots/%s_Profiled_afb_bin%d.pdf",outfile,iBin));
    delete c;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotFCN_Profiled_Fh( int iBin, const char outfile[] = "FCN")
{
    setTDRStyle();
    TGraph *gr_afb = new TGraph();
    int n = -1;
    int nnbins = 100; 
    if (iBin == 0 ) nnbins = 200;
    for (int idex = 0; idex < nnbins; idex++) {	
        double FCN = 999, fcn = 999;
        int Index = 0, index = 0, NIndex = 1000;
        double Iafb, Ifh, afb, fh;
        do {
            index+=1;
            while ( ! (fcn = readEffPOutput(outfile, iBin, idex, index, "FCN_fh", 0)) && index < NIndex) index+=1;
            afb = readEffPOutput(outfile, iBin, idex, index, "F_Fix_afb", 0);  
            fh  = readEffPOutput(outfile, iBin, idex, index, "F_Pro_fh", 0);
//	          if (fcn < 0 && (afberr > 0.001) && (fherr > 0.001) ) {
            if ( FCN > fcn ) { FCN = fcn; Index = index; }
//	          if ( FCN > fcn && (fabs(afb) < (fh / 2.))) { FCN = fcn; Index = index; }
        } while ( index < NIndex);
        Iafb   = readEffPOutput(outfile, iBin, idex, Index, "Fix_afb", 0);
        Ifh   = readEffPOutput(outfile, iBin, idex, Index, "Pro_fh", 0);
        afb    = readEffPOutput(outfile, iBin, idex, Index, "F_Fix_afb", 0);  
        fh    = readEffPOutput(outfile, iBin, idex, Index, "F_Pro_fh", 0);  
//        cout<<"rm -rf OutputValues/angular2D_Profiled_Afb/bin8/OutputValues8_EffP"<<idex<<"_Index"<<Index<<".txt"<<endl;
        cout<<"afb_index_"<<idex<<endl;
        cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
        cout<<"afb = "<<afb<<endl;
        cout<<"fh  = "<<fh<<endl;
        double val[3]={0,0,0};
        val[0] = Iafb; val[1] = readEffPOutput(outfile, iBin, idex, Index, "Fix_afb", 1);
        writeEffP(iBin, TString::Format("Iafb_%s_Fix_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = Ifh; val[1] = readEffPOutput(outfile, iBin, idex, Index, "Pro_fh", 1);
        writeEffP(iBin, TString::Format("Ifh_%s_Pro_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = afb;
        writeEffP(iBin, TString::Format("afb_%s_Fix_%d",outfile,idex),val);
        val[1]=0; val[2]=0;
        val[0] = fh;
        writeEffP(iBin, TString::Format("fh_%s_Pro_%d",outfile,idex),val);
        if(Index == 1000) continue;  
        n+=1;
        gr_afb->SetPoint(n, afb, fh);
    }
    const double upAfb[11] = { 0.20, 0.60, 0.15, 0.00, 0.05, 0.00, 0.10, 0.20, 0.20, 0.40, 0.10};
    const double dnAfb[11] = {-0.20,-0.60,-0.15, 0.00,-0.05, 0.00,-0.10,-0.20,-0.20,-0.40,-0.10};
    const double upFh[11]  = { 0.60, 2.00, 0.30, 0.00, 0.15, 0.00, 0.20, 0.40, 0.40, 1.00, 0.20};
    TCanvas *c = new TCanvas("c","c",800,600);
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    TLine *line = new TLine();
//    c->SetTitle("FCN distribution of A_{FB}");
    gr_afb->GetYaxis()->SetTitle("Profiled F_{H}");
    gr_afb->GetXaxis()->SetTitle("True A_{FB}");
    gr_afb->GetXaxis()->SetLimits(dnAfb[iBin], upAfb[iBin]);
    gr_afb->GetYaxis()->SetRangeUser(0.,upFh[iBin]);
    gr_afb->Draw("AP");
    line->SetLineStyle(1);
    line->SetLineColor(2);
    line->DrawLine(0,0, upAfb[iBin], 2*upAfb[iBin]);
    line->DrawLine(0,0, dnAfb[iBin],-2*dnAfb[iBin]);
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->SetTextFont(12);
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    TArrow *ar2 = new TArrow(upAfb[iBin]-0.01,upFh[iBin]*0.70,upAfb[iBin]-0.01,upFh[iBin]*0.9,0.05,"|>");
    ar2->SetAngle(60);
    ar2->SetLineWidth(8);
    ar2->SetLineColor(3);
    ar2->SetFillColor(3);
    ar2->Draw();
    const double upg0[11] = {-0.200,-0.600,-0.100, 0.00,-0.050, 0.00,-0.100,-0.150,-0.200,-0.400,-0.100};
    const double upg1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double upg2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    const double upg3[11] = { 0.200, 0.600, 0.100, 0.00, 0.050, 0.00, 0.100, 0.150, 0.200, 0.400, 0.100};
    TF1 *g1    = new TF1("g1","[0]+[1]*x", upg0[iBin], upg1[iBin]);
    TF1 *g2    = new TF1("g2","[0]+[1]*x+[2]*x**2", upg1[iBin], upg2[iBin]);
    TF1 *g3    = new TF1("g3","[0]+[1]*x", upg2[iBin], upg3[iBin]);
    g1->SetLineColor(4);
    g2->SetLineColor(4);
    g3->SetLineColor(4);
    gr_afb->Fit(g1,"RS","",upg0[iBin],upg1[iBin]);
    gr_afb->Fit(g2,"RS","",upg1[iBin],upg2[iBin]);
    gr_afb->Fit(g3,"RS","",upg2[iBin],upg3[iBin]);
    g1->Draw("SAME");
    g2->Draw("SAME");
    g3->Draw("SAME");
    Double_t par[7];
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[2]);
    g3->GetParameters(&par[5]);
    writeParam(iBin,"Pro_Fh_shape", par, 7);
    c->Print(TString::Format("./plots/%s_Profiled_fh_bin%d.pdf",outfile,iBin));
    delete c;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotAfbFh_NLL( int iBin, int Index, const char outfile[] = "AfbFh2D")
{
  	setTDRStyle();
    const double stepSizeAfb[11] = { 0.0200,  0.0300, 0.0060, 0.00, 0.0060, 0.00, 0.0060, 0.0040, 0.0080, 0.0200, 0.0060};
    const double stepSizeFh[11]  = { 0.0200,  0.0200, 0.0100, 0.00, 0.0060, 0.00, 0.0080, 0.0100, 0.0100, 0.0200, 0.0080};
    const double dataAfbLo[11]  = {-0.550,-0.550,-0.025, 0.00,-0.025, 0.00,-0.150,-0.150,-0.150,-0.350,-0.025};
    const double dataFhBand[11] = { 1.20,   1.20,  0.15, 0.00,  0.15, 0.00,  0.40,  0.40,  0.40,  0.80,  0.15};
//    const double ToysAfb[5] = { -0.0050, -0.0070,  0.0010, 0.0398, 0.0700 };
//    const double ToysFh[5]  = {  0.0200,  0.0706,  0.1700, 0.1900, 0.1650 };
//  	TGraph *gr_cov = new TGraph(5, ToysAfb, ToysFh);
    const double ToysAfb[2] = { -0.0070,  0.0398 };
    const double ToysFh[2]  = {  0.0706,  0.1900 };
  	TGraph *gr_cov = new TGraph(2, ToysAfb, ToysFh);
  	TGraphAsymmErrors *gr_afb = new TGraphAsymmErrors();
  	TGraphAsymmErrors *gr_fh  = new TGraphAsymmErrors();
    double  fh  = readParam(iBin,"fh", 0);
    double  afb = readParam(iBin,"afb", 0); 
    int iiafb=0, iifh=0;
    // FH 
    const double dnA1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    double par_Afb[4]; 
    par_Afb[0] = readParam(iBin,"Pro_Afb_shape", 0);
    par_Afb[1] = readParam(iBin,"Pro_Afb_shape", 1);
    par_Afb[2] = readParam(iBin,"Pro_Afb_shape", 2);
    par_Afb[3] = readParam(iBin,"Pro_Afb_shape", 3);
    // Loop over phase space
    double thisAfb = afb;
    double thisFh  = fh;
    for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
        thisFh = iFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        gr_fh->SetPoint(iifh, thisAfb, thisFh);
        iifh++;
    }
    // AFB 
    const double dnF1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double dnF2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    double par_Fh[7]; 
    par_Fh[0] = readParam(iBin,"Pro_Fh_shape", 0);
    par_Fh[1] = readParam(iBin,"Pro_Fh_shape", 1);
    par_Fh[2] = readParam(iBin,"Pro_Fh_shape", 2);
    par_Fh[3] = readParam(iBin,"Pro_Fh_shape", 3);
    par_Fh[4] = readParam(iBin,"Pro_Fh_shape", 4);
    par_Fh[5] = readParam(iBin,"Pro_Fh_shape", 5);
    par_Fh[6] = readParam(iBin,"Pro_Fh_shape", 6);
    thisAfb = afb;
    thisFh  = fh;
    cout<<"  ---    1  ---"<<endl;
    for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        gr_afb->SetPoint(iiafb, thisAfb, thisFh);
        iiafb++;
    }
//  	TGraph2D *gr   = new TGraph2D("gr","#Delta NLL distribution for data; A_{FB}; F_{H}; #Delta NLL");
  	TGraph2D *gr   = new TGraph2D();
  	TGraphAsymmErrors *gr_0 = new TGraphAsymmErrors();
  	TGraphAsymmErrors *gr_1 = new TGraphAsymmErrors();
  	TGraphAsymmErrors *gr_2 = new TGraphAsymmErrors();
  	TGraphAsymmErrors *gr_3 = new TGraphAsymmErrors();
  	TGraphAsymmErrors *gr_data = new TGraphAsymmErrors();
  	TCanvas *c = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    double outFCErrFh[2], outFCErrAfb[2];
    outFCErrFh[0]  = readParam(iBin,"FCErrFh" , 0);
    outFCErrFh[1]  = readParam(iBin,"FCErrFh" , 1);
    outFCErrAfb[0] = readParam(iBin,"FCErrAfb", 0);
    outFCErrAfb[1] = readParam(iBin,"FCErrAfb", 1);
  	gr_data->SetPoint(0, afb, fh);
//		gr_data->SetPointError(0, fabs(outFCErrAfb[0]), fabs(outFCErrAfb[1]), fabs(outFCErrFh[0]), fabs(outFCErrFh[1]) );
	
  	double fcn = 999, dataFCN = 999;
    dataFCN = readParam(iBin,"FCN", 0); 
  	int index = 0, n = -1, NIndex = 2000;
    int n0=-1, n1=-1, n2=-1, n3=-1;
  	double aa, bb;
  	for (index=0; index<NIndex; index++) {
  			if ( !(fcn = readOutput(outfile, iBin, index, "FCN", 0)) ) continue;
  			n+=1;
  			aa = readOutput(outfile, iBin, index, "F_Fix_afb", 0);  // best fitted values
  			bb  = readOutput(outfile, iBin, index, "F_Fix_fh", 0);
  			if (fcn < 0 ) {
  				  gr->SetPoint(n, aa, bb, fcn-dataFCN);
            if (fcn-dataFCN < 0) { n0+=1; gr_0->SetPoint(n0, aa, bb); cout<<fcn-dataFCN<<endl;}
            else if (fcn-dataFCN < 0.5) { n1+=1; gr_1->SetPoint(n1, aa, bb); }
            else if (fcn-dataFCN < 2.0) { n2+=1; gr_2->SetPoint(n2, aa, bb); } 
            else { n3+=1; gr_3->SetPoint(n3, aa, bb); } 
//            cout<<fcn-dataFCN<<endl;
//            if (bb<0.015) cout<<index<<"  "<<n<<endl;
  			}
    }
  	gr->Draw("surf1");	
    c->Print(TString::Format("./plots/%s_afb_fh_2D_bin%d.pdf",outfile,iBin));  // fitted
  	c->Clear();

  	gr_1->SetMarkerColor(394);
  	gr_1->GetXaxis()->SetTitle("A_{FB}");
  	gr_1->GetYaxis()->SetTitle("F_{H}");
  	gr_1->GetXaxis()->SetLimits(  dataAfbLo[iBin]-0.05, -dataAfbLo[iBin]+0.05);  
  	gr_1->GetYaxis()->SetRangeUser(0.00, dataFhBand[iBin]);
  	gr_1->Draw("AP");	
//  	gr_0->SetMarkerColor(1);
//  	gr_0->Draw("P SAME");	
  	gr_2->SetMarkerColor(8);
  	gr_2->Draw("P SAME");	
  	gr_3->SetMarkerColor(886);
  	gr_3->Draw("P SAME");	

  	gr_fh->SetMarkerStyle(34);
  	gr_fh->SetMarkerColor(2);
//  	gr_fh->Draw("P SAME");	
  	gr_afb->SetMarkerStyle(34);
  	gr_afb->SetMarkerColor(890);
//  	gr_afb->Draw("P SAME");	
  	gr_data->SetMarkerStyle(29);
  	gr_data->SetMarkerColor(1);
  	gr_data->SetMarkerSize(2);
  	gr_data->SetLineWidth(2);
  	gr_data->Draw("P SAME");	
  	gr_cov->SetMarkerStyle(21);
  	gr_cov->SetMarkerColor(1);
//  	gr_cov->Draw("P SAME");	

    TLine *line = new TLine();
    line->SetLineStyle(1);
    line->SetLineColor(2);
    line->DrawLine(0,0,  dataAfbLo[iBin]-0.05, dataFhBand[iBin]);
    line->DrawLine(0,0, -dataAfbLo[iBin]+0.05, dataFhBand[iBin]);
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->SetTextFont(12);
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    latex->DrawLatexNDC(0.15,0.28,TString::Format("A_{FB}=%.4f",afb).Data());
    latex->DrawLatexNDC(0.15,0.18,TString::Format("F_{H}  =%.4f",fh).Data());
//    latex->DrawLatexNDC(0.15,0.28,TString::Format("A_{FB}=%.4f^{%+.4f}_{%+.4f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
//    latex->DrawLatexNDC(0.15,0.18,TString::Format("F_{H}  =%.4f^{%+.4f}_{%+.4f}",fh,outFCErrFh[1],outFCErrFh[0]).Data());
//		TLegend *leg =new TLegend(0.16,0.75,0.36,0.88,NULL,"brNDC");
		TLegend *leg =new TLegend(0.69,0.15,0.89,0.40,NULL,"brNDC");
//		leg->AddEntry(gr_0," #Delta NLL < 0","P");
		leg->AddEntry(gr_1," 0.0 < #Delta NLL < 0.5","P");
		leg->AddEntry(gr_2," 0.5 < #Delta NLL < 2.0","P");
		leg->AddEntry(gr_3," #Delta NLL > 2.0","P");
		leg->AddEntry(gr_fh," Profiled A_{FB}","P");
		leg->AddEntry(gr_afb," Profiled F_{H}","P");
		leg->AddEntry(gr_data," Best Fit","P");
		leg->SetLineColor(1);
		leg->SetFillColor(0);
		leg->SetTextSize(0.02);
		leg->Draw();					
    c->Print(TString::Format("./plots/%s_afb_fh_bin%d.pdf",outfile,iBin));  // fitted
  	c->Clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D_vali( const char outfile[] = "angular2D")
{//{{{
	setTDRStyle();
	double conX[2]    = { 9.385, 13.52};
	double conXerr[2] = { 0.705,  0.66};
	double conYA[2]   = { 0.00,   0.00};
	double conYAerr[2]= { 0.50,   0.50};
	double conYF[2]   = { 0.70,   0.70};
	double conYFerr[2]= { 0.80,   0.80};

	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={ 0.5, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
  double aa[9]  ={0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10};
  double aa1[9] ={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  double aa2[9] ={-0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30, -0.30};
  double ff[9]  ={0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40};
  double ff1[9] ={0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  double ff2[9] ={0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80};
	double yfh[9],  yuerrfh[9], yderrfh[9],  yafb[9],  yuerrafb[9],  yderrafb[9];
	double yfh1[9], yuerrfh1[9],yderrfh1[9], yafb1[9], yuerrafb1[9], yderrafb1[9];
	double yfh2[9], yuerrfh2[9],yderrfh2[9], yafb2[9], yuerrafb2[9], yderrafb2[9];
	for (int i = 0; i < 9; i++) {
		yfh[i]  = -2; yuerrfh[i] = -2; yderrfh[i] = -2;
		yafb[i] = -2; yuerrafb[i]= -2; yderrafb[i]= -2;
		yfh1[i]  = -2; yuerrfh1[i] = -2; yderrfh1[i] = -2;
		yafb1[i] = -2; yuerrafb1[i]= -2; yderrafb1[i]= -2;
		yfh2[i]  = -2; yuerrfh2[i] = -2; yderrfh2[i] = -2;
		yafb2[i] = -2; yuerrafb2[i]= -2; yderrafb2[i]= -2;
	}
// Checkout input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]      = readParam(ibin,"afb_v0",0);
		yuerrafb[i]  = readParam(ibin,"afb_v0",1);
		yderrafb[i]  = -readParam(ibin,"afb_v0",2);
		yfh[i]       = readParam(ibin,"fh_v0",0);
		yuerrfh[i]   = readParam(ibin,"fh_v0",1);
		yderrfh[i]   = -readParam(ibin,"fh_v0",2);
		//printf("Afb0[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		//printf("Fh0 [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i], yuerrfh[i], yderrfh[i]);
		yafb1[i]      = readParam(ibin,"afb_v1",0);
		yuerrafb1[i]  = readParam(ibin,"afb_v1",1);
		yderrafb1[i]  = -readParam(ibin,"afb_v1",2);
		yfh1[i]       = readParam(ibin,"fh_v1",0);
		yuerrfh1[i]   = readParam(ibin,"fh_v1",1);
		yderrfh1[i]   = -readParam(ibin,"fh_v1",2);
    if (ibin == 0) {
       yuerrafb1[i] = 0;
       yderrafb1[i] = 0;
       yuerrfh1[i] = 0;
       yderrfh1[i] = 0;
    }
		//printf("Afb1[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb1[i],yuerrafb1[i],yderrafb1[i]);
		//printf("Fh1 [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh1[i], yuerrfh1[i], yderrfh1[i]);
		yafb2[i]      = readParam(ibin,"afb_v2",0);
		yuerrafb2[i]  = readParam(ibin,"afb_v2",1);
		yderrafb2[i]  = -readParam(ibin,"afb_v2",2);
		yfh2[i]       = readParam(ibin,"fh_v2",0);
		yuerrfh2[i]   = readParam(ibin,"fh_v2",1);
		yderrfh2[i]   = -readParam(ibin,"fh_v2",2);
    if (ibin == 0) {
       yafb2[i]     = -2;
       yfh2[i]      = -2;
       yuerrafb2[i] = 0;
       yderrafb2[i] = 0;
       yuerrfh2[i]  = 0;
       yderrfh2[i]  = 0;
    }
		//printf("Afb2[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb2[i],yuerrafb2[i],yderrafb2[i]);
		//printf("Fh2 [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh2[i], yuerrfh2[i], yderrfh2[i]);
	}
//	Draw
	TCanvas *c = new TCanvas();
	TH1F *frame = new TH1F("frame","",22,0.,22);
	frame->SetStats(kFALSE);
	frame->SetTitleOffset(1.1,"Y");
	frame->SetTitle("");
	frame->Draw("2");	
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("F_{H}");
	frame->SetAxisRange(-0.1,1.2,"Y");
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	d_fh->SetMarkerColor(3);
	d_fh->SetLineColor(3);
	d_fh->SetMarkerStyle(24);
	d_fh->Draw("P");
	TGraphAsymmErrors *d_fh1  = new TGraphAsymmErrors(7,x,yfh1,xerr,xerr,yderrfh1,yuerrfh1);
	d_fh1->SetMarkerColor(2);
	d_fh1->SetLineColor(2);
	d_fh1->SetMarkerStyle(25);
	d_fh1->Draw("P");
	TGraphAsymmErrors *d_fh2  = new TGraphAsymmErrors(7,x,yfh2,xerr,xerr,yderrfh2,yuerrfh2);
	d_fh2->SetMarkerColor(4);
	d_fh2->SetLineColor(4);
	d_fh2->SetMarkerStyle(30);
	d_fh2->Draw("P");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();

	TGraph *f_fh  = new TGraph(7,x,ff);
	f_fh->SetMarkerColor(3);
  f_fh->Draw("P");
	TGraph *f1_fh  = new TGraph(7,x,ff1);
	f1_fh->SetMarkerColor(2);
  f1_fh->Draw("P");
	TGraph *f2_fh  = new TGraph(7,x,ff2);
	f2_fh->SetMarkerColor(4);
  f2_fh->Draw("P");

	TLegend *leg =new TLegend(0.69,0.69,0.90,0.86);
	leg->AddEntry(d_fh, " Input: F_{H} = 0.40 ","lep");
	leg->AddEntry(d_fh1," Input: F_{H} = 0.05 ","lep");
	leg->AddEntry(d_fh2," Input: F_{H} = 0.80 ","lep");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.02);
	leg->Draw();
	
	TLatex *t1 = new TLatex();
	t1->SetNDC();
	t1->SetTextFont(12);
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.62,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/%s_fh.pdf",outfile));
	c->Clear();
	
	frame->SetXTitle("q^{2} [(GeV)^{2}]");
	frame->SetYTitle("A_{FB}");
	frame->SetAxisRange(-0.5,0.5,"Y");
	frame->Draw();
	TGraphAsymmErrors *d_afb  = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	d_afb->SetMarkerColor(3);
	d_afb->SetLineColor(3);
	d_afb->SetMarkerStyle(24);
	d_afb->Draw("P");
	TGraphAsymmErrors *d_afb1  = new TGraphAsymmErrors(7,x,yafb1,xerr,xerr,yderrafb1,yuerrafb1);
	d_afb1->SetMarkerColor(2);
	d_afb1->SetLineColor(2);
	d_afb1->SetMarkerStyle(25);
	d_afb1->Draw("P");
	TGraphAsymmErrors *d_afb2  = new TGraphAsymmErrors(7,x,yafb2,xerr,xerr,yderrafb2,yuerrafb2);
	d_afb2->SetMarkerColor(4);
	d_afb2->SetLineColor(4);
	d_afb2->SetMarkerStyle(30);
	d_afb2->Draw("P");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2");

	TGraph *f_afb  = new TGraph(7,x,aa);
	f_afb->SetMarkerColor(3);
  f_afb->Draw("P");
	TGraph *f1_afb  = new TGraph(7,x,aa1);
	f1_afb->SetMarkerColor(2);
  f1_afb->Draw("P");
	TGraph *f2_afb  = new TGraph(7,x,aa2);
	f2_afb->SetMarkerColor(4);
  f2_afb->Draw("P");
   
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	TLegend *leg_1 =new TLegend(0.69,0.69,0.90,0.86);
	leg_1->AddEntry(d_afb, " Input: A_{FB} = 0.10 ","lep");
	leg_1->AddEntry(d_afb1," Input: A_{FB} = 0.01 ","lep");
	leg_1->AddEntry(d_afb2," Input: A_{FB} =-0.30 ","lep");
	leg_1->SetLineColor(0);
	leg_1->SetFillColor(0);
	leg_1->SetTextSize(0.02);
	leg_1->Draw();
	
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.62,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_Toy_bin(int iBin, const char outfile[] = "angular2D_Toy")
{//{{{
	setTDRStyle();
	printf("INFO\t\t: Processing angular2D_Toy_bin(iBin=%d)\n",iBin);
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
	double iafb = readParam(iBin, "Iafb_angular2D", 0);  // Selected initial values!
	double ifh  = readParam(iBin, "Ifh_angular2D", 0);
	double Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
	double Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
	// // Angular parameters
	RooRealVar fh("fh", "F_{H}", Ifh, 0., 3. );
	RooRealVar afb("afb", "A_{FB}", Iafb, -1., 1.);
		
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
	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10)); 
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
	effP10.setConstant(kTRUE);
/*	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
/////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
*/
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
	// // Signal double gaussian
	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.26,5.30);
	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
	sigGauss1_sigma.setConstant(kTRUE);
	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
	sigGauss2_sigma.setConstant(kTRUE);
	RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
	sigM_frac.setConstant(kTRUE);
	// // mass distro of signal
	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
	RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
	
	RooArgSet f_sigA_argset(CosThetaL);
	f_sigA_argset.add(RooArgSet(fh,afb));
	TString f_sigA_format;
	TString f_rec_format;
	TString f_ang_format;
  f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	if (iBin != 0 && iBin != 1 && iBin != 9) {
		f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
		f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
	} else {
		f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
	}
	f_sigA_argset.add(RooArgSet(f_rec_format));
	f_sigA_argset.add(RooArgSet(f_ang_format));
	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
	RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);
	// Create signal distribution
	RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_sig prepared. <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	
	// Create combinatorial background distribution
	RooRealVar bkgCombM_c("bkgCombM_c","c",-0.5,-20.,1.);
	RooRealVar offset("offset","offset",-5.);
	RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
	RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
	
	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",readParam(iBin,"bkgCombL_c0",0));
	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0));
	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0));
	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0)); 
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
		default:
		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
				bkgCombL_c0.setConstant(kTRUE);
				bkgCombL_c1.setConstant(kTRUE);
				bkgCombL_c2.setConstant(kTRUE);
				if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
				f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
				bkgCombL_c3.setConstant(kTRUE);
				}
				break;
	}
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	RooRealVar bkgGauss_mean("bkgGauss_mean","cos#theta_{l}", readParam(iBin,"bkgGauss_mean",0));
	bkgGauss_mean.setConstant(kTRUE);
   RooRealVar bkgGauss_sigma("bkgGauss_sigma","#sigma",  readParam(iBin,"bkgGauss_sigma",0));
	bkgGauss_sigma.setConstant(kTRUE);
   RooRealVar bkg_frac("bkg_frac","bkg_frac",readParam(iBin,"bkg_frac",0));
	bkg_frac.setConstant(kTRUE);
	
	RooGaussian f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", CosThetaL, bkgGauss_mean, bkgGauss_sigma);
	RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombL_P), bkg_frac);
	
	RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL, f_bkgCombM);
	cout<<">>>>>>>>>>>>>>>> INFO: f_bkgComb prepared. <<<<<<<<<<<<<<<<<<<<<<"<<endl;
	
  RooRealVar nsig("nsig","nsig",50,0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",100,0,6E3);
	
	RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_sig),RooArgList(nbkgComb,nsig));// no penalty term
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	

	// Get data and apply unbinned fit
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	
  f_fitresult->Print();
  delete data;
  odatacardpath=summarypath;
  double val[2]={0,0};
  val[0] = f_fitresult->status();
  writeParam(iBin,  "migrad",  val, 1);
  val[0] = f_fitresult->minNll();
  writeParam(iBin, "FC_FCN", val);
  val[0] = fh.getVal(); val[1] = fh.getError();
  writeParam(iBin,  "FC_fh",  val, 2);
  val[0] = afb.getVal(); val[1] = afb.getError();
  writeParam(iBin,  "FC_afb",  val, 2);
/*  // map of signal
  bkgCombM_c.setConstant(kTRUE);
  offset.setConstant(kTRUE);
  fh.setConstant(kTRUE);
  afb.setConstant(kTRUE);
  nsig.setConstant(kTRUE);
  nbkgComb.setConstant(kTRUE);
  sigGauss_mean.setConstant(kTRUE);
  RooWorkspace *wspace2 = new RooWorkspace("wspace","wspace");
  wspace2->import(f);
  wspace2->writeToFile(TString::Format("%s/wspace_FC_bin%d.root",summarypath.Data(),iBin),true);
*/
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_Toy_unCons_bin(int iBin, const char outfile[] = "angular2D_Toy_unCons")
{//{{{
  	setTDRStyle();
  	printf("INFO\t\t: Processing angular2D_Toy_unCons_bin(iBin=%d)\n",iBin);
  // Create parameters and PDFs
  	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
  	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
  	RooRealVar Q2("Q2","q^{2}",1.0,22.);
  	double iafb = readParam(iBin, "Iafb_angular2D", 0);  // Selected initial values!
  	double ifh  = readParam(iBin, "Ifh_angular2D", 0);
  	// // Angular parameters
  	RooRealVar fh("fh", "F_{H}", ifh, -1000., 1000. );
  	RooRealVar afb("afb", "A_{FB}", iafb, -1000., 1000.);
  		
  /////////////////////////////////////////////////////////  Acc X Reco  ///////////////////////////////////////
  //	Acceptance
  	RooRealVar accP0("accP0","accP0",readParam(iBin,"acc", 0));
  	RooRealVar accP1("accP1","accP1",readParam(iBin,"acc", 1));
  	RooRealVar accP2("accP2","accP2",readParam(iBin,"acc", 2));
  	RooRealVar accP3("accP3","accP3",readParam(iBin,"acc", 3));   
//  	accP0.setConstant(kTRUE);
//  	accP1.setConstant(kTRUE);
//  	accP2.setConstant(kTRUE);
//  	accP3.setConstant(kTRUE);
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
  	RooRealVar effP10("effP10","effP10",readParam(iBin,"accXrecoEff", 10)); 
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
  	effP10.setConstant(kTRUE);
//  	effP0.setError(readParam(iBin,"accXrecoEffErr", 0));
//  	effP1.setError(readParam(iBin,"accXrecoEffErr", 1));
//  	effP2.setError(readParam(iBin,"accXrecoEffErr", 2));
//  	effP3.setError(readParam(iBin,"accXrecoEffErr", 3));  
//  	effP4.setError(readParam(iBin,"accXrecoEffErr", 4));  
//  	effP5.setError(readParam(iBin,"accXrecoEffErr", 5)); 
//  	effP6.setError(readParam(iBin,"accXrecoEffErr", 6)); 
//  	effP7.setError(readParam(iBin,"accXrecoEffErr", 7)); 
//  	effP8.setError(readParam(iBin,"accXrecoEffErr", 8)); 
//  	effP9.setError(readParam(iBin,"accXrecoEffErr", 9)); 
  /////////////////////////////////////////////////////////  Total Efficiency  ///////////////////////////////////////
  
  ///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
  	// // Signal double gaussian
  	RooRealVar sigGauss_mean("sigGauss_mean","M_{K^{#pm}#Mu#Mu}",5.279,5.26,5.30);
  	RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0));
  	sigGauss1_sigma.setConstant(kTRUE);
  	RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0));
  	sigGauss2_sigma.setConstant(kTRUE);
  	RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0));
  	sigM_frac.setConstant(kTRUE);
  	// // mass distro of signal
  	RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
  	RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
  	RooAddPdf f_sigM("f_sigM","f_sigM", RooArgList(f_sigMGauss1, f_sigMGauss2), sigM_frac);
  	
  	RooArgSet f_sigA_argset(CosThetaL);
  	f_sigA_argset.add(RooArgSet(fh,afb));
  	TString f_sigA_format;
  	TString f_rec_format;
  	TString f_ang_format;
    //f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
  	f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
  	if (iBin != 0 && iBin != 1 && iBin != 9) {
  		  f_sigA_argset.add(RooArgSet(effP0, effP1, effP2, effP3, effP4, effP5, effP6));
  		  f_rec_format = "( effP0+effP1*CosThetaL+effP2*CosThetaL**2+effP3*CosThetaL**3+effP4*CosThetaL**4+effP5*CosThetaL**5+effP6*CosThetaL**6 )";
  	} else {
  		  f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
  		  f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
  		  f_rec_format = "( accP0 + accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
  	}
  	f_sigA_argset.add(RooArgSet(f_rec_format));
  	f_sigA_argset.add(RooArgSet(f_ang_format));
  	f_sigA_format = TString::Format("%s * %s",f_rec_format.Data(),f_ang_format.Data());
  	RooGenericPdf f_sigA("f_sigA", f_sigA_format,f_sigA_argset);
  	// Create signal distribution
  	RooProdPdf f_sig("f_sig","f_sig",f_sigM,f_sigA);
  	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_sig prepared. <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  	
  	// Create combinatorial background distribution
  	RooRealVar bkgCombM_c("bkgCombM_c","c",-0.5,-20.,1.);
  	RooRealVar offset("offset","offset",-5.);
  	RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgList(Bmass,offset));
  	RooExponential f_bkgCombM("f_bkgCombM","f_bkgCombM",Bmass_offset,bkgCombM_c);// exponential decay
  	
  	RooRealVar bkgCombL_c0("bkgCombL_c0","c0",readParam(iBin,"bkgCombL_c0",0));
  	RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0));
  	RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0));
  	RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0)); 
  	RooArgSet f_bkgCombL_argset;
  	switch (iBin) {
  		  default:
  		      f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
  				  bkgCombL_c0.setConstant(kTRUE);
  				  bkgCombL_c1.setConstant(kTRUE);
  				  bkgCombL_c2.setConstant(kTRUE);
  				  if (iBin == 8 || iBin == 7 || iBin == 2 || iBin ==9) {
  				  f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
  				  bkgCombL_c3.setConstant(kTRUE);
        }
        break;
  	}
  	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
  	RooRealVar bkgGauss_mean("bkgGauss_mean","cos#theta_{l}", readParam(iBin,"bkgGauss_mean",0));
  	bkgGauss_mean.setConstant(kTRUE);
    RooRealVar bkgGauss_sigma("bkgGauss_sigma","#sigma",  readParam(iBin,"bkgGauss_sigma",0));
  	bkgGauss_sigma.setConstant(kTRUE);
    RooRealVar bkg_frac("bkg_frac","bkg_frac",readParam(iBin,"bkg_frac",0));
  	bkg_frac.setConstant(kTRUE);
  	
  	RooGaussian f_bkgCombLGauss("f_bkgCombLGauss","f_bkgCombLGauss", CosThetaL, bkgGauss_mean, bkgGauss_sigma);
  	RooAddPdf f_bkgCombL("f_bkgCombL","f_bkgCombL", RooArgList(f_bkgCombLGauss, f_bkgCombL_P), bkg_frac);
  	
  	RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb",f_bkgCombL, f_bkgCombM);
  	cout<<">>>>>>>>>>>>>>>> INFO: f_bkgComb prepared. <<<<<<<<<<<<<<<<<<<<<<"<<endl;
  	
    RooRealVar nsig("nsig","nsig",500,0,4E3);
  	RooRealVar nbkgComb("nbkgComb","nbkgComb",1000,0,6E3);
  	
  	RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_sig),RooArgList(nbkgComb,nsig));// no penalty term
  	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
  ///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
  
  	// Get data and apply unbinned fit
  	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
//  	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)));	
//  	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//  	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//    f_fitresult->Print();

    // Fitting procedure in TMinuit
    double isMigradConverge[2] = {-1,0};
    double isMinosValid = -1;
    RooAbsReal *nll = f.createNLL(*data,Extended(kTRUE),Offset(kFALSE),NumCPU(1));// Minos and Save are unknown.
    RooMinuit minuit(*nll);
    printf("INFO\t\t: Start MIGRAD loop\n");
    for(int iLoop = 0; iLoop < 10; iLoop++){
        isMigradConverge[0] = minuit.migrad();
        printf("INFO\t\t: MIGRAD return code=%.0f\n",isMigradConverge[0]);
        if (isMigradConverge[0] == 0) break;
    }
    isMigradConverge[1] = minuit.save()->minNll();
    int gKeepParam = 1;
    odatacardpath=summarypath;
    if (gKeepParam) {
        writeParam(iBin, "migrad", isMigradConverge);
        double val[4]={0,0,0,0};
        val[0] = toBoundedFh( fh.getVal() );
//        val[1] = fh.getError();val[2]=fh.getErrorLo();val[3]=fh.getErrorHi();
        writeParam(iBin, "fh_migrad", val, 4);
        val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//        val[1] = afb.getError();val[2]=afb.getErrorLo();val[3]=afb.getErrorHi();
        writeParam(iBin, "afb_migrad",val, 4);
    }
    double isHesseValid = minuit.hesse();
    // Keep HESSE result as preliminary
    if (gKeepParam) {
        writeParam(iBin, "hesse", &isHesseValid, 1);
        minuit.save();
        double val[4]={0,0,0,0};
        val[0] = toBoundedFh( fh.getVal() );
//        val[1] = fh.getError();val[2]=fh.getErrorLo();val[3]=fh.getErrorHi();
        writeParam(iBin, "fh_hesse", val, 4);
        val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//        val[1] = afb.getError();val[2]=afb.getErrorLo();val[3]=afb.getErrorHi();
        writeParam(iBin, "afb_hesse",val, 4);
    }
//    printf("INFO\t\t: Start MINOS loop\n");
//    for(int iLoop = 0; iLoop < 3; iLoop++){
//        isMinosValid = minuit.minos(RooArgSet(afb,fh,nsig));
//        printf("INFO\t\t: MINOS return code=%.0f\n",isMinosValid);
//        if (isMinosValid == 0) break;
//    }
//    if (gKeepParam) {
//        writeParam(iBin, "minos", &isMinosValid, 1);
//    }
    minuit.save();

    delete data;
    odatacardpath=summarypath;
    double val[4]={0,0,0,0};
    val[0] = toBoundedFh( fh.getVal() );
//    val[1] = fh.getError();val[2]=fh.getErrorLo();val[3]=fh.getErrorHi();
    writeParam(iBin, "FC_fh", val, 4);
    val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//    val[1] = afb.getError();val[2]=afb.getErrorLo();val[3]=afb.getErrorHi();
    writeParam(iBin, "FC_afb",val, 4);
    val[0] = minuit.save()->minNll();
		writeTXT(iBin, "FC_FCN", val, 2);

//    val[0] = f_fitresult->status();
//    writeParam(iBin,  "migrad",  val, 1);
//    val[0] = f_fitresult->minNll();
//    writeParam(iBin, "FC_FCN", val);
////    val[0] = 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
//    val[0] = toBoundedFh( fh.getVal() );
//    writeParam(iBin,  "FC_fh",  val, 2);
////    val[0] = (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() ); 
//    val[0] = toBoundedAfb( afb.getVal(), fh.getVal() );
//    writeParam(iBin,  "FC_afb",  val, 2);

//    // map of signal
//    bkgCombM_c.setConstant(kTRUE);
//    offset.setConstant(kTRUE);
//    fh.setConstant(kTRUE);
//    afb.setConstant(kTRUE);
//    nsig.setConstant(kTRUE);
//    nbkgComb.setConstant(kTRUE);
//    sigGauss_mean.setConstant(kTRUE);
//    RooWorkspace *wspace2 = new RooWorkspace("wspace","wspace");
//    wspace2->import(f);
//    wspace2->writeToFile(TString::Format("%s/wspace_FC_bin%d.root",summarypath.Data(),iBin),true);
  
  	std::vector<double> output;
  	output.push_back(fh.getVal());
  	output.push_back(fh.getError());
  	output.push_back(afb.getVal());
  	output.push_back(afb.getError());
  	return output;
  	
}//}}}

std::vector<double> angular2D_toy_bin(int iBin, const char outfile[] = "angular2D_toy")
{//{{{
	setTDRStyle();
	cout<<endl<<"iBin = "<<iBin<<endl<<endl; 
	// Create parameters and PDFs
	RooRealVar CosThetaL("CosThetaL", "cos#theta_{l}", -1., 1.);
	RooRealVar Bmass("Bmass","M_{K^{#pm}#Mu#Mu}",5.10,5.60);
	RooRealVar Q2("Q2","q^{2}",1.0,22.);
//	 double Iafb = readParam(iBin, "afb", 0);  //  mean values!
//   double Ifh  = readParam(iBin, "fh", 0);
//	 double iafb = readParam(iBin, "Iafb_angular2D", 0);  // Selected initial values!
//   double ifh  = readParam(iBin, "Ifh_angular2D", 0);
//   double Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
//   double Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
//   cout<<Iafb<<"  "<<Ifh<<endl;
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    //TFile *f_wspace_sigM = new TFile(TString::Format("%s/wspace_Sm_bin%d.root",iwspacepath.Data(),iBin));
    //TFile *f_wspace_comb = new TFile(TString::Format("%s/wspace_comb_bin%d.root",iwspacepath.Data(),iBin));
    TFile *f_wspace_pdf  = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    //RooWorkspace *wspace_sigM = (RooWorkspace*)f_wspace_sigM->Get("wspace");
    //RooWorkspace *wspace_comb = (RooWorkspace*)f_wspace_comb->Get("wspace");
    RooWorkspace *wspace_pdf  = (RooWorkspace*)f_wspace_pdf->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    RooGenericPdf *f_sigM = 0;
    RooGenericPdf *f_bkgComb = 0;
    //RooGenericPdf *f = 0;
    RooRealVar *fh = 0;
    RooRealVar *afb = 0;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        fh = (RooRealVar*)wspace_sigA->var("fh");
        afb = (RooRealVar*)wspace_sigA->var("afb");
        fh->setRange(0,3);
        fh->setVal(readParam(iBin, "fh", 0));
        afb->setRange(-1,1);
        afb->setVal(readParam(iBin, "afb", 0));
    }
/*    if (wspace_sigM){
       f_sigM = (RooGenericPdf*)wspace_sigM->pdf("f_sigM");
    }
    if (wspace_comb){
       f_bkgComb = (RooGenericPdf*)wspace_comb->pdf("f_bkgComb");
    }
*/    
    if (wspace_pdf){
       f_sigM = (RooGenericPdf*)wspace_pdf->pdf("f_sigM");
       f_bkgComb = (RooGenericPdf*)wspace_pdf->pdf("f_bkgComb");
       //f = (RooGenericPdf*)wspace_pdf->pdf("kernel");
    }
    fh->Print();
    afb->Print();
	RooRealVar nsig("nsig","nsig",50,0,4E3);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",100,0,6E3);
	RooProdPdf *f_sig=new RooProdPdf("f_sig","f_sig",*f_sigM,*f_sigA);
   RooAddPdf *f=0;
   f = new RooAddPdf("kernel","kernel",RooArgList(*f_bkgComb,*f_sig),RooArgList(nbkgComb,nsig));// no penalty term
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	

	// Get data and apply unbinned fit
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
//	ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
	
	RooFitResult *f_fitresult = f->fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3) {
	std::vector<double> output;
	//output.push_back(fh->getVal());
	//output.push_back(fh->getError());
	output.push_back(afb->getVal());
	output.push_back(afb->getError());
	return output;
	}

	std::vector<double> output;
	output.push_back(afb->getVal());
	return output;
	
}//}}}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void genToySignal(int iBin, int nEvents = 100, double newAfb = -1., double newFh = -1.) // For validation.
{//{{{
    //int defNEvents  = 174;
    //double defAfb   = readParam(iBin,"genafb",0);
    //double defFh    = readParam(iBin,"genfh",0);
    double defAfb   = -1;
    double defFh    = -1;
    double defAfbErr= fabs(readParam(iBin,"genafb",1));
    double defFhErr = fabs(readParam(iBin,"genfh",1));

    if (newAfb != -1.) defAfb = newAfb;
    if (newFh  != -1.) defFh  = newFh;

    // Get PDF and parameters
    TFile *f_wspace_Sm = new TFile(TString::Format("%s/wspace_Sm_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_Sm = (RooWorkspace*)f_wspace_Sm->Get("wspace");
    RooAddPdf *f_sigM = 0;
    RooRealVar *Bmass = 0;
    if (wspace_Sm){
        f_sigM = (RooAddPdf*)wspace_Sm->pdf("f_sigM");
        Bmass = (RooRealVar*)wspace_Sm->var("Bmass");
    }else{
        printf("ERROR\t\t: Please have wsapce_Sm_bin?.root prepared.\n");
        return;
    }
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    RooRealVar *CosThetaL = 0;
    RooRealVar *afb;
    RooRealVar *fh;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        CosThetaL = (RooRealVar*)wspace_sigA->var("CosThetaL");
        afb       = (RooRealVar*)wspace_sigA->var("afb");
        fh        = (RooRealVar*)wspace_sigA->var("fh");
    }else{
        printf("ERROR\t\t: Please have wsapce_sigA_bin?.root prepared.\n");
        return;
    }
    // Random generator for afb/fh values
    RooGaussian gaus_afb("gaus_afb","",*afb,RooConst(defAfb),RooConst(defAfbErr));
    RooGaussian gaus_fh ("gaus_fh" ,"" ,*fh ,RooConst(defFh) ,RooConst(2.*defFhErr));
    gaus_afb.Print();
    gaus_fh.Print();
    //if (nEvents == 0) nEvents = defNEvents;
    RooProdPdf *f = new RooProdPdf("f", "f", RooArgSet(*f_sigM,*f_sigA));

    double rndBmass(0);
    double rndCosL(0);
    double Q2((Q2rangeup[iBin]+Q2rangedn[iBin])/2);
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    int    Triggers(1);

    TFile fout(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");
    TH1D *h_Afb     =new TH1D("h_Afb","",100,-1,1);
    TH1D *h_Fh      =new TH1D("h_Fh","",100,0,3);
    TH2D *h2_AfbFh  =new TH2D("h2_AfbFh","",100,-1,1,100,0,3);
    //TH1D *h_ubdAfb     =new TH1D("h_ubdAfb","",1000,-100,100);
    //TH1D *h_ubdFh      =new TH1D("h_ubdFh","",1000,-100,100);
    //TH2D *h2_ubdAfbFh  =new TH2D("h2_ubdAfbFh","",1000,-100,100,1000,-100,100);

    printf("INFO\t\t: Enter generating loop\n");
    bool needAfbFhError = false; // VERY SLOW if this is true
    if (needAfbFhError){
        const RooArgSet *dataInEntry = 0;
        RooDataSet* data_fh = gaus_fh .generate(*fh ,nEvents);
        RooDataSet* data_afb = gaus_afb.generate(*afb,nEvents);
        for (int iEvt = 0; iEvt < nEvents; iEvt++) {
            double temp_fh  = fabs( ((RooRealVar*)(data_fh->get(iEvt)->find("fh")))->getVal() );
            double temp_afb = ((RooRealVar*)(data_afb->get(iEvt)->find("afb")))->getVal();
            h_Fh        ->Fill(temp_fh);
            h_Afb       ->Fill(temp_afb);
            h2_AfbFh    ->Fill(temp_afb,temp_fh);
            fh ->setVal(toUnboundedFh(temp_fh));
            afb->setVal(toUnboundedAfb(temp_afb,temp_fh));
            //h_ubdFh     ->Fill(fh->getVal());
            //h_ubdAfb    ->Fill(afb->getVal());
            //h2_ubdAfbFh ->Fill(afb->getVal(),fh->getVal());
            dataInEntry = f->generate(RooArgSet(*Bmass,*CosThetaL),1)->get(0);
            rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
            rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
            tree->Fill();
        }
    }else{
        //fh ->setVal(toUnboundedFh(defFh));
        //afb->setVal(toUnboundedAfb(defAfb,defFh));
        fh ->setVal(defFh);
        afb->setVal(defAfb);
        const RooDataSet *dataInEntry = 0;
        dataInEntry = f->generate(RooArgSet(*Bmass,*CosThetaL),nEvents);
        for (int iEvt = 0; iEvt < nEvents; iEvt++) {
            //h_Fh        ->Fill(defFh);
            //h_Afb       ->Fill(defAfb);
            //h2_AfbFh    ->Fill(defAfb,defFh);
            h_Fh     ->Fill(fh->getVal());
            h_Afb    ->Fill(afb->getVal());
            h2_AfbFh ->Fill(afb->getVal(),fh->getVal());
            rndBmass=((RooRealVar*)(dataInEntry->get(iEvt))->find("Bmass"    ))->getVal();
            rndCosL =((RooRealVar*)(dataInEntry->get(iEvt))->find("CosThetaL"))->getVal();
            tree->Fill();
        }
    }
    fout.Write();
    fout.Close();

    f_wspace_Sm->Close();
    f_wspace_sigA->Close();
    delete f_wspace_Sm;
    delete f_wspace_sigA;
    return;
}//}}}

void genToyCombBkg(int iBin, int nEvents = 0, const char outfile[] = "combBkgToy") // For validation.
{//{{{
    TFile *f_wspace_comb_M = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_comb_M = (RooWorkspace*)f_wspace_comb_M->Get("wspace");
    RooRealVar *Bmass = 0;
    RooGenericPdf *f_bkgCombM = 0;
    RooRealVar *CosThetaL = 0;
    RooGenericPdf *f_bkgCombA = 0;
    if (wspace_comb_M){
        Bmass = (RooRealVar*)wspace_comb_M->var("Bmass");
        f_bkgCombM = (RooGenericPdf*)wspace_comb_M->pdf("f_bkgCombM");
        f_bkgCombA = (RooGenericPdf*)wspace_comb_M->pdf("f_bkgCombL");
        CosThetaL = (RooRealVar*)wspace_comb_M->var("CosThetaL");
    }
    RooProdPdf f("f","f",RooArgSet(*f_bkgCombM,f_bkgCombA));
    
    double rndBmass(0);
    double rndCosL(0);
    double Q2((Q2rangeup[iBin]+Q2rangedn[iBin])/2);
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    int    Triggers(1);
    
    TFile fout(TString::Format("%s/%s_Bin%d.root",owspacepath.Data(),outfile,iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");

    cout<<" >>>>>>>>>>>>>>>>>>>>>> DONE 1 "<<nEvents<<endl<<endl;
    RooDataSet *outdata = f.generate(RooArgSet(*Bmass,*CosThetaL),nEvents);
    cout<<" >>>>>>>>>>>>>>>>>>>>>> DONE 2 "<<nEvents<<endl<<endl;
    const RooArgSet *dataInEntry = 0;
    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
        dataInEntry = outdata->get(iEvt);
        rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
        rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
        tree->Fill();
    }
    fout.Write();
    fout.Close();

    f_wspace_comb_M->Close();
    delete f_wspace_comb_M;
}//}}}

void splitMCSamples(double oLumis=datasetLumi[0]) // Split MC samples into small parts
{//{{{
    static char decmode[20];
    while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s")*strcmp(decmode, "signal")*strcmp(decmode, "manual") != 0){
        printf("Please insert type [ signal / jpsi / psi2s / manual ]:");
        scanf("%19s",decmode);
    }
    
    double iLumis = 0.;
    if (strcmp(decmode, "signal") == 0){
        iLumis=datasetLumi[1];
    }else if (strcmp(decmode, "jpsi") == 0){
        iLumis = datasetLumi[2];
    }else if (strcmp(decmode, "psi2s") == 0){
        iLumis = datasetLumi[3];
    }else{
        printf("Input luminosity: ");
        scanf("%lf",&iLumis);
    }
    printf("INFO\t\t: Splitting %f/fb into %.0f parts\n",iLumis, floor(iLumis/oLumis));
    
    if (oLumis > iLumis){
        printf("ERROR\t\t: Luminosity of output should be larger than the reservoir. EXIT 1.");
        return;
    }

    int mumuMassWindowBin = 0;
    for(int iPart=0; iPart < floor(iLumis/oLumis); iPart++){
        TFile *fout = new TFile(TString::Format("%s/splitMC_%s_part%04d.root",owspacepath.Data(),decmode,iPart+1), "RECREATE");
        TTree *tout = ch->CopyTree("","",(int)floor(oLumis/iLumis*ch->GetEntries()),iPart*(int)floor(oLumis/iLumis*ch->GetEntries()))->CopyTree(mumuMassWindow[mumuMassWindowBin]);
        fout->Write();
        fout->Close();
        tout->Delete();
    }
    return;
}//}}} 

// void splitToySamples(const char iFile[], const int nEvtsPerSet, const char oFileLabel[]="splitToy" )
// {//{{{
//     TFile *fin = new TFile(iFile);
//     TTree *tree = (TTree*)fin->Get("tree");
//     for(int iPart=0; iPart < ceil(tree->GetEntries())/nEvtsPerSet; iPart++){
//         TFile *fout = new TFile(TString::Format("%s/%s_set%04d.root",owspacepath.Data(),oFileLabel,iPart+1), "RECREATE");
//         TTree *tout = tout->CopyTree("","",nEvtsPerSet,iPart*nEvtsPerSet);
//         fout->Write();
//         fout->Close();
//     }
//     fin->Close();
//     return;
// }//}}}

void splitToySamples(const char iFile[], const int nSets, const int *nEvtsPerSet, const char oFileLabel[]="splitToy" )
{//{{{
    TFile *fin = new TFile(iFile);
    TTree *tree = (TTree*)fin->Get("tree");
    long int evtCounter = 0;
    for(int iPart=0; iPart < nSets; iPart++){
        TFile *fout = new TFile(TString::Format("%s/%s_set%04d.root",owspacepath.Data(),oFileLabel,iPart+1), "RECREATE");
        TTree *tout = tree->CopyTree("","",nEvtsPerSet[iPart],evtCounter);
        evtCounter+=nEvtsPerSet[iPart];
        fout->Write();
        fout->Close();
        //tout->Delete();
    }
    fin->Close();
    return;
}//}}}


void createFCToys(int iBin, int nToy=500)// Create toys for contour scanning
{//{{{
    // Remark: keep extension possibility for 2-D contour.
    // Buffers for checking directory/FILE
    //bool contourMode=false;// In case contour is needed.
    struct stat fiBuff;
    //FILE *fi = 0;

    // Setting
    TString otoyspath = TString::Format("./ToyMC/bin%d",iBin);
    const double stepSizeAfb[11] = { 0.0200,  0.0300, 0.0060, 0.00, 0.0060, 0.00, 0.0060, 0.0060, 0.0080, 0.0200, 0.0040};
    const double stepSizeFh[11]  = { 0.0200,  0.0200, 0.0060, 0.00, 0.0060, 0.00, 0.0080, 0.0060, 0.0100, 0.0200, 0.0040};
    const double dataAfbLo[11]  = {-0.500,-0.600,-0.150, 0.00,-0.150, 0.00,-0.200,-0.200,-0.250,-0.400,-0.100};
    const double dataFhBand[11] = { 1.40,   1.50,  0.30, 0.00,  0.30, 0.00,  0.40,  0.50,  0.60,  1.00,  0.20};

    // Create output directory
    if (stat(TString::Format("./ToyMC"),&fiBuff) != 0){
        mkdir(TString::Format("./ToyMC"),0755);
        if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
            mkdir(TString::Format("%s",otoyspath.Data()),0755);
        }
    }

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
//    double  fh   = wspace->var("fh")->getVal();
//    double  afb  = wspace->var("afb")->getVal();
    double  nsig = wspace->var("nsig")->getVal();
    double  nbkgComb = wspace->var("nbkgComb")->getVal();

    TRandom3 *rndGenerator = new TRandom3();
    int nsigInToy[nToy];
    int nbkgCombInToy[nToy];
    int nsigInToys = 0;
    int nbkgCombInToys = 0;
    for(int iToy=0; iToy<nToy; iToy++){
        nsigInToy[iToy]=rndGenerator->Poisson(nsig);
        nbkgCombInToy[iToy]=rndGenerator->Poisson(nbkgComb);
        nsigInToys+=nsigInToy[iToy];
        nbkgCombInToys+=nbkgCombInToy[iToy];
    }

    const double dnA1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    double par_Afb[4]; 
    par_Afb[0] = readParam(iBin,"Pro_Afb_shape", 0);
    par_Afb[1] = readParam(iBin,"Pro_Afb_shape", 1);
    par_Afb[2] = readParam(iBin,"Pro_Afb_shape", 2);
    par_Afb[3] = readParam(iBin,"Pro_Afb_shape", 3);
    // Loop over phase space
    double thisAfb = 0.;
    double thisFh =0.01;
    for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
        thisFh = iFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
        }
        owspacepath=TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/set0001",owspacepath.Data()),&fiBuff) == 0) { cout<<thisAfb<<"  "<<thisFh<<endl; continue;}

        genToySignal(iBin, nsigInToys,thisAfb,thisFh);
        splitToySamples(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
        genToyCombBkg(iBin,nbkgCombInToys);
        cout<<" >>>>>>>>>>>>>>>>>>>>>> DONE "<<nbkgCombInToys<<endl<<endl;
        splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
        for(int iToy = 0; iToy<nToy; iToy++){
            if (stat(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),&fiBuff) != 0){
                mkdir(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),0755);
            }
            rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
            rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
        }
    }// Afb loop

    const double dnF1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double dnF2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    double par_Fh[7]; 
    par_Fh[0] = readParam(iBin,"Pro_Fh_shape", 0);
    par_Fh[1] = readParam(iBin,"Pro_Fh_shape", 1);
    par_Fh[2] = readParam(iBin,"Pro_Fh_shape", 2);
    par_Fh[3] = readParam(iBin,"Pro_Fh_shape", 3);
    par_Fh[4] = readParam(iBin,"Pro_Fh_shape", 4);
    par_Fh[5] = readParam(iBin,"Pro_Fh_shape", 5);
    par_Fh[6] = readParam(iBin,"Pro_Fh_shape", 6);
    thisAfb = 0.;
    thisFh =0.01;
    cout<<"  ---    1  ---"<<endl;
    for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        cout<<thisAfb<<"  "<<thisFh<<endl;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
        }
        owspacepath=TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/set0001",owspacepath.Data()),&fiBuff) == 0) { cout<<thisAfb<<"  "<<thisFh<<endl; continue;}

        genToySignal(iBin, nsigInToys,thisAfb,thisFh);
        splitToySamples(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
        genToyCombBkg(iBin,nbkgCombInToys);
        splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
        for(int iToy = 0; iToy<nToy; iToy++){
            if (stat(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),&fiBuff) != 0){
                mkdir(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),0755);
            }
            rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
            rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
        }
    }// Fh loop
    cout<<"  ---    2  ---"<<endl;
    return;
}//}}}

void harvestFCFitResults(int iBin, int nToy=500, bool contourMode=false)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    contourMode=false;// In case contour is needed.
    // Setting
    TString otoyspath = TString::Format("./ToyMC/bin%d",iBin);
    const int    nHistBins   = 1000;
//    const int    nHistBins   = 100;
    const double stepSizeAfb[11] = { 0.0200,  0.0300, 0.0060, 0.00, 0.0060, 0.00, 0.0060, 0.0060, 0.0080, 0.0200, 0.0040};
    const double stepSizeFh[11]  = { 0.0200,  0.0200, 0.0060, 0.00, 0.0060, 0.00, 0.0080, 0.0060, 0.0100, 0.0200, 0.0040};
    const double dataAfbLo[11]  = {-0.500,-0.600,-0.150, 0.00,-0.150, 0.00,-0.200,-0.200,-0.250,-0.400,-0.100};
    const double dataFhBand[11] = { 1.40,   1.50,  0.30, 0.00,  0.30, 0.00,  0.40,  0.50,  0.60,  1.00,  0.20};

    // Get parameters for the q2 Bin
    //iwspacepath="./RootFiles";
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
    printf("INFO\t\t: bounded fh=%+.3f, bounded afb=%+.3f\n",fh,afb);

    const double dnA1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    double par_Afb[4]; 
    par_Afb[0] = readParam(iBin,"Pro_Afb_shape", 0);
    par_Afb[1] = readParam(iBin,"Pro_Afb_shape", 1);
    par_Afb[2] = readParam(iBin,"Pro_Afb_shape", 2);
    par_Afb[3] = readParam(iBin,"Pro_Afb_shape", 3);
    // Loop over phase space
    double thisAfb = 0.;
    double thisFh =0.01;
    for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
        thisFh = iFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        owspacepath=TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;

        TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
        TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
        TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);

        for(int iToy = 0; iToy<nToy; iToy++){
            //iwspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            idatacardpath=TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            if (readParam(iBin,"migrad", 0) == 0.){
                if (readParam(iBin,"FC_fh",0) > 3e-3) {  // bin4 
//                if (readParam(iBin,"fh_migrad",0) > 3e-3) { 
                        //if (readParam(iBin,"FC_afb",0) < -0.6) continue;
                        //cout<<idatacardpath<<"  "<<readParam(iBin,"FC_afb",0)<<endl;
                    h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));  // bin4
                    h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0)); // bin4
//                    h_setSummaryFh  ->Fill(readParam(iBin,"fh_migrad",0));
//                    h_setSummaryAfb ->Fill(readParam(iBin,"afb_migrad",0));
                }
            }
        }
        h_setSummaryAfb->Draw();
        h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb.cc",owspacepath.Data()));
        h_setSummaryFh->Draw();
        h_setSummaryFh->SaveAs(TString::Format("%s/h_setSummaryFh.cc",owspacepath.Data()));
        fout->Write();
        fout->Close();
    }// Afb loop
    const double dnF1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double dnF2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    double par_Fh[7]; 
    par_Fh[0] = readParam(iBin,"Pro_Fh_shape", 0);
    par_Fh[1] = readParam(iBin,"Pro_Fh_shape", 1);
    par_Fh[2] = readParam(iBin,"Pro_Fh_shape", 2);
    par_Fh[3] = readParam(iBin,"Pro_Fh_shape", 3);
    par_Fh[4] = readParam(iBin,"Pro_Fh_shape", 4);
    par_Fh[5] = readParam(iBin,"Pro_Fh_shape", 5);
    par_Fh[6] = readParam(iBin,"Pro_Fh_shape", 6);
    thisAfb = 0.;
    thisFh =0.01;
    cout<<"  ---    1  ---"<<endl;
    for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
//        cout<<thisAfb<<"  "<<thisFh<<endl;
        if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        owspacepath=TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;

        TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
        TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
        TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);

        for(int iToy = 0; iToy<nToy; iToy++){
            idatacardpath=TString::Format("%s/FH_afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            if (readParam(iBin,"migrad", 0) == 0.){
                if (readParam(iBin,"FC_fh",0) > 3e-3) {  // bin4 
//                if (readParam(iBin,"fh_migrad",0) > 3e-3) { 
                    h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));  // bin4
                    h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0)); // bin4
//                    h_setSummaryFh  ->Fill(readParam(iBin,"fh_migrad",0));
//                    h_setSummaryAfb ->Fill(readParam(iBin,"afb_migrad",0));
                }
            }
        }
        h_setSummaryAfb->Draw();
        h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb.cc",owspacepath.Data()));
        h_setSummaryFh->Draw();
        h_setSummaryFh->SaveAs(TString::Format("%s/h_setSummaryFh.cc",owspacepath.Data()));
        fout->Write();
        fout->Close();
    }// Fh loop
    
    return;

}//}}}
void getCIFromTH1F(TH1F *hin, double &lowerBd, double &upperBd, double coverage=0.684)
{//{{{
    unsigned int nEntries = 0;
    for(int Ibin=1; Ibin <= hin->GetNbinsX(); Ibin++){
       nEntries+=hin->GetBinContent(Ibin);
    }

    unsigned int nMinIntervals = 0;
    float lowerBin = 1;
    float upperBin = hin->GetNbinsX();

    unsigned int minBinInterval = hin->GetNbinsX();
    unsigned int buffEntries = 0;
    for(int lBin=1; lBin <= hin->GetNbinsX(); lBin++){
        buffEntries = 0;
        for(int Ibin=lBin; Ibin <= hin->GetNbinsX(); Ibin++){
            buffEntries += hin->GetBinContent(Ibin);
            if (buffEntries > coverage*nEntries){
                if (Ibin-lBin < minBinInterval){
                    minBinInterval = Ibin-lBin;
                    lowerBin = lBin;
                    upperBin = Ibin;
                    nMinIntervals = 1;
                }else if (Ibin-lBin == minBinInterval){
                    lowerBin = lowerBin*nMinIntervals+lBin;
                    upperBin = upperBin*nMinIntervals+Ibin;
                    nMinIntervals++;
                    lowerBin /= nMinIntervals;
                    upperBin /= nMinIntervals;
                }
                break;
            }
        }
    }
    upperBd = hin->GetBinCenter(floor(upperBin));
    lowerBd = hin->GetBinCenter(floor(lowerBin));
    return;
}//}}}
void getFCInterval(int iBin, int nToy=500)
{//{{{
	  setTDRStyle();
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Settings
    TString otoyspath = TString::Format("./ToyMC/bin%d",iBin);
    const double stepSizeAfb[11] = { 0.0200,  0.0300, 0.0060, 0.00, 0.0060, 0.00, 0.0060, 0.0060, 0.0080, 0.0200, 0.0040};
    const double stepSizeFh[11]  = { 0.0200,  0.0200, 0.0060, 0.00, 0.0060, 0.00, 0.0080, 0.0060, 0.0100, 0.0200, 0.0040};
    const double dataAfbLo[11]  = {-0.500,-0.600,-0.150, 0.00,-0.150, 0.00,-0.200,-0.200,-0.250,-0.400,-0.100};
    const double dataFhBand[11] = { 1.40,   1.50,  0.30, 0.00,  0.30, 0.00,  0.40,  0.50,  0.60,  1.00,  0.20};
    double dataAfbLower[11], dataFhLower[11];

    cout<<"  ---    2  ---"<<endl;
    // Get parameters for the q2 Bin
    //iwspacepath="./RootFiles";
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
    dataAfbLower[iBin] = -0.5*fh;
    dataFhLower[iBin]  = fabs(afb)*2.0;
    // output
    int maxPoints = 1024;
    int nFhPoints = 0;
    int nAfbPoints = 0;
    double fhMeasuredHi[maxPoints];
    double fhMeasuredLo[maxPoints];
    double fhTruth[maxPoints];
    double afbMeasuredHi[maxPoints];
    double afbMeasuredLo[maxPoints];
    double afbTruth[maxPoints];
    double outFCErrFh[2];
    double outFCErrAfb[2];
    memset(fhMeasuredHi,-2,maxPoints*sizeof(double));
    memset(fhMeasuredLo,-2,maxPoints*sizeof(double));
    memset(fhTruth,-2,maxPoints*sizeof(double));
    memset(afbMeasuredHi,-2,maxPoints*sizeof(double));
    memset(afbMeasuredLo,-2,maxPoints*sizeof(double));
    memset(afbTruth,-2,maxPoints*sizeof(double));

    // Loop over phase space
    double thisAfb = 0.;
    double thisFh =0.01;
    TFile *fin = 0;
    TFile *fin1 = 0;
    TH1F  *h_setSummaryAfb = 0;
    TH1F  *h_setSummaryFh  = 0;
    // loop for errAfb
    cout<<"0"<<endl;
    const double dnF1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double dnF2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    double par_Fh[7]; 
    par_Fh[0] = readParam(iBin,"Pro_Fh_shape", 0);
    par_Fh[1] = readParam(iBin,"Pro_Fh_shape", 1);
    par_Fh[2] = readParam(iBin,"Pro_Fh_shape", 2);
    par_Fh[3] = readParam(iBin,"Pro_Fh_shape", 3);
    par_Fh[4] = readParam(iBin,"Pro_Fh_shape", 4);
    par_Fh[5] = readParam(iBin,"Pro_Fh_shape", 5);
    par_Fh[6] = readParam(iBin,"Pro_Fh_shape", 6);
    for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        fin = new TFile(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
        h_setSummaryAfb = (TH1F*)fin->Get("h_setSummaryAfb");
        //cout<<"1"<<TString::Format("%s/afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data()<<endl;
        getCIFromTH1F(h_setSummaryAfb,afbMeasuredLo[nAfbPoints],afbMeasuredHi[nAfbPoints]);
        //cout<<"3"<<endl;
        afbTruth[nAfbPoints] = thisAfb;
        nAfbPoints++;
        fin->Close();
    }
    
    thisAfb = 0.;
    thisFh =0.01;
    const double dnA1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    double par_Afb[4]; 
    par_Afb[0] = readParam(iBin,"Pro_Afb_shape", 0);
    par_Afb[1] = readParam(iBin,"Pro_Afb_shape", 1);
    par_Afb[2] = readParam(iBin,"Pro_Afb_shape", 2);
    par_Afb[3] = readParam(iBin,"Pro_Afb_shape", 3);
    for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
        thisFh = iFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        fin1 = new TFile(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
        h_setSummaryFh  = (TH1F*)fin1->Get("h_setSummaryFh");
        getCIFromTH1F(h_setSummaryFh,fhMeasuredLo[nFhPoints],fhMeasuredHi[nFhPoints]);
        fhTruth[nFhPoints] = thisFh;
        nFhPoints++;
        fin1->Close();
    }

    // Make plots, fit, and find confidence interval
    TFile *fout = new TFile(TString::Format("%s/wspace_FCConfInterval_bin%d.root",owspacepath.Data(),iBin),"RECREATE");
    TCanvas *canvas = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    TLine *line = new TLine();
    TGraph *g_fhIntervalHi = new TGraph(nFhPoints, fhMeasuredHi,fhTruth);
    TGraph *g_fhIntervalLo = new TGraph(nFhPoints, fhMeasuredLo,fhTruth);
    TGraph *g_afbIntervalHi = new TGraph(nAfbPoints, afbMeasuredHi,afbTruth);
    TGraph *g_afbIntervalLo = new TGraph(nAfbPoints, afbMeasuredLo,afbTruth);
    TF1 *f1_fhIntervalHi  = new TF1("f1_fhIntervalHi",  "[0]+[1]*x+[2]*x**2",0,3);
    TF1 *f1_fhIntervalLo  = new TF1("f1_fhIntervalLo",  "[0]+[1]*x+[2]*x**2",0,3);
    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x+[2]*x**2",-1,1);
    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x+[2]*x**2",-1,1);
    double Hi_Fh[11]  = { 0.00, 0.10, 0.00, 3.00, 0.10, 5.00, 0.00, 0.00, 0.00, 0.00, 0.00};
    double Lo_Fh[11]  = { 0.10, 0.20, 0.00, 3.00, 0.00, 5.00, 0.00, 0.00, 0.00, 0.00,-0.01};
    double Hi_Afb[11] = {-0.10, 0.00, 0.00, 3.00,-0.00, 5.00, 0.00, 0.00, 0.00, 0.00, 0.00};
    double Lo_Afb[11] = { 0.00, 0.00,-0.10, 3.00, 0.00, 5.00, 0.00, 0.00, 0.00, 0.00, 0.10};
    TFitResultPtr r_fhIntervalHi = g_fhIntervalHi->Fit("f1_fhIntervalHi","S","",fhMeasuredHi[0]+Hi_Fh[iBin]*abs(fhMeasuredHi[0]),fhMeasuredHi[nFhPoints-1]);
    TFitResultPtr r_fhIntervalLo = g_fhIntervalLo->Fit("f1_fhIntervalLo","S","",fhMeasuredLo[0]+Lo_Fh[iBin]*abs(fhMeasuredHi[0]),fhMeasuredLo[nFhPoints-1]);       
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",afbMeasuredHi[0]+Hi_Afb[iBin]*abs(afbMeasuredHi[0]),afbMeasuredHi[nAfbPoints-1]); 
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",afbMeasuredLo[0]+Lo_Afb[iBin]*abs(afbMeasuredLo[0]),afbMeasuredLo[nAfbPoints-1]);
    cout<<afbMeasuredHi[0]<<"  "<<afbMeasuredLo[0]<<endl;
    // Calculate the error, you must think about the case in which f1 not defined!
    double XPFh[2], XPAfb[2];
    XPAfb[0]=f1_afbIntervalHi->Eval(afb);
    XPAfb[1]=f1_afbIntervalLo->Eval(afb);
    if (f1_fhIntervalHi->Eval(fh) > 0.) XPFh[0]=f1_fhIntervalHi->Eval(fh);
    else XPFh[0]=0.;
    XPFh[1]=f1_fhIntervalLo->Eval(fh);
//    if (fabs(f1_afbIntervalHi->Eval(afb)) <= 0.5*fh) XPAfb[0]=f1_afbIntervalHi->Eval(afb);
//    else XPAfb[0]=-0.5*fh;
//    if (fabs(f1_afbIntervalLo->Eval(afb)) <= 0.5*fh) XPAfb[1]=f1_afbIntervalLo->Eval(afb);
//    else XPAfb[1]= 0.5*fh;
//    if (f1_fhIntervalHi->Eval(fh) > dataFhLower[iBin] && f1_fhIntervalHi->Eval(fh) < fh) XPFh[0]=f1_fhIntervalHi->Eval(fh);
//    else XPFh[0]=dataFhLower[iBin];
//    if (f1_fhIntervalLo->Eval(fh) > dataFhLower[iBin] && f1_fhIntervalLo->Eval(fh) < 3.) XPFh[1]=f1_fhIntervalLo->Eval(fh);
//    else XPFh[1]=3.;
    cout<<" FH_L = "<<XPFh[0]<<"   FH_H = "<<XPFh[1]<<endl;
    cout<<"AFB_L = "<<XPAfb[0]<<"  AFB_H = "<<XPAfb[1]<<endl;
    outFCErrFh[0] = r_fhIntervalHi.Get()  != 0 ? XPFh[0]-fh : 0.-fh;
    outFCErrFh[1] = r_fhIntervalLo.Get()  != 0 ? XPFh[1]-fh : 3.-fh;
    outFCErrAfb[0]= r_afbIntervalHi.Get() != 0 ? XPAfb[0]-afb : -0.5*fh-afb;
    outFCErrAfb[1]= r_afbIntervalLo.Get() != 0 ? XPAfb[1]-afb : 0.5*fh-afb;
    if (iBin == 2 || iBin == 4 ) {
        outFCErrAfb[0] = -0.5*fh-afb;
        outFCErrAfb[1] = 0.5*fh-afb;
    }

    const double upAfb[11] = { 0.60,  0.80, 0.10, 0.00, 0.12, 0.00, 0.10, 0.15, 0.10, 0.30, 0.10};
    const double upFh[11]  = { 1.40,  1.50, 0.40, 0.00, 0.50, 0.00, 0.40, 0.50, 0.50, 0.90, 0.40};
    // Fh plots
    gStyle->SetOptFit(0);
    g_fhIntervalLo->SetTitle("");
    g_fhIntervalLo->GetXaxis()->SetTitle("Measured F_{H}");
    g_fhIntervalLo->GetXaxis()->SetLimits(0.,upFh[iBin]);
    g_fhIntervalLo->GetYaxis()->SetTitle("True F_{H}");
//    g_fhIntervalLo->GetYaxis()->SetRangeUser(0.,1.2*fhTruth[nFhPoints-1]);
    g_fhIntervalLo->GetYaxis()->SetRangeUser(0., dataFhBand[iBin]);
    g_fhIntervalLo->Draw("AP");
    g_fhIntervalHi->Draw("P SAME");
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->SetLineColor(2);
//    line->DrawLine(fh, 0.,fh,fhTruth[nFhPoints-1]*1.1);
    line->DrawLine(fh, 0.,fh,dataFhBand[iBin]);
    //line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(0,fabs(afb*2.),upFh[iBin],fabs(afb*2.));
    latex->DrawLatexNDC(0.15,0.81,TString::Format("F_{H}=%.4f^{%+.4f}_{%+.4f}",fh,outFCErrFh[1],outFCErrFh[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->SetTextFont(12);
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_fh_bin%d.pdf",plotpath.Data(),iBin));
    // Afb plots
    g_afbIntervalHi->SetTitle("");
    g_afbIntervalHi->GetXaxis()->SetTitle("Measured A_{FB}");
    g_afbIntervalHi->GetXaxis()->SetLimits(-upAfb[iBin],upAfb[iBin]);
    g_afbIntervalHi->GetYaxis()->SetTitle("True A_{FB}");
//    g_afbIntervalHi->GetYaxis()->SetRangeUser(-1.3*afbTruth[nAfbPoints-1],1.3*afbTruth[nAfbPoints-1]);
    g_afbIntervalHi->GetYaxis()->SetRangeUser(dataAfbLo[iBin],-dataAfbLo[iBin]);
    g_afbIntervalHi->Draw("AP");
    g_afbIntervalLo->Draw("P SAME");
    line->SetLineStyle(2);
    line->SetLineColor(2);
//    line->DrawLine(afb,-1*afbTruth[nAfbPoints-1]*1.2,afb,afbTruth[nAfbPoints-1]*1.2);
    line->DrawLine(afb,dataAfbLo[iBin],afb,-dataAfbLo[iBin]);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(-upAfb[iBin], dataAfbLower[iBin],upAfb[iBin], dataAfbLower[iBin]);
    line->DrawLine(-upAfb[iBin],-dataAfbLower[iBin],upAfb[iBin],-dataAfbLower[iBin]);
    latex->DrawLatexNDC(0.15,0.81,TString::Format("A_{FB}=%.4f^{%+.4f}_{%+.4f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_afb_bin%d.pdf",plotpath.Data(),iBin));

    fout->Close();
    cout<<"  ---    3  ---"<<endl;

//    writeParam(iBin,"FCErrFh" ,outFCErrFh );
//    writeParam(iBin,"FCErrAfb",outFCErrAfb);

    return;

}//}}}
void getFCInterval2(int iBin, int nToy=500)
{//{{{
	  setTDRStyle();
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Settings
    TString otoyspath = TString::Format("./ToyMC/bin%d",iBin);
    const double targetCoverage = 0.6827;
//    double targetCoverage = 0.683;
//    double targetCoverage = 0.85;
    const double stepSizeAfb[11] = { 0.0200,  0.0300, 0.0060, 0.00, 0.0060, 0.00, 0.0060, 0.0060, 0.0080, 0.0200, 0.0040};
    const double stepSizeFh[11]  = { 0.0200,  0.0200, 0.0060, 0.00, 0.0060, 0.00, 0.0080, 0.0060, 0.0100, 0.0200, 0.0040};
//    const double dataAfbLo[11]  = {-0.500,-0.600,-0.150, 0.00,-0.150, 0.00,-0.200,-0.200,-0.250,-0.400,-0.100};
//    const double dataFhBand[11] = { 1.40,   1.50,  0.50, 0.00,  0.30, 0.00,  0.60,  0.50,  0.60,  1.00,  0.20};
    const double dataAfbLo[11]  = {-1.000,-1.000,-1.000, 0.00,-1.000, 0.00,-1.000,-1.000,-1.000,-1.000,-1.000};
    const double dataFhBand[11] = { 3.00,   3.00,  3.00, 0.00,  3.00, 0.00,  3.00,  3.00,  3.00,  3.00,  3.00};
    double dataAfbLower[11], dataFhLower[11];

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
    f_wspace->Close();
    dataAfbLower[iBin] = -0.5*fh;
    dataFhLower[iBin]  = fabs(afb)*2.0;

    // output
    const int maxPoints = 1000;
    int nFhPoints = 0;
    int nAfbPoints = 0;
    double fhTruth[maxPoints];
    double afbTruth[maxPoints];
    double outFCErrFh[2];
    double outFCErrAfb[2];
    double fhTrueValues[maxPoints];// True values for each toy set
    double afbTrueValues[maxPoints];
    TF1    *fhPDFs[maxPoints];      // P(measured|true);
    TF1    *afbPDFs[maxPoints];
    for (int i = 0; i < maxPoints; i++) {
        fhTruth[i] = -2.;
        afbTruth[i] = -2.;
        fhTrueValues[i] = -2.;
        afbTrueValues[i] = -2.;
    }
//    TLatex *lala = new TLatex();
//    lala->SetNDC();

// Loop over phase space
    TLatex *latex2 = new TLatex();
    latex2->SetNDC();
    latex2->SetTextColor(2);
    latex2->SetTextSize(0.04);
    double thisAfb = 0.;
    double thisFh =0.01;
    TFile *finAfb[maxPoints];
    TH1F  *h_setSummaryAfb[maxPoints];
    // loop for errAfb
    cout<<">>>>> step 0 <<<<<"<<endl;
    const double dnF1[11] = {-0.100,-0.400,-0.000, 0.00,-0.004, 0.00,-0.021,-0.050,-0.135,-0.200,-0.005};
    const double dnF2[11] = { 0.098, 0.370, 0.000, 0.00, 0.004, 0.00, 0.010, 0.035, 0.035, 0.160, 0.005};
    double par_Fh[7]; 
    par_Fh[0] = readParam(iBin,"Pro_Fh_shape", 0);
    par_Fh[1] = readParam(iBin,"Pro_Fh_shape", 1);
    par_Fh[2] = readParam(iBin,"Pro_Fh_shape", 2);
    par_Fh[3] = readParam(iBin,"Pro_Fh_shape", 3);
    par_Fh[4] = readParam(iBin,"Pro_Fh_shape", 4);
    par_Fh[5] = readParam(iBin,"Pro_Fh_shape", 5);
    par_Fh[6] = readParam(iBin,"Pro_Fh_shape", 6);
    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb[iBin] < 2.; iTrueAfb++){
        thisAfb = iTrueAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        afbTruth[iTrueAfb] = thisAfb;
        if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) {
            finAfb[iTrueAfb] = 0;
            continue;
        }
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        //    cout<<iTrueAfb<<" "<<thisAfb<<" "<<thisFh<<" "<<finAfb[iTrueAfb]<<endl;
        finAfb[iTrueAfb] = new TFile(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),"READ");
        h_setSummaryAfb[iTrueAfb] = (TH1F*)finAfb[iTrueAfb]->Get("h_setSummaryAfb");
        h_setSummaryAfb[iTrueAfb]->SetName(TString::Format("h_setSummaryAfb%+05.0f_fh%+05.0f",thisAfb*10000,thisFh*10000).Data());
        
        afbTrueValues[nAfbPoints]=thisAfb;
        nAfbPoints++;// for TGraph, connecting dots

        // Fit h_setSummaryAfb and keep the TF1
        const double upA[11] = { 0.60, 0.65, 0.18, 0.00, 0.20, 0.00, 0.20, 0.30, 0.18, 0.30, 0.10};
        const double dnA[11] = {-0.60,-0.65,-0.18, 0.00,-0.20, 0.00,-0.20,-0.30,-0.18,-0.30,-0.10};
        //afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+05.0f",thisAfb*10000).Data(),"gaus(0)", -1, 1);
//        afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+05.0f",thisAfb*10000).Data(),"[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4", dnA[iBin],upA[iBin]);
//        afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+05.0f",thisAfb*10000).Data(),"gaus(0)+[3]*x+[4]*x**2", dnA[iBin],upA[iBin]);
        afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+05.0f",thisAfb*10000).Data(),"gaus(0)", dnA[iBin],upA[iBin]);
//        afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+05.0f",thisAfb*10000).Data(),"[0]*exp(-0.5*((x-[1])/[2])**2)", dnA[iBin],upA[iBin]);
        afbPDFs[iTrueAfb]->SetParameter(1,0.);
        afbPDFs[iTrueAfb]->SetParameter(2,0.01);
//        afbPDFs[iTrueAfb]->SetParameter(4,0.02);
        afbPDFs[iTrueAfb]->SetParLimits(0,0,50);
        afbPDFs[iTrueAfb]->SetParLimits(2,0.00,1);
//        afbPDFs[iTrueAfb]->SetParLimits(3,0.,200);//spike at boundary
//        afbPDFs[iTrueAfb]->SetParLimits(4,0.001,0.02);//spike at boundary

	      TCanvas canvas("canvas");
        //h_setSummaryAfb[iTrueAfb]->Rebin(5);
        h_setSummaryAfb[iTrueAfb]->Scale(1./(h_setSummaryAfb[iTrueAfb]->GetSumOfWeights()*h_setSummaryAfb[iTrueAfb]->GetBinWidth(1)));//Normalize
//        h_setSummaryAfb[iTrueAfb]->Scale(1./h_setSummaryAfb[iTrueAfb]->Integral("width"));//Normalize
//        h_setSummaryAfb[iTrueAfb]->Scale(1./h_setSummaryAfb[iTrueAfb]->GetEntries());
        h_setSummaryAfb[iTrueAfb]->GetXaxis()->SetTitle("Measured A_{FB}");
        h_setSummaryAfb[iTrueAfb]->Draw("PE");
//	      lala->DrawLatexNDC(0.15,0.90,TString::Format("h_setSummaryAfb%+05.0f_fh%+05.0f",thisAfb*10000,thisFh*10000).Data());
        int r = h_setSummaryAfb[iTrueAfb]->Fit(afbPDFs[iTrueAfb],"RLQ");// Likelihood fit
        afbPDFs[iTrueAfb]->Draw("SAME");
        latex2->DrawLatexNDC(0.61,0.45,TString::Format("Width(#sigma) = %.4f#pm%.4f",afbPDFs[iTrueAfb]->GetParameter(2), afbPDFs[iTrueAfb]->GetParError(2)).Data());
        canvas.SaveAs(TString::Format("%s/h_setSummaryAfb%+05.0f_fit.png",otoyspath.Data(),thisAfb*10000));
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        wspace->import(*afbPDFs[iTrueAfb]);
        wspace->writeToFile(TString::Format("%s/wspace_afb%+05.0f.root",otoyspath.Data(),thisAfb*10000),true);
        canvas.Close();
        finAfb[iTrueAfb]->Close();
        if (r != 0) { cout<<"status = "<<r<<endl; }
    }
    cout<<">>>>> step 1 <<<<<"<<endl;
    TFile *finFh[maxPoints];
    TH1F  *h_setSummaryFh[maxPoints] ;
    thisAfb = 0.;
    thisFh =0.01;
    const double dnA1[11] = { 0.165, 0.078, 0.164, 0.00, 0.022, 0.00, 0.024, 0.075, 0.140, 0.270, 0.022};
    double par_Afb[4]; 
    par_Afb[0] = readParam(iBin,"Pro_Afb_shape", 0);
    par_Afb[1] = readParam(iBin,"Pro_Afb_shape", 1);
    par_Afb[2] = readParam(iBin,"Pro_Afb_shape", 2);
    par_Afb[3] = readParam(iBin,"Pro_Afb_shape", 3);
    for(int iTrueFh = 0; iTrueFh*stepSizeFh[iBin] < dataFhBand[iBin]; iTrueFh++){
        thisFh = iTrueFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        fhTruth[iTrueFh] = thisFh;
        if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) {
            finFh[iTrueFh] = 0;
            continue;
        }
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        finFh[iTrueFh] = new TFile(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),"READ");
        h_setSummaryFh[iTrueFh] = (TH1F*)finFh[iTrueFh]->Get("h_setSummaryFh");
        h_setSummaryFh[iTrueFh]->SetName(TString::Format("h_setSummaryFh%+05.0f_afb%+05.0f",thisFh*10000,thisAfb*10000).Data());

        fhTrueValues[nFhPoints]=thisFh;
        nFhPoints++;// for TGraph, connecting dots
        
        // Fit h_setSummaryFh and keep the TF1
        const double upF[11] = { 2.50, 2.50, 0.80, 0.00, 0.90, 0.00, 0.90, 0.90, 0.90, 1.50, 0.40};
        const double dnF[11] = { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
        //fhPDFs[iTrueFh] = new TF1(TString::Format("fhPdf_%+05.0f",thisFh*10000).Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)+[5]*exp(-0.5*((x-0.001)/[6])**2)", 0,3);
        fhPDFs[iTrueFh] = new TF1(TString::Format("fhPdf_%+05.0f",thisFh*10000).Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)+[5]*exp(-0.5*((x-0.001)/[6])**2)", dnF[iBin],upF[iBin]);
        //fhPDFs[iTrueFh] = new TF1(TString::Format("fhPdf_%+05.0f",thisFh*10000).Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)", dnF[iBin],upF[iBin]);
        //fhPDFs[iTrueFh] = new TF1(TString::Format("fhPdf_%+05.0f",thisFh*10000).Data(),"gaus(0)+[3]", dnF[iBin],upF[iBin]);
        fhPDFs[iTrueFh]->SetParameter(2,0.1);
        fhPDFs[iTrueFh]->SetParameter(3,10);
        fhPDFs[iTrueFh]->SetParameter(4,0.002);
        fhPDFs[iTrueFh]->SetParameter(5,10);
        fhPDFs[iTrueFh]->SetParameter(6,0.002);
        fhPDFs[iTrueFh]->SetParLimits(0,0,10);
        fhPDFs[iTrueFh]->SetParLimits(2,0.02,1);
        fhPDFs[iTrueFh]->SetParLimits(3,0.,200);//spike at boundary
        fhPDFs[iTrueFh]->SetParLimits(4,0.001,0.02);//spike at boundary
        fhPDFs[iTrueFh]->SetParLimits(5,0.,200);//spike at boundary
        fhPDFs[iTrueFh]->SetParLimits(6,0.001,0.02);//spike at boundary
        
	      TCanvas canvas("canvas");
        h_setSummaryFh[iTrueFh]->Scale(1./(h_setSummaryFh[iTrueFh]->GetSumOfWeights()*h_setSummaryFh[iTrueFh]->GetBinWidth(1)));//Normalize
        //h_setSummaryFh[iTrueFh]->Scale(1./h_setSummaryFh[iTrueFh]->Integral("width"));//Normalize
        //h_setSummaryFh[iTrueFh]->Scale(1./h_setSummaryFh[iTrueFh]->GetEntries());
        h_setSummaryFh[iTrueFh]->GetXaxis()->SetTitle("Measured F_{H}");
        h_setSummaryFh[iTrueFh]->Draw("PE");
//	      lala->DrawLatexNDC(0.15,0.90,TString::Format("h_setSummaryFh%+05.0f_afb%+05.0f",thisFh*10000,thisAfb*10000).Data());
        int r = h_setSummaryFh[iTrueFh]->Fit(fhPDFs[iTrueFh],"RLQ");// Likelihood fit
        fhPDFs[iTrueFh]->Draw("SAME");
        latex2->DrawLatexNDC(0.61,0.45,TString::Format("Width(#sigma) = %.4f#pm%.4f",fhPDFs[iTrueFh]->GetParameter(2), fhPDFs[iTrueFh]->GetParError(2)).Data());
        canvas.SaveAs(TString::Format("%s/h_setSummaryfh%+05.0f_fit.png",otoyspath.Data(),thisFh*10000));
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        wspace->import(*fhPDFs[iTrueFh]);
        wspace->writeToFile(TString::Format("%s/wspace_fh%+05.0f.root",otoyspath.Data(),thisFh*10000),true);
        canvas.Close();
        finFh[iTrueFh]->Close();
        if (r != 0) { cout<<"status = "<<r<<endl; }
        
    }
    cout<<finAfb[1]<<"         1 "<<endl;
    cout<<">>>>> step 2 <<<<<"<<endl;

        // As PDFs are settled, build F&C interval by P(x_test|true)/P(x_test|best_x)
    thisAfb = 0.;
    thisFh =0.01;
    double buffMaxLikelihoodAfb[maxPoints];// find P(x_test|best_x) as a function of measured afb.
    double buffMaxLikelihoodFh[maxPoints];// find P(x_test|best_x) as a function of measured fh.
    for (int i = 0; i < maxPoints; i++) {
        buffMaxLikelihoodAfb[i] = -2;
        buffMaxLikelihoodFh[i] = -2;
    }
    for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        afbTruth[iAfb] = thisAfb;
        //if (afbTruth[iAfb] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu_best) for x=thisAfb=%+04f\n",afbTruth[iAfb]);
//        double test_a = 0;
        for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb[iBin] < 2.; iTrueAfb++){
            thisAfb = iTrueAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
            if ( thisAfb <  dataAfbLo[iBin] ) continue;
            if ( thisAfb > -dataAfbLo[iBin] ) continue;
            if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
            else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
            else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
            //if (finAfb[iTrueAfb] == 0 ) continue;
            if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),afbTruth[iTrueAfb]*10000,thisFh*10000).Data(),&fiBuff) != 0) {
                finAfb[iTrueAfb] = 0;
                continue;
            }
            TFile *f_wspace = new TFile(TString::Format("%s/wspace_afb%+05.0f.root",otoyspath.Data(),afbTruth[iTrueAfb]*10000));
            RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
            afbPDFs[iTrueAfb] = (TF1*)wspace->obj(TString::Format("afbPdf_%+05.0f",afbTruth[iTrueAfb]*10000).Data());
//            test_a += (*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin];
//            cout<<afbTruth[iAfb]<<"   "<<stepSizeAfb[iBin]<<"  "<<(*afbPDFs[iTrueAfb])(afbTruth[iAfb])<<"   ";
//            cout<<(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin]<<"   "<<buffMaxLikelihoodAfb[iAfb]<<endl;
            if ((*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin] > buffMaxLikelihoodAfb[iAfb]) buffMaxLikelihoodAfb[iAfb] = (*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin];
            f_wspace->Close();
        }
//        cout<<test_a<<" is 1 ?"<<endl;
    }
    cout<<">>>>> step 3 <<<<<"<<endl;
    thisAfb = 0.;
    thisFh =0.01;
    for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
        thisFh = iFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        fhTruth[iFh] = thisFh;
        //if (fhTruth[iFh] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu_best) for x=thisFh=%+04f\n",fhTruth[iFh]);
//        double test_b = 0;
        for(int iTrueFh = 0; iTrueFh*stepSizeFh[iBin] < dataFhBand[iBin]; iTrueFh++){
            thisFh = iTrueFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
            if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
            if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
            //if (finFh[iTrueFh] == 0 ) continue;
            if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,fhTruth[iTrueFh]*10000).Data(),&fiBuff) != 0) {
                finFh[iTrueFh] = 0;
                continue;
            }
            if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
            TFile *f_wspace = new TFile(TString::Format("%s/wspace_fh%+05.0f.root",otoyspath.Data(),fhTruth[iTrueFh]*10000));
            RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
            fhPDFs[iTrueFh] = (TF1*)wspace->obj(TString::Format("fhPdf_%+05.0f",fhTruth[iTrueFh]*10000).Data());
//            test_b += (*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin];
//            cout<<(*fhPDFs[iTrueFh])(fhTruth[iFh])<<"   ";
//            cout<<buffMaxLikelihoodFh[iFh]<<endl;
            if ((*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin] > buffMaxLikelihoodFh[iFh]) buffMaxLikelihoodFh[iFh] = (*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin];
            f_wspace->Close();
        }
//        cout<<test_b<<" is 1 ?"<<endl;
    }
    cout<<">>>>> step 4 <<<<<"<<endl;
    
        // Create likelihood ratio ordering and find F-C interval
    double afbMeasuredHi[maxPoints];
    double afbMeasuredLo[maxPoints];
    for (int i = 0; i < maxPoints; i++) {
        afbMeasuredHi[i] = -2.;
        afbMeasuredLo[i] = -2.;
    }
    cout<<"nAfbPoints = "<<nAfbPoints<<endl;
    thisAfb = 0.;
    thisFh =0.01;
    TH1F *h_LRAfb[maxPoints];
    TH1F *h_buffProbAsLRAfb = new TH1F("buffProbAsLRAfb","",1000,0,1);
    int countAfbMeasuredHiLo = 0;
    nAfbPoints = 0;
    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb[iBin] < 2.; iTrueAfb++){
        thisAfb = iTrueAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
        if ( thisAfb <  dataAfbLo[iBin] ) continue;
        if ( thisAfb > -dataAfbLo[iBin] ) continue;
        if (thisAfb < dnF1[iBin]) thisFh = par_Fh[0] + par_Fh[1]*thisAfb;
        else if (thisAfb < dnF2[iBin]) thisFh = par_Fh[2] + par_Fh[3]*thisAfb + par_Fh[4]*thisAfb*thisAfb;
        else thisFh = par_Fh[5] + par_Fh[6]*thisAfb;
        //if (finAfb[iTrueAfb] == 0) continue;
        afbTruth[iTrueAfb] = thisAfb;
        //if (afbTruth[iTrueAfb] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu)/P(x,mu_best) for mu=thisAfb=%+04f\n",afbTruth[iTrueAfb]);
        //if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/FH_afb%+05.0f_fh%+05.0f",otoyspath.Data(),afbTruth[iTrueAfb]*10000,thisFh*10000).Data(),&fiBuff) != 0) {
            finAfb[iTrueAfb] = 0;
            continue;
        }
        afbTrueValues[nAfbPoints]=afbTruth[iTrueAfb];
        nAfbPoints++;// for TGraph, connecting dots
        
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_afb%+05.0f.root",otoyspath.Data(),afbTruth[iTrueAfb]*10000));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        afbPDFs[iTrueAfb] = (TF1*)wspace->obj(TString::Format("afbPdf_%+05.0f",afbTruth[iTrueAfb]*10000).Data());
        // Create LR ordering if finAfb exists
        h_LRAfb[iTrueAfb] = new TH1F(TString::Format("h_LRAfb_trueAfb%+05.0f",afbTruth[iTrueAfb]*10000).Data(),"",2/stepSizeAfb[iBin],-1,1);
//        h_LRAfb[iTrueAfb] = new TH1F(TString::Format("h_LRAfb_trueAfb%+05.0f",afbTruth[iTrueAfb]*10000).Data(),"",1000,-1,1);
        //double buffMaxLikelihoodAfbRatio[maxPoints];
        //memset(buffMaxLikelihoodAfbRatio,-2,maxPoints*sizeof(double));
        for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
            thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
            if ( thisAfb <  dataAfbLo[iBin] ) continue;
            if ( thisAfb > -dataAfbLo[iBin] ) continue;
            if (afbTruth[iAfb]< -1 ) continue;
            h_LRAfb[iTrueAfb]->SetBinContent(iAfb,(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin]/buffMaxLikelihoodAfb[iAfb]);
////            h_LRAfb[iTrueAfb]->SetBinContent(iAfb,(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeFh[iBin]/buffMaxLikelihoodAfb[iAfb]);
        }
        RooWorkspace *wspace1 = new RooWorkspace("wspace1","wspace1");
        wspace1->import(*h_LRAfb[iTrueAfb]);
        wspace1->writeToFile(TString::Format("%s/LRAFB_afb%+05.0f.root",otoyspath.Data(),afbTruth[iTrueAfb]*10000),true);

        TFile *f_wspace2 = new TFile(TString::Format("%s/LRAFB_afb%+05.0f.root",otoyspath.Data(),afbTruth[iTrueAfb]*10000));
        RooWorkspace *wspace2 = (RooWorkspace*)f_wspace2->Get("wspace1");
        h_LRAfb[iTrueAfb] = (TH1F*)wspace2->obj(TString::Format("h_LRAfb_trueAfb%+05.0f",afbTruth[iTrueAfb]*10000).Data());
        // Find interval
        double buffCov = 0;
        h_buffProbAsLRAfb->Reset();
        for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
            thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
            if ( thisAfb <  dataAfbLo[iBin] ) continue;
            if ( thisAfb > -dataAfbLo[iBin] ) continue;
            if (afbTruth[iAfb] < -1) continue;
//            cout<<iAfb<<"   "<<h_LRAfb[iTrueAfb]->GetBinContent(iAfb)<<"  "<<(*afbPDFs[iTrueAfb])(afbTruth[iAfb])<<"  "<<stepSizeAfb[iBin]<<endl;
            h_buffProbAsLRAfb->Fill(h_LRAfb[iTrueAfb]->GetBinContent(iAfb),(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb[iBin]);
        }
        f_wspace->Close();
        for(int tbin = h_buffProbAsLRAfb->GetNbinsX(); tbin > 0; tbin--){
            buffCov += h_buffProbAsLRAfb->GetBinContent(tbin);
            if (buffCov > targetCoverage ){
                buffCov = h_buffProbAsLRAfb->GetBinCenter(tbin);// serve as lower bound of LR
                break;
            }
        }
	      TCanvas canvas("canvas");
        TLine *line0 = new TLine();
        canvas.Divide(1,2);
        canvas.cd(1);
        gStyle->SetOptStat(0);
        h_LRAfb[iTrueAfb]->GetXaxis()->SetTitle("Measured A_{FB}");
        h_LRAfb[iTrueAfb]->GetYaxis()->SetTitle("Likelihood Ratio");
        h_LRAfb[iTrueAfb]->Draw("P");
        line0->SetLineWidth(3);
        line0->SetLineStyle(2);
        line0->SetLineColor(2);
        line0->DrawLine( -1., buffCov, 1., buffCov );
        canvas.cd(2);
        h_buffProbAsLRAfb->GetXaxis()->SetTitle("Likelihood Ratio");
        h_buffProbAsLRAfb->GetYaxis()->SetTitle("Probability");
        h_buffProbAsLRAfb->Draw("P");
        line0->DrawLine( buffCov, 0, buffCov, h_buffProbAsLRAfb->GetMaximum());
        canvas.Update();
        canvas.SaveAs(TString::Format("%s/h_afb%+05.0f_buffProbAsLRAfb.png",otoyspath.Data(),afbTruth[iTrueAfb]*10000));
        canvas.Close();
        //printf("DEBUG\t\t: lowest LR for iTrueAfb=%d = %f\n", iTrueAfb, buffCov);
        afbMeasuredLo[countAfbMeasuredHiLo] = -2.;
        afbMeasuredHi[countAfbMeasuredHiLo] = -2.;
//        cout<<countAfbMeasuredHiLo<<" --> "<<afbMeasuredLo[countAfbMeasuredHiLo]<<endl;
        for(int iAfb = 0; iAfb*stepSizeAfb[iBin] < 2.; iAfb++){
            thisAfb = iAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
            if ( thisAfb <  dataAfbLo[iBin] ) continue;
            if ( thisAfb > -dataAfbLo[iBin] ) continue;
            if (afbTruth[iAfb] < -1) continue;
            if ( h_LRAfb[iTrueAfb]->GetBinContent(iAfb) < buffCov ) continue;
            //printf("DEBUG\t\t: h_LRAfb(%d) = %f\n",iAfb,h_LRAfb[iTrueAfb]->GetBinContent(iAfb));
            if ( afbMeasuredLo[countAfbMeasuredHiLo] > -1 && (h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-buffCov)*(h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)-buffCov) < 0  ){
                // #2, 4,
//            if ( afbMeasuredLo[countAfbMeasuredHiLo] > -1 &&  afbMeasuredHi[countAfbMeasuredHiLo] < -1. && (h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-buffCov)*(h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)-buffCov) < 0  ){
//                cout<<buffCov<<"  "<<h_LRAfb[iTrueAfb]->GetBinContent(iAfb)<<"  "<<h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)<<"   "<<afbTruth[iAfb]<<endl;
                // #9
//                if ( fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1))+fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)) > 0.2 ) continue;
                // #6 
//                if ( fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1))+fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)) > 0.60 ) continue;
//                if (afbTruth[iAfb] > 0.09) continue;
//                printf("INFO\t\t: Upper limit found at %f!\n",afbTruth[iAfb]);
                afbMeasuredHi[countAfbMeasuredHiLo] = afbTruth[iAfb];
            }
            if ( afbMeasuredLo[countAfbMeasuredHiLo] < -1 && h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1) > buffCov ){
//                printf("INFO\t\t:   Lower limit found at boundary at %f!\n", afbTruth[iAfb-1]);
                afbMeasuredLo[countAfbMeasuredHiLo] = afbTruth[iAfb-1];// boundary case
            } 
            if ( afbMeasuredLo[countAfbMeasuredHiLo] < -1 && (h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-buffCov)*(h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1)-buffCov) < 0  ){
//                printf("INFO\t\t: Lower limit found at %f!\n", afbTruth[iAfb]);
                afbMeasuredLo[countAfbMeasuredHiLo] = afbTruth[iAfb];
            }
        }
//        cout<<countAfbMeasuredHiLo<<" --> "<<afbMeasuredLo[countAfbMeasuredHiLo]<<" ; afbTrueValues = "<<afbTrueValues[countAfbMeasuredHiLo]<<endl<<endl;
        f_wspace2->Close();
        //if (afbMeasuredHi[countAfbMeasuredHiLo] < -1) afbMeasuredHi[countAfbMeasuredHiLo] = 0.75;// Usually not used, simply push to the limit.
        if (afbMeasuredHi[countAfbMeasuredHiLo] < -1) continue;// Usually not used, simply push to the limit.
        countAfbMeasuredHiLo+=1;
        //countAfbMeasuredHiLo = countAfbMeasuredHiLo + 1;
    }
    cout<<">>>>> step 5 <<<<<"<<endl;
    TH1F *h_LRFh[maxPoints];
    TH1F *h_buffProbAsLRFh = new TH1F("buffProbAsLRFh","",1000,0,1);
    int countFhMeasuredHiLo = 0;
    thisAfb = 0.;
    thisFh =0.01;
    double fhMeasuredHi[maxPoints];
    double fhMeasuredLo[maxPoints];
    for (int i = 0; i < maxPoints; i++) {
        fhMeasuredHi[i] = -2.;
        fhMeasuredLo[i] = -2.;
    }
    nFhPoints = 0;
    for(int iTrueFh = 0; iTrueFh*stepSizeFh[iBin] < dataFhBand[iBin]; iTrueFh++){
        thisFh = iTrueFh*stepSizeFh[iBin]+stepSizeFh[iBin]/2.;
        if (thisFh < dnA1[iBin]) thisAfb = par_Afb[0] + par_Afb[1]*thisFh;
        if (thisFh > dnA1[iBin]) thisAfb = par_Afb[2] + par_Afb[3]*thisFh;
        fhTruth[iTrueFh] = thisFh;
        //if (fhTruth[iTrueFh] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu)/P(x,mu_best) for mu=thisFh=%+04f\n",fhTruth[iTrueFh]);
        //cout<<"iTrueFh = "<<iTrueFh<<";  fhTruth = "<<fhTruth[iTrueFh]<<endl;
        if (stat(TString::Format("%s/AFB_afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,fhTruth[iTrueFh]*10000).Data(),&fiBuff) != 0) {
            finFh[iTrueFh] = 0;
            continue;
        }
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        fhTrueValues[nFhPoints]=fhTruth[iTrueFh];
        nFhPoints++;// for TGraph, connecting dots
        //cout<<"5.4"<<endl;
        
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_fh%+05.0f.root",otoyspath.Data(),fhTruth[iTrueFh]*10000));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        fhPDFs[iTrueFh] = (TF1*)wspace->obj(TString::Format("fhPdf_%+05.0f",fhTruth[iTrueFh]*10000).Data());
        // Create LR ordering if finFh exists
        h_LRFh[iTrueFh] = new TH1F(TString::Format("h_LRFh_trueFh%+05.0f",fhTruth[iTrueFh]*10000).Data(),"",3/stepSizeFh[iBin],0,3);
        //double buffMaxLikelihoodFhRatio[maxPoints];
        //memset(buffMaxLikelihoodFhRatio,-2,maxPoints*sizeof(double));
        for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
            if (fhTruth[iFh]< -1 ) continue;
            h_LRFh[iTrueFh]->SetBinContent(iFh,(*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin]/buffMaxLikelihoodFh[iFh]);
            //cout<<iFh<<" "<<(*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin]<<"  "<<buffMaxLikelihoodFh[iFh]<<endl;
        }
        RooWorkspace *wspace1 = new RooWorkspace("wspace1","wspace1");
        wspace1->import(*h_LRFh[iTrueFh]);
        wspace1->writeToFile(TString::Format("%s/LRFH_fh%+05.0f.root",otoyspath.Data(),fhTruth[iTrueFh]*10000),true);

        TFile *f_wspace2 = new TFile(TString::Format("%s/LRFH_fh%+05.0f.root",otoyspath.Data(),fhTruth[iTrueFh]*10000));
        RooWorkspace *wspace2 = (RooWorkspace*)f_wspace2->Get("wspace1");
        h_LRFh[iTrueFh] = (TH1F*)wspace2->obj(TString::Format("h_LRFh_trueFh%+05.0f",fhTruth[iTrueFh]*10000).Data());

        // Find interval
        double buffCov = 0;
        h_buffProbAsLRFh->Reset();
        for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
            if (fhTruth[iFh] < -1) continue;
                //cout<<iFh<<"   "<<h_LRFh[iTrueFh]->GetBinContent(iFh)<<"  "<<(*fhPDFs[iTrueFh])(fhTruth[iFh])<<"  "<<stepSizeFh[iBin]<<endl;
            h_buffProbAsLRFh->Fill(h_LRFh[iTrueFh]->GetBinContent(iFh),(*fhPDFs[iTrueFh])(fhTruth[iFh])*stepSizeFh[iBin]);
        }
        f_wspace->Close();
//        h_LRFh[iTrueFh]->SaveAs(TString::Format("%s/h_LRFh_fh%+05.0f.cc",otoyspath.Data(),fhTruth[iTrueFh]*10000));
//        h_buffProbAsLRFh->SaveAs(TString::Format("%s/h_buffProbAsLRFh_fh%+05.0f.cc",otoyspath.Data(),fhTruth[iTrueFh]*10000));
//        TFile *fff = new TFile(TString::Format("%s/hsimple_fh%+05.0f.root",otoyspath.Data(),fhTruth[iTrueFh]*10000));
//        h_LRFh[iTrueFh]->Write(TString::Format("h_LRFh_fh%+05.0f",fhTruth[iTrueFh]*10000));
//        h_buffProbAsLRFh->Write(TString::Format("h_buffProbAsLRFh_fh%+05.0f",fhTruth[iTrueFh]*10000));
//        fff->Close();
        for(int tbin = h_buffProbAsLRFh->GetNbinsX(); tbin > 0; tbin--){
            buffCov += h_buffProbAsLRFh->GetBinContent(tbin);
            if (buffCov > targetCoverage){
                buffCov = h_buffProbAsLRFh->GetBinCenter(tbin);// serve as lower bound of LR
                break;
            }
        }
	      TCanvas canvas("canvas");
        TLine *line0 = new TLine();
        canvas.Divide(1,2);
        canvas.cd(1);
        gStyle->SetOptStat(0);
        h_LRFh[iTrueFh]->GetXaxis()->SetTitle("Measured F_{H}");
        h_LRFh[iTrueFh]->GetYaxis()->SetTitle("Likelihood Ratio");
        h_LRFh[iTrueFh]->Draw("P");
        line0->SetLineWidth(3);
        line0->SetLineStyle(2);
        line0->SetLineColor(2);
        line0->DrawLine( 0., buffCov, 3., buffCov );
        canvas.cd(2);
        h_buffProbAsLRFh->GetXaxis()->SetTitle("Likelihood Ratio");
        h_buffProbAsLRFh->GetYaxis()->SetTitle("Probability");
        h_buffProbAsLRFh->Draw("P");
        line0->DrawLine( buffCov, 0, buffCov, h_buffProbAsLRFh->GetMaximum());
        canvas.Update();
        canvas.SaveAs(TString::Format("%s/h_fh%+05.0f_buffProbAsLRFh.png",otoyspath.Data(),fhTruth[iTrueFh]*10000));
        canvas.Close();
        //cout<<buffCov<<endl;
        //printf("DEBUG\t\t: lowest LR for iTrueFh=%d = %f\n", iTrueFh, buffCov);
        fhMeasuredLo[countFhMeasuredHiLo] = -2.;
        fhMeasuredHi[countFhMeasuredHiLo] = -2.;
//        cout<<countFhMeasuredHiLo<<" --> "<<fhMeasuredLo[countFhMeasuredHiLo]<<endl;
        for(int iFh = 0; iFh*stepSizeFh[iBin] < dataFhBand[iBin]; iFh++){
            if (fhTruth[iFh] < -1) continue;
            if ( h_LRFh[iTrueFh]->GetBinContent(iFh) < buffCov ) continue;
            //printf("DEBUG\t\t: h_LRFh(%d) = %f\n",iFh,h_LRFh[iTrueFh]->GetBinContent(iFh));
            if ( fhMeasuredLo[countFhMeasuredHiLo] > -1 && (h_LRFh[iTrueFh]->GetBinContent(iFh)-buffCov)*(h_LRFh[iTrueFh]->GetBinContent(iFh+1)-buffCov) < 0  ){
//                if ( fabs(h_LRFh[iTrueFh]->GetBinContent(iFh)-h_LRFh[iTrueFh]->GetBinContent(iFh-1))+fabs(h_LRFh[iTrueFh]->GetBinContent(iFh)-h_LRFh[iTrueFh]->GetBinContent(iFh+1)) > 0.6 ) continue;
//                printf("INFO\t\t: Upper limit found at %f!\n",fhTruth[iFh]);
                fhMeasuredHi[countFhMeasuredHiLo] = fhTruth[iFh];
            }
            if ( fhMeasuredLo[countFhMeasuredHiLo] < -1 && h_LRFh[iTrueFh]->GetBinContent(iFh-1) > buffCov ){
                printf("INFO\t\t:   Lower limit found at boundary!\n");
                cout<<buffCov<<"   "<<fhTruth[iFh]<<"   "<<fhTruth[iFh-1]<<endl;
                fhMeasuredLo[countFhMeasuredHiLo] = fhTruth[iFh-1];// boundary case
            }
            if ( fhMeasuredLo[countFhMeasuredHiLo] < -1 && (h_LRFh[iTrueFh]->GetBinContent(iFh)-buffCov)*(h_LRFh[iTrueFh]->GetBinContent(iFh-1)-buffCov) < 0  ){
//                printf("INFO\t\t: Lower limit found at %f!\n", fhTruth[iFh]);
                fhMeasuredLo[countFhMeasuredHiLo] = fhTruth[iFh];
            }
        }
        f_wspace2->Close();
//        cout<<countFhMeasuredHiLo<<" --> "<<fhMeasuredLo[countFhMeasuredHiLo]<<" ; fhTrueValues = "<<fhTrueValues[countFhMeasuredHiLo]<<endl<<endl;
        //if (fhMeasuredHi[countFhMeasuredHiLo] < -1) fhMeasuredHi[countFhMeasuredHiLo] = 1.;// Usually not used, simply push to the limit.
        if (fhMeasuredHi[countFhMeasuredHiLo] < -1) continue;// Usually not used, simply push to the limit.
        countFhMeasuredHiLo+=1;
    }
    cout<<">>>>> step 6 <<<<<"<<endl;
//    cout<<"nAfbPoints = "<<nAfbPoints<<endl;
//    thisAfb = afb;
//    thisFh = fh;
//    int n_test = 0;
//    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb[iBin] < 2.; iTrueAfb++){
//        thisAfb = iTrueAfb*stepSizeAfb[iBin]-1.+stepSizeAfb[iBin]/2.;
//        afbTruth[iTrueAfb] = thisAfb;
//        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) {
//            continue;
//        }
//        afbTrueValues[n_test]=thisAfb;
//        cout<<afbTrueValues[n_test]<<endl;
//        n_test+=1; 
//    }
//    for (int ii = 0; ii < nAfbPoints; ii++) {
//        cout<<"nAfbPoints = "<<ii<<",  afbTrueValues = "<<afbTrueValues[ii]<<"; afbMeasuredHi = "<<afbMeasuredHi[ii]<<"; afbMeasuredLo = "<<afbMeasuredLo[ii]<<endl;
//    }
//    cout<<">>>>> step 7 <<<<<"<<endl;
//    cout<<"iTrueFh = 1"<<";  fhMeasuredHi = "<<fhMeasuredHi[1]<<endl;
//    cout<<"iTrueAfb = 1"<<";  afbMeasuredHi = "<<afbMeasuredHi[1]<<endl;
    // Make plots, fit, and find confidence interval
    TCanvas *canvas = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLatex *latex1 = new TLatex();
    latex1->SetNDC();
    TLine *line = new TLine();
    TGraph *g_fhIntervalHi = new TGraph(nFhPoints,fhMeasuredHi,   fhTrueValues);
    TGraph *g_fhIntervalLo = new TGraph(nFhPoints,fhMeasuredLo,   fhTrueValues);
    TGraph *g_afbIntervalHi = new TGraph(nAfbPoints,afbMeasuredHi,afbTrueValues);
    TGraph *g_afbIntervalLo = new TGraph(nAfbPoints,afbMeasuredLo,afbTrueValues);
    TF1 *f1_fhIntervalHi  = new TF1("f1_fhIntervalHi",  "[0]+[1]*x+[2]*x**2",0,3);
    TF1 *f1_fhIntervalLo  = new TF1("f1_fhIntervalLo",  "[0]+[1]*x+[2]*x**2",0,3);
    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x+[2]*x**2",-1,1);
    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x+[2]*x**2",-1,1);
    double Hi_Fhd[11]  = { 0.60, 0.60, 0.15, 3.00, 0.10, 5.00, 0.16, 0.14, 0.10, 0.20, 0.07};
    double Hi_Fhu[11]  = { 1.40, 1.20, 0.45, 3.00, 0.40, 5.00, 0.60, 0.60, 0.60, 1.00, 0.24};
    double Lo_Fhd[11]  = { 0.05, 0.40, 0.00, 3.00, 0.00, 5.00, 0.00, 0.02, 0.01, 0.05, 0.00};
    double Lo_Fhu[11]  = { 0.80, 1.00, 0.20, 3.00, 0.15, 5.00, 0.25, 0.30, 0.30, 0.55, 0.15};
    double Hi_Afbd[11] = {-0.05,-0.50,-0.05, 3.00,-0.08, 5.00,-0.07,-0.02,-0.06,-0.22,-0.04};
    double Hi_Afbu[11] = { 0.35, 0.40, 0.11, 3.00, 0.12, 5.00, 0.15, 0.13, 0.10, 0.24, 0.05};
    double Lo_Afbd[11] = {-0.30,-0.50,-0.15, 3.00,-0.17, 5.00,-0.15,-0.10,-0.10,-0.26,-0.05};
    double Lo_Afbu[11] = { 0.15, 0.40, 0.04, 3.00, 0.10, 5.00, 0.07, 0.05, 0.10, 0.10, 0.05};
    TFitResultPtr r_fhIntervalHi = g_fhIntervalHi->Fit("f1_fhIntervalHi","RS","",Hi_Fhd[iBin],Hi_Fhu[iBin]);
    TFitResultPtr r_fhIntervalLo = g_fhIntervalLo->Fit("f1_fhIntervalLo","RS","",Lo_Fhd[iBin],Lo_Fhu[iBin]);       
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","RS","",Hi_Afbd[iBin],Hi_Afbu[iBin]); 
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","RS","",Lo_Afbd[iBin],Lo_Afbu[iBin]);
//    TFitResultPtr r_fhIntervalHi = g_fhIntervalHi->Fit("f1_fhIntervalHi","RS","",Hi_Fh[iBin],fhMeasuredHi[nFhPoints-1]);
//    TFitResultPtr r_fhIntervalLo = g_fhIntervalLo->Fit("f1_fhIntervalLo","RS","",Lo_Fh[iBin],fhMeasuredLo[nFhPoints-1]);       
//    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","RS","",Hi_Afb[iBin],afbMeasuredHi[nAfbPoints-1]); 
//    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","RS","",Lo_Afb[iBin],afbMeasuredLo[nAfbPoints-1]);
//    int r_fhIntervalHi_1 = g_fhIntervalHi->Fit("f1_fhIntervalHi","RS","",Hi_Fh[iBin],fhMeasuredHi[nFhPoints-1]);
//    int r_fhIntervalLo_1 = g_fhIntervalLo->Fit("f1_fhIntervalLo","RS","",Lo_Fh[iBin],fhMeasuredLo[nFhPoints-1]);       
//    int r_afbIntervalHi_1= g_afbIntervalHi->Fit("f1_afbIntervalHi","RS","",Hi_Afb[iBin],afbMeasuredHi[nAfbPoints-1]); 
//    int r_afbIntervalLo_1= g_afbIntervalLo->Fit("f1_afbIntervalLo","RS","",Lo_Afb[iBin],afbMeasuredLo[nAfbPoints-1]);
//    cout<<afbMeasuredHi[0]<<"  "<<afbMeasuredLo[0]<<endl;
    double ParFCC_AFB_L[3], ParFCC_AFB_R[3], ParFCC_FH_L[3], ParFCC_FH_R[3];
    f1_fhIntervalHi->GetParameters(&ParFCC_FH_R[0]);
    f1_fhIntervalLo->GetParameters(&ParFCC_FH_L[0]);
    f1_afbIntervalHi->GetParameters(&ParFCC_AFB_R[0]);
    f1_afbIntervalLo->GetParameters(&ParFCC_AFB_L[0]);
    writeParam(iBin,"ParFCC_FH_R" , ParFCC_FH_R, 3 );
    writeParam(iBin,"ParFCC_FH_L" , ParFCC_FH_L, 3 );
    writeParam(iBin,"ParFCC_AFB_R" , ParFCC_AFB_R, 3 );
    writeParam(iBin,"ParFCC_AFB_L" , ParFCC_AFB_L, 3 );
    // Calculate the error, you must think about the case in which f1 not defined!
    double XPFh[2], XPAfb[2];
    XPAfb[0]=f1_afbIntervalHi->Eval(afb);
    XPAfb[1]=f1_afbIntervalLo->Eval(afb);
    if (f1_fhIntervalHi->Eval(fh) > 0.) XPFh[0]=f1_fhIntervalHi->Eval(fh);
    else XPFh[0]=0.;
    XPFh[1]=f1_fhIntervalLo->Eval(fh);
/////    if (fabs(f1_afbIntervalHi->Eval(afb)) <= 0.5*fh) XPAfb[0]=f1_afbIntervalHi->Eval(afb);
/////    else XPAfb[0]=-0.5*fh;
/////    if (fabs(f1_afbIntervalLo->Eval(afb)) <= 0.5*fh) XPAfb[1]=f1_afbIntervalLo->Eval(afb);
/////    else XPAfb[1]= 0.5*fh;
/////    if (f1_fhIntervalHi->Eval(fh) > dataFhLower[iBin] && f1_fhIntervalHi->Eval(fh) < fh) XPFh[0]=f1_fhIntervalHi->Eval(fh);
/////    else XPFh[0]=dataFhLower[iBin];
/////    if (f1_fhIntervalLo->Eval(fh) > dataFhLower[iBin] && f1_fhIntervalLo->Eval(fh) < 3.) XPFh[1]=f1_fhIntervalLo->Eval(fh);
/////    else XPFh[1]=3.;
    cout<<" FH_L = "<<XPFh[0]<<"   FH_H = "<<XPFh[1]<<endl;
    cout<<"AFB_L = "<<XPAfb[0]<<"  AFB_H = "<<XPAfb[1]<<endl;
//    cout<<" FH_Hi = "<<r_fhIntervalHi_1<<";  FH_Lo = "<<r_fhIntervalLo_1<<endl;
//    cout<<"AFB_Hi = "<<r_afbIntervalHi_1<<"; AFB_Lo = "<<r_afbIntervalLo_1<<endl;
    outFCErrFh[0] = r_fhIntervalHi.Get()  != 0 ? XPFh[0]-fh : 0.-fh;
    outFCErrFh[1] = r_fhIntervalLo.Get()  != 0 ? XPFh[1]-fh : 3.-fh;
    outFCErrAfb[0]= r_afbIntervalHi.Get() != 0 ? XPAfb[0]-afb : -0.5*fh-afb;
    outFCErrAfb[1]= r_afbIntervalLo.Get() != 0 ? XPAfb[1]-afb : 0.5*fh-afb;
////////  test  //////////
//    double XPFh_test[2];
//    XPFh_test[0]=f1_fhIntervalHi->Eval(0.2);
//    XPFh_test[1]=f1_fhIntervalLo->Eval(0.2);
//    double outFCErrFh_test[2];
//    outFCErrFh_test[0] = r_fhIntervalHi.Get()  != 0 ? XPFh_test[0]-0.2 : 0.-0.2;
//    outFCErrFh_test[1] = r_fhIntervalLo.Get()  != 0 ? XPFh_test[1]-0.2 : 3.-0.2;
////////  test  //////////
    cout<<r_afbIntervalLo.Get()<<endl;

//    const double dataAfbLo[11]  = {-0.500,-0.600,-0.150, 0.00,-0.150, 0.00,-0.200,-0.200,-0.250,-0.400,-0.100};
//    const double dataFhBand[11] = { 1.40,   1.50,  0.50, 0.00,  0.30, 0.00,  0.60,  0.50,  0.60,  1.00,  0.20};
    const double upAfb[11] = { 0.40,  0.80, 0.20, 0.00, 0.20, 0.00, 0.15, 0.15, 0.10, 0.40, 0.05};
    const double upFh[11]  = { 1.70,  1.60, 0.50, 0.00, 0.50, 0.00, 0.60, 0.70, 0.70, 1.10, 0.26};
    
    const double upAfbY[11]= { 0.50,  0.60, 0.10, 0.00, 0.15, 0.00, 0.12, 0.15, 0.20, 0.40, 0.10};
    const double upFhY[11] = { 1.40,  1.50, 0.30, 0.00, 0.25, 0.00, 0.40, 0.50, 0.60, 1.00, 0.20};
        // Fh plots
    gStyle->SetOptFit(0);
    g_fhIntervalHi->SetTitle("");
    g_fhIntervalHi->GetXaxis()->SetTitle("Measured F_{H}");
    g_fhIntervalHi->GetXaxis()->SetLimits(0.,upFh[iBin]);
    g_fhIntervalHi->GetYaxis()->SetTitle("True F_{H}");
    g_fhIntervalHi->GetYaxis()->SetRangeUser(0.,upFhY[iBin]);
//    g_fhIntervalLo->SetMarkerColor(2);
    g_fhIntervalHi->Draw("AP");
    g_fhIntervalLo->Draw("P SAME");
////////  Xcheck  //////////
//    double fhMeasuredTest[16] = {0.039, 0.201, 0.069, 0.255, 0.081, 0.279, 0.093, 0.315, 0.141, 0.369, 0.195, 0.441, 0.225, 0.501, 0.297, 0.567};
//    double fhTest[16]         = {0.081, 0.081, 0.123, 0.123, 0.153, 0.153, 0.183, 0.183, 0.243, 0.243, 0.303, 0.303, 0.351, 0.351, 0.423, 0.423};
//    TGraph *g_fhTest = new TGraph(16,fhMeasuredTest, fhTest);
//    g_fhTest->SetMarkerColor(2);
//    g_fhTest->Draw("P SAME");
////////  Xcheck  //////////
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(fh, 0.,fh,upFhY[iBin]);
////////  test  //////////
//    line->SetLineColor(4);
//    line->DrawLine(0.2, 0.,0.2,upFhY[iBin]);
//    latex->DrawLatexNDC(0.65,0.25,TString::Format("F_{H}_{t}=%.4f^{%+.4f}_{%+.4f}",0.2,outFCErrFh_test[1],outFCErrFh_test[0]).Data());
////////  test  //////////
    latex->DrawLatexNDC(0.65,0.25,TString::Format("F_{H}=%.4f^{%+.4f}_{%+.4f}",fh,outFCErrFh[1],outFCErrFh[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->SetTextFont(12);
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval2_fh_bin%d.pdf",plotpath.Data(),iBin));
//    canvas->Print(TString::Format("%s/FCConfInterval2_fh_bin%d_test.pdf",plotpath.Data(),iBin));
//    canvas->Print(TString::Format("%s/FCConfInterval2_fh_bin%d_Xcheck.pdf",plotpath.Data(),iBin));

    // Afb plots, remark true afb could be negative
    g_afbIntervalHi->SetTitle("");
    g_afbIntervalHi->GetXaxis()->SetTitle("Measured A_{FB}");
    g_afbIntervalHi->GetXaxis()->SetLimits(-upAfb[iBin],upAfb[iBin]);
    g_afbIntervalHi->GetYaxis()->SetTitle("True A_{FB}");
    g_afbIntervalHi->GetYaxis()->SetRangeUser(-upAfbY[iBin], upAfbY[iBin]);
    g_afbIntervalHi->Draw("AP");
//    g_afbIntervalLo->SetMarkerColor(2);
    g_afbIntervalLo->Draw("P SAME");
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(afb,-upAfbY[iBin],afb,upAfbY[iBin]);
    latex->DrawLatexNDC(0.65,0.25,TString::Format("A_{FB}=%.4f^{%+.4f}_{%+.4f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
	  latex1->DrawLatexNDC(0.15,0.90,TString::Format("CMS Preliminary"));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval2_afb_bin%d.pdf",plotpath.Data(),iBin));


    writeParam(iBin,"FCErrFh" ,outFCErrFh );
    writeParam(iBin,"FCErrAfb",outFCErrAfb);
//    writeParam(iBin,"FCErrFh_test" ,outFCErrFh_test );

//    return;

}//}}}



