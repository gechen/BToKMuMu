// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 
//
// -----------------------------------------------
//       Author: Geng CHEN <geng.chen@cern.ch> 
//       Created:   [2014-09-15 Mon 13:14] 
// -----------------------------------------------
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
//#include <regex>
#include "Math/MinimizerOptions.h"
#include "TROOT.h"
#include "tdrstyle.cc"

#include "RooGlobalFunc.h"
#include <stdlib.h>
#include "TMatrixDSym.h"
#include "RooMultiVarGaussian.h"
#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TPad.h> 
#include <TGaxis.h> 
#include <TLegend.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TPaveText.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TLorentzVector.h>
#include <TChainElement.h>

#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooChebychev.h> 
#include <RooGenericPdf.h> 
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooExtendPdf.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddition.h>
#include <RooMinuit.h>
#include <RooHist.h>
#include <RooWorkspace.h>

#include <iostream>
#include <fstream>

#include "tools.cc" 

using namespace std; 
using namespace RooFit;

// Tags configration
bool is7TeVCheck = false; // Using 2011 efficiency map.
TChain *ch=new TChain("tree");
TString plotpath="./plots";
TString iwspacepath="./RootFiles";
TString iCombBkgWspacepath=".";
TString owspacepath="./RootFiles";
TString idatacardpath="./fitParameters";
TString odatacardpath="./fitParameters";
TString ologpath=".";
TString summarypath;
double  scaleFactor=1.;

char LHCb_files[19][200] ={ "q2_0_to_1.dat",
							"q2_1_to_2.dat",
							"q2_2_to_3.dat",
							"q2_3_to_4.dat",
							"q2_4_to_5.dat",
							"q2_5_to_6.dat",
							"q2_6_to_7.dat",
							"q2_7_to_8.dat",
							"q2_11_to_11.dat",
							"q2_11_to_12.dat",
							"q2_15_to_16.dat",
							"q2_16_to_17.dat",
							"q2_17_to_18.dat",
							"q2_18_to_19.dat",
							"q2_19_to_20.dat",
							"q2_20_to_21.dat",
							"q2_21_to_22.dat",
							"q2_1_to_6.dat",
							"q2_15_to_22.dat"};
// Lumis.
////const float Lumi_JPSIX = 9.81, Lumi_Sig = 3296.81, Lumi_Jpsi = 18.58, Lumi_Psi = 212.50, Lumi_Data = 20.47;   // fb-1
////const float Lumi_Scale = 2.940;
double datasetLumi[6] = {20.47, 3296.81, 18.58, 212.50, 5951.1, 9.81};//data, signal MC, JpsiK MC, Psi'K MC, K*0mumu MC, JpsiX MC.

//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[15][200] = {"genQ2 <=  2.00 && genQ2 >  1.00",
                            "genQ2 <=  4.30 && genQ2 >  2.00",
                            "genQ2 <=  8.68 && genQ2 >  4.30",
                                "genQ2 <= 10.09 && genQ2 >  8.68",
							"genQ2 <= 12.86 && genQ2 > 10.09",
                                "genQ2 <= 14.18 && genQ2 > 12.86",
							"genQ2 <= 16.00 && genQ2 > 14.18",
                            "genQ2 <= 18.00 && genQ2 > 16.00",
                            "genQ2 <= 22.00 && genQ2 > 18.00",
                            "genQ2 <=  6.00 && genQ2 >  1.00",
							"(genQ2 > 1. && genQ2 < 8.68) || (genQ2 >= 14.18  && genQ2 < 22.) || (genQ2 >= 10.09  && genQ2 < 12.86)", //};
                            "genQ2 <=  4.30 && genQ2 >  1.00",
                            "genQ2 <=  8.68 && genQ2 >  6.00",
                            "genQ2 <=  6.00 && genQ2 >  4.30",
							"genQ2 <= 22.00 && genQ2 >  1.00"};
char Q2range[15][200] ={"Q2 <=  2.00 && Q2 >  1.00",
                        "Q2 <=  4.30 && Q2 >  2.00",
                        "Q2 <=  8.68 && Q2 >  4.30",
                            "Q2 <= 10.09 && Q2 > 8.68",
						"Q2 <= 12.86 && Q2 > 10.09",
                            "Q2 <= 14.18 && Q2 > 12.86",
						"Q2 <= 16.00 && Q2 > 14.18",
                        "Q2 <= 18.00 && Q2 > 16.00",
                        "Q2 <= 22.00 && Q2 > 18.00",
                        "Q2 <=  6.00 && Q2 >  1.00",
						"(Q2 > 1.  && Q2 < 8.68) || (Q2 >= 14.18 && Q2 < 22.) || (Q2 >= 10.09  && Q2 < 12.86)", //};
                        "Q2 <=  4.30 && Q2 >  1.00",
                        "Q2 <=  8.68 && Q2 >  6.00",
                        "Q2 <=  6.00 && Q2 >  4.30",
						"Q2 <= 22.00 && Q2 >  1.00"};
////char Q2String[15][200]={": 1.00 - 2.00 GeV^{2}",
////                        ": 2.00 - 4.30 GeV^{2}",
////                        ": 4.30 - 8.68 GeV^{2}",
////                            ": 8.68 - 10.09 GeV^{2}",
////						": 10.09 - 12.86 GeV^{2}",
////                            ": 12.86 - 14.18 GeV^{2}",
////						": 14.18 - 16.00 GeV^{2}",
////                        ": 16.00 - 18.00 GeV^{2}",
////                        ": 18.00 - 22.00 GeV^{2}",
////                        ": 1.00 - 6.00 GeV^{2}",
////						": 1.00 - 8.68, 10.09 - 12.86,", //};
////                        ": 1.00 - 4.30 GeV^{2}",
////                        ": 6.00 - 8.68 GeV^{2}",
////                        ": 4.30 - 6.00 GeV^{2}",
////						": 1.00 - 22.00 GeV^{2}"};
////char Q2String[15][200]={" (1.00-2.00) GeV^{2}",
////                        " (2.00-4.30) GeV^{2}",
////                        " (4.30-8.68) GeV^{2}",
////                            " (8.68-10.09) GeV^{2}",
////						" (10.09-12.86) GeV^{2}",
////                            " (12.86-14.18) GeV^{2}",
////						" (14.18-16.00) GeV^{2}",
////                        " (16.00-18.00) GeV^{2}",
////                        " (18.00-22.00) GeV^{2}",
////                        " (1.00-6.00) GeV^{2}",
////						" (1.00-8.68, 10.09-12.86,", //};
////                        "1.00 - 4.30 GeV^{2}",
////                        "6.00 - 8.68 GeV^{2}",
////                        "4.30 - 6.00 GeV^{2}",
////						"1.00 - 22.00 GeV^{2}"};
char Q2String[15][200]={"1 < #it{q^{2}} < 2 GeV^{2}",
                        "2 < #it{q^{2}} < 4.3 GeV^{2}",
                        "4.3 < #it{q^{2}} < 8.68 GeV^{2}",
                            "8.68 < #it{q^{2}} < 10.09 GeV^{2}",
                        "10.3 < #it{q^{2}} < 12.86 GeV^{2}",
                            "12.86 < #it{q^{2}} < 14.18 GeV^{2}",
                        "14.18 < #it{q^{2}} < 16 GeV^{2}",
                        "16 < #it{q^{2}} < 18 GeV^{2}",
                        "18 < #it{q^{2}} < 22 GeV^{2}",
                        "1 < #it{q^{2}} < 6 GeV^{2}",
                        "1 < #it{q^{2}} < 22 GeV^{2}",
                        "1.00 - 4.30 GeV^{2}",
                        "6.00 - 8.68 GeV^{2}",
                        "4.30 - 6.00 GeV^{2}",
						"1.00 - 22.00 GeV^{2}"};
char FIgString[12][20]={"(a)",
                        "(b)",
                        "(c)",
                        "()",
                        "(d)",
                        "()",
                        "(e)",
                        "(f)",
                        "(g)",
                        "(h)",
                        "(i)",
                        "()"};
char Q2Check_1[22][200] ={"Q2  <=  2.00 && Q2 >  1.00",
						"Q2  <=  3.00 && Q2 >  2.00",
						"Q2  <=  4.00 && Q2 >  3.00",
						"Q2  <=  5.00 && Q2 >  4.00",
						"Q2  <=  6.00 && Q2 >  5.00",
						"Q2  <=  7.00 && Q2 >  6.00",
						"Q2  <=  8.00 && Q2 >  7.00",
						"Q2  <=  8.68 && Q2 >  8.00",
						"Q2  <= 10.09 && Q2 >  8.68",
						"Q2  <= 11.00 && Q2 > 10.09",
						"Q2  <= 12.00 && Q2 > 11.00",
						"Q2  <= 12.86 && Q2 > 12.00",
						"Q2  <= 14.18 && Q2 > 12.86",
						"Q2  <= 15.00 && Q2 > 14.18",
						"Q2  <= 16.00 && Q2 > 15.00",
						"Q2  <= 17.00 && Q2 > 16.00",
						"Q2  <= 18.00 && Q2 > 17.00",
						"Q2  <= 19.00 && Q2 > 18.00",
						"Q2  <= 20.00 && Q2 > 19.00",
						"Q2  <= 21.00 && Q2 > 20.00",
						"Q2  <= 22.00 && Q2 > 21.00",
						"Q2  <= 23.00 && Q2 > 22.00"};
char Q2Check[15][200] ={"Q2  >  4.30 && Q2 <=  7.88",
						"Q2  >  4.30 && Q2 <=  8.28",
						"Q2  >  4.30 && Q2 <=  8.68",
						"Q2  >  4.30 && Q2 <=  9.08",
						"Q2  >  4.30 && Q2 <=  9.48",
						"Q2  > 10.68 && Q2 <= 12.26",
						"Q2  > 10.39 && Q2 <= 12.56",
						"Q2  > 10.09 && Q2 <= 12.86",
						"Q2  >  9.79 && Q2 <= 13.16",
						"Q2  >  9.49 && Q2 <= 13.46",
						"Q2  > 14.78 && Q2 <= 16.00",
						"Q2  > 14.48 && Q2 <= 16.00",
						"Q2  > 14.18 && Q2 <= 16.00",
						"Q2  > 13.88 && Q2 <= 16.00",
						"Q2  > 13.58 && Q2 <= 16.00"};
char CTL[2][100] = {"CosThetaL != -999", "0"};
double Q2rangedn[14] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 1.00 ,  1.00,  1.00, 6.00, 4.30};
double Q2rangeup[14] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 22.00 , 6.00 , 22.00,  4.30, 8.68, 6.00};
double genAfb   [11]={-0.00302  , -0.066   , 0.182    , 0.317    , 0.374    , 0.412    , 0.421    , 0.376    , -0.098   , 0.398    , 0.069};
double genAfberr[11]={ 0.000764 , 0.000285 , 0.000223 , 0.000383 , 0.000277 , 0.000421 , 0.000395 , 0.000422 , 0.000333 , 0.000146 , 0.000292};
double genFh    [11]={ 0.064596 , 0.791    , 0.649    , 0.524    , 0.454    , 0.399    , 0.369    , 0.341    , 0.7620   , 0.355    , 0.694};
double genFherr [11]={ 0.00185  , 0.000409 , 0.000271 , 0.000420 , 0.000287 , 0.000415 , 0.000369 , 0.000361 , 0.000240 , 0.000132 , 0.000258};
char mumuMassWindow[5][300] = { 
	" Mumumass > 0 ",
	" (Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)", // No JPsi && Psi(2S)
	" (Mumumass < 3.096916+3.*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr) \
	  || (Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr)",  // Only JPsi && Psi(2S)
	" Mumumass < 3.096916+2*Mumumasserr && Mumumass > 3.096916-2*Mumumasserr",  // Only JPsi
	" Mumumass < 3.686109+2*Mumumasserr && Mumumass > 3.686109-2*Mumumasserr"}; //Only Psi(2S)
char bmassWindow[3][100] = {
	" Bmass > 5.10 ",
	" Bmass > 5.10 ",
	" Bmass > 5.10 "};

double toUnboundedFh(double fh){
    return TMath::Tan( (fh - 1.5) * TMath::Pi() / 3. );
}
double toBoundedFh(double fh_ubd){
    return 3./2. + 3. * TMath::ATan(fh_ubd) / TMath::Pi() ;
}
double toUnboundedAfb(double afb, double fh){
    return TMath::Tan( afb * TMath::Pi() / tan( (fh - 1.5) * TMath::Pi() / 3.) ) ;
}
double toBoundedAfb(double afb_ubd, double fh_ubd){
    return ( 1. * TMath::ATan(afb_ubd) / TMath::Pi()) * ( 3./2. + 3. * TMath::ATan(fh_ubd) / TMath::Pi() );
}
double toFhUpErr(double fh, double fherr){
		return 3./2. + 3. * TMath::ATan( fh + fherr ) / TMath::Pi() - toBoundedFh(fh); 
}
double toFhDnErr(double fh, double fherr){
		return 3./2. + 3. * TMath::ATan( fh - fherr ) / TMath::Pi() - toBoundedFh(fh); 
}
double toAfbUpErr(double afb, double afberr, double fh, double fherr, double coAfbFh){
		return  sqrt( abs ( TMath::Power( ( 1. * toBoundedFh(fh) * TMath::ATan( afb + afberr ) / TMath::Pi() - toBoundedAfb(afb, fh) ), 2) + TMath::Power( ( 1. * TMath::ATan( afb ) * toFhUpErr(fh, fherr) / TMath::Pi() ), 2)  + 2. * (( 1. * toBoundedFh(fh) * TMath::ATan( afb + afberr ) / TMath::Pi() - toBoundedAfb(afb, fh) ) / afberr) * (( 1. * TMath::ATan( afb ) * toFhUpErr(fh, fherr) / TMath::Pi() ) / fherr) * coAfbFh ));
}
double toAfbDnErr(double afb, double afberr, double fh, double fherr, double coAfbFh){
		return  -sqrt( abs ( TMath::Power( ( 1. * toBoundedFh(fh) * TMath::ATan( afb - afberr ) / TMath::Pi() - toBoundedAfb(afb, fh) ), 2) + TMath::Power( ( 1. * TMath::ATan( afb ) * toFhDnErr(fh, fherr) / TMath::Pi() ), 2)  + 2. * (( 1. * toBoundedFh(fh) * TMath::ATan( afb - afberr ) / TMath::Pi() - toBoundedAfb(afb, fh) ) / afberr) * (( 1. * TMath::ATan( afb ) * toFhDnErr(fh, fherr) / TMath::Pi() ) / fherr) * coAfbFh ));
}
// Other tools
//_________________________________________________________________________________
std::string GetDirectory (const std::string& path)
{
   size_t found = path.find_last_of("/");
   return(path.substr(0, found));
}
bool scanAfbFhPositivePdf(double afb, double fh, bool fineScan=false)
{//{{{
    if ( fh < 1.5 && fabs(afb)<(fh*0.5) ) return true;

    // Create test function to find possible domain for fh/afb. CosThetaL as x.
    TString f1_format = "0.75*(1-[0])*(1-x*x) + 0.5*[0] + [1]*x";
    TF1 *f1_model = new TF1("f1_model", f1_format.Data(),-1.,1.);
    f1_model->SetTitle("PDF value;cos#theta_{L};");

    f1_model->SetParameter(0,fh);
    f1_model->SetParameter(1,afb);
    int nScanSteps = 1000;
    if (fineScan) nScanSteps *= 10;
    bool isPositivePDF=true;
    if (f1_model->Eval(1.) < 0 || f1_model->Eval(-1.) < 0){
        isPositivePDF = false;
    }else{
        for (int i = 0; i < nScanSteps; ++i) {//cosThetaL
            if (afb*2./nScanSteps*i-1. > 0) continue;
                if (f1_model->Eval(2./nScanSteps*i-1.) < 0.){
                    isPositivePDF = false;
                    break;
                }
            if (!isPositivePDF) break;
        }
    }

    return isPositivePDF;
}//}}}

void scanAfbFhPositivePdf()
{//{{{
    // Create test function to find possible domain for fh/afb. CosThetaL as x, CosThetaK as y.
    TCanvas *canvas= new TCanvas("canvas");
    TH2F *h2_minPdfValue = new TH2F("h2_minPdfValue","",200,-1,1,200,0,3);
    h2_minPdfValue->SetStats(false);
    h2_minPdfValue->SetXTitle("A_{FB}");
    h2_minPdfValue->SetYTitle("F_{H}");

    for( int xBin = 1; xBin <= h2_minPdfValue->GetNbinsX(); xBin++){//afb
        for( int yBin = 1; yBin <= h2_minPdfValue->GetNbinsY(); yBin++){//fh
            bool isPositivePDF = scanAfbFhPositivePdf((yBin-0.5)/h2_minPdfValue->GetNbinsY(),(2.*xBin-1.)/h2_minPdfValue->GetNbinsX()-1);
            if(isPositivePDF){
                h2_minPdfValue->SetBinContent(xBin,yBin,1);
            }else{
            }
        }
    }

    // Draw contour
    h2_minPdfValue->Draw("COL");
    canvas->Update();
    canvas->SaveSource(TString::Format("%s/scanAfbFhPositivePdf.cc",plotpath.Data()));
    canvas->Print(TString::Format("%s/scanAfbFhPositivePdf.pdf",plotpath.Data()));
    return;
}//}}}

TF2  *f2_fcn = NULL;
double model_2D(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(xx,yy);
}//}}}

TH2F *h2_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h2_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_fcn->GetNbinsY(); j++) {
            int gBin = h2_fcn->GetBin(i,j);
            //double x[2] = {h2_fcn->GetXaxis()->GetBinCenter(i),h2_fcn->GetYaxis()->GetBinCenter(j)};
            double measure  = h2_fcn->GetBinContent(gBin);
            double error    = h2_fcn->GetBinError(gBin);
            
            //// Naively test using center value
            //double func     = model_2D(x, par);//Take center value
            //double delta    = (measure-func)/error;
            //if (measure != 0) 
            //    f+=delta*delta;
            
            //// Real run using integral
            for (int k = 0; k < f2_fcn->GetNpar(); k++){//nPar MUST be the same value as f2_fcn
                f2_fcn->SetParameter(k,par[k]);
            }
            double xi = h2_fcn->GetXaxis()->GetBinLowEdge(i);
            double xf = h2_fcn->GetXaxis()->GetBinUpEdge(i);
            double yi = h2_fcn->GetYaxis()->GetBinLowEdge(j);
            double yf = h2_fcn->GetYaxis()->GetBinUpEdge(j);
            //f2_fcn->SetRange(xi,xf,yi,yf);
            //double minX, minY;
            //f2_fcn->GetMinimumXY(minX,minY);
            //if (f2_fcn->Eval(minX,minY) < 0){
            //    f += 100;
            //}else{
                f += pow( (f2_fcn->Integral(xi,xf,yi,yf)/(xf-xi)/(yf-yi)-measure)/error,2);
            //}
        }
    }
    //printf("FCN in calls = %f\n",f);
}//}}}
 
double readParam(int iBin, const char parName[], int iColumn)
{//{{{
	std::vector<double> output;
	char lineBuff[512];
	char *valBuff;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",idatacardpath.Data(),iBin),"r");
	while(fgets(lineBuff,512,fp) != NULL ){
		valBuff = strtok(lineBuff," ");
		if ( strcmp(valBuff,parName) == 0 ){
		//	printf("INFO: readParam, matched %15s!\n",valBuff);
			valBuff = strtok(NULL," ");
			while(valBuff != NULL){
			//	output.push_back(stof(valBuff));//stof if c++11 function, use other function
				output.push_back(std::atof(valBuff));
				valBuff = strtok(NULL," ");
			}
			break;
		}
		memset(lineBuff,' ',512*sizeof(char));
	}
	fclose(fp);
	
	if (iColumn < output.size() ){
	//	cout<<"INFO: readParam, get "<<parName<<"["<<iColumn<<"] = "<<output.at(iColumn)<<endl;
	//	printf("INFO: readParam, get %s[%d]= %f \n",parName,iColumn, output.at(iColumn));
		return output.at(iColumn);
	}else{
	//	printf("ERROR: readParam, empty column! Return 0.\n");
		return 0.;
	}
}//}}}
void writeParam(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
	struct stat fiBuff;
	FILE *fi = 0;
	if (stat(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),&fiBuff) == 0){
		rename(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
		fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"r");
	}else{
        fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"w");
	}
	
	bool parExist = false;
	char lineBuff[512];
	char *valBuff = 0;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),"w");
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
	remove(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
double readOutput(const char outfile[], int iBin, int Index, const char parName[], int iColumn)
{//{{{
	std::vector<double> output;
	char lineBuff[512];
	char *valBuff;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = NULL;
	fp = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt",outfile,iBin,iBin,Index),"r");
	if (fp == NULL) { return 0; }
	while(fgets(lineBuff,512,fp) != NULL ){
		valBuff = strtok(lineBuff," ");
		if ( strcmp(valBuff,parName) == 0 ){
			valBuff = strtok(NULL," ");
			while(valBuff != NULL){
				output.push_back(std::atof(valBuff));
				valBuff = strtok(NULL," ");
			}
			break;
		}
		memset(lineBuff,' ',512*sizeof(char));
	}
	fclose(fp);
	
	if (iColumn < output.size() ){
		return output.at(iColumn);
	}else{
		return 0.;
	}
}//}}}
void writeOutput(const char outfile[], int iBin, int Index, const char parName[], double *val, int nVal=3, bool overwrite=true)
{//{{{
	struct stat fiBuff;
	FILE *fi = 0;
	if (stat(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt",outfile,iBin,iBin,Index),&fiBuff) == 0){
		rename(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt",outfile,iBin,iBin,Index),TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt.temp",outfile,iBin,iBin,Index));
		fi = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt.temp",outfile,iBin,iBin,Index),"r");
	}else{
		fi = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt.temp",outfile,iBin,iBin,Index),"w");
	}
	
	bool parExist = false;
	char lineBuff[512];
	char *valBuff = 0;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt",outfile,iBin,iBin,Index),"w");
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
	remove(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_%d.txt.temp",outfile,iBin,iBin,Index));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
double readEffP(int iBin, const char parName[], int iColumn)
{//{{{
	std::vector<double> output;
	char lineBuff[512];
	char *valBuff;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("./fitParameters/EffPara%d.txt",iBin),"r");
	while(fgets(lineBuff,512,fp) != NULL ){
		valBuff = strtok(lineBuff," ");
		if ( strcmp(valBuff,parName) == 0 ){
			valBuff = strtok(NULL," ");
			while(valBuff != NULL){
				output.push_back(std::atof(valBuff));
				valBuff = strtok(NULL," ");
			}
			break;
		}
		memset(lineBuff,' ',512*sizeof(char));
	}
	fclose(fp);
	
	if (iColumn < output.size() ){
		return output.at(iColumn);
	}else{
		return 0.;
	}
}//}}}
void writeEffP(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
	struct stat fiBuff;
	FILE *fi = 0;
	if (stat(TString::Format("./fitParameters/EffPara%d.txt",iBin),&fiBuff) == 0){
		rename(TString::Format("./fitParameters/EffPara%d.txt",iBin),TString::Format("./fitParameters/EffPara%d.txt.temp",iBin));
		fi = fopen(TString::Format("./fitParameters/EffPara%d.txt.temp",iBin),"r");
	}else{
		fi = fopen(TString::Format("./fitParameters/EffPara%d.txt.temp",iBin),"w");
	}
	
	bool parExist = false;
	char lineBuff[512];
	char *valBuff = 0;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("./fitParameters/EffPara%d.txt",iBin),"w");
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
	remove(TString::Format("./fitParameters/EffPara%d.txt.temp",iBin));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
double readEffPOutput(const char outfile[], int iBin, int idex, int Index, const char parName[], int iColumn)
{//{{{
	std::vector<double> output;
	char lineBuff[512];
	char *valBuff;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = NULL;
	fp = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt",outfile,iBin,iBin,idex,Index),"r");
	if (fp == NULL) { return 0; }
	while(fgets(lineBuff,512,fp) != NULL ){
		valBuff = strtok(lineBuff," ");
		if ( strcmp(valBuff,parName) == 0 ){
			valBuff = strtok(NULL," ");
			while(valBuff != NULL){
				output.push_back(std::atof(valBuff));
				valBuff = strtok(NULL," ");
			}
			break;
		}
		memset(lineBuff,' ',512*sizeof(char));
	}
	fclose(fp);
	
	if (iColumn < output.size() ){
		return output.at(iColumn);
	}else{
		return 0.;
	}
}//}}}
void writeEffPOutput(const char outfile[], int iBin, int idex, int Index, const char parName[], double *val, int nVal=4, bool overwrite=true)
{//{{{
	struct stat fiBuff;
	FILE *fi = 0;
	if (stat(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt",outfile,iBin,iBin,idex,Index),&fiBuff) == 0){
		rename(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt",outfile,iBin,iBin,idex,Index),TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt.temp",outfile,iBin,iBin,idex,Index));
		fi = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt.temp",outfile,iBin,iBin,idex,Index),"r");
	}else{
		fi = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt.temp",outfile,iBin,iBin,idex,Index),"w");
	}
	
	bool parExist = false;
	char lineBuff[512];
	char *valBuff = 0;
	memset(lineBuff,' ',512*sizeof(char));
	FILE *fp = fopen(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt",outfile,iBin,iBin,idex,Index),"w");
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
	remove(TString::Format("./OutputValues/%s/bin%d/OutputValues%d_EffP%d_%d.txt.temp",outfile,iBin,iBin,idex,Index));
}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
void MultivariateGaussianTest(int iBin,  const char outfile[] = "accXrecoEff")
{
if (iBin != 0 && iBin != 1 && iBin != 9) {
//if (iBin != 20) {
    int dim = 7;
	RooArgList yVec;
    RooArgList muVec;
    int i,j;
    RooRealVar* y;
    RooRealVar* mu_y;
    for (i = 0; i < dim; i++) {
        char* name = Form("y%d", i);
		double mmean;
		mmean = readParam(iBin,"accXrecoEff",i);
		double emean;
		emean = readParam(iBin,"accXrecoEffErr",i);
        y = new RooRealVar(name, name, mmean-emean, mmean+emean);
        yVec.add(*y);
        char* mu_name = Form("mu_y%d",i);
        mu_y = new RooRealVar(mu_name, mu_name, readParam(iBin,"accXrecoEff",i));
        muVec.add(*mu_y);
    }
    TMatrixDSym cov(dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            double errPar = readParam(iBin,TString::Format("EffErr%d",i+1), j);
            cov(i,j) = errPar;
        }
    }
    RooMultiVarGaussian mvg("mvg", "mvg", yVec, muVec, cov);
	RooDataSet* meee = mvg.generate(yVec, 200);
	meee->RooDataSet::write(TString::Format("fitParameters/eff_par_%d.txt",iBin));
} else {
    int dim = 7;
	RooArgList yVec;
    RooArgList muVec;
    int i,j;
    RooRealVar* y;
    RooRealVar* mu_y;
    for (i = 0; i < dim; i++) {
        char* name = Form("y%d", i);
		double mmean;
		mmean = readParam(iBin,"reco",i);
		double emean;
		emean = readParam(iBin,"recoErr",i);
        y = new RooRealVar(name, name, mmean-emean, mmean+emean);
        yVec.add(*y);
        char* mu_name = Form("mu_y%d",i);
        mu_y = new RooRealVar(mu_name, mu_name, readParam(iBin,"reco",i));
        muVec.add(*mu_y);
    }
    TMatrixDSym cov(dim);
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
			double errPar = readParam(iBin,TString::Format("RecoErr%d",i+1), j);
			cov(i,j) = errPar;
        }
    }
    RooMultiVarGaussian mvg("mvg", "mvg", yVec, muVec, cov);
	RooDataSet* meee = mvg.generate(yVec, 200);
	meee->RooDataSet::write(TString::Format("fitParameters/reco_par_%d.txt",iBin));

    int dim_a = 3;
	RooArgList xVec;
    RooArgList muVec_a;
    RooRealVar* x;
    RooRealVar* mu_x;
    for (i = 1; i <= dim_a; i++) {
        char* name = Form("x%d", i);
		double mmean;
		mmean = readParam(iBin,"acc",i);
		double emean;
		emean = readParam(iBin,"accErr",i);
        x = new RooRealVar(name, name, mmean-emean, mmean+emean);
        xVec.add(*x);
        char* mu_name = Form("mu_x%d",i);
        mu_x = new RooRealVar(mu_name, mu_name, readParam(iBin,"acc",i));
        muVec_a.add(*mu_x);
    }
    TMatrixDSym cov_a(dim_a);
    for (i = 0; i < dim_a; i++) {
        for (j = 0; j < dim_a; j++) {
			double errPar = readParam(iBin,TString::Format("AccErr%d",i+1), j);
			cov_a(i,j) = errPar;
        }
    }
    RooMultiVarGaussian mvg_a("mvg_a", "mvg_a", xVec, muVec_a, cov_a);
	RooDataSet* meee_a = mvg_a.generate(xVec, 200);
	meee_a->RooDataSet::write(TString::Format("fitParameters/acc_par_%d.txt",iBin));
}

}
////////////////////////////////////////////////////////////////////////////////////////////
void PlotFCN( int iBin, const char outfile[] = "FCN")
{
	setTDRStyle();
	TGraph *gr_afb = new TGraph();
	TGraph *gr_fh  = new TGraph();
	double FCN = 999, fcn = 999;
	int Index = 0, index = -1, n = -1, NIndex = 1000;
	double Iafb, Ifh, afb, fh;
	do {
		index+=1;
		while ( ! (fcn = readOutput(outfile, iBin, index, "FCN", 0)) && index < NIndex) index+=1;
		afb = readOutput(outfile, iBin, index, "F_afb", 0);  
		fh  = readOutput(outfile, iBin, index, "F_fh", 0);
////		while ( ! (fcn = readOutput(outfile, iBin, index, "migrad", 1)) && index < NIndex) index+=1;
////		afb = readOutput(outfile, iBin, index, "afb_migrad", 0);  
////		fh  = readOutput(outfile, iBin, index, "fh_migrad", 0);
        if (fh > 1 ) 
            cout<<fh<<"  "<<index<<endl;
        if (afb > 0.5 || afb < -0.5 ) 
            cout<<afb<<"  "<<index<<endl;
		if (fcn < 0 && fcn > -1.e8) {
			n+=1;
			gr_afb->SetPoint(n, afb, fcn);
			gr_fh->SetPoint(n, fh, fcn);
		}
		if ( FCN > fcn ) { FCN = fcn; Index = index; }
	//	if ( FCN > fcn && (fabs(afb) < (fh / 2.))) { FCN = fcn; Index = index; }
	} while ( index < NIndex);
	Iafb   = readOutput(outfile, iBin, Index, "afb", 0);
	Ifh   = readOutput(outfile, iBin, Index, "fh", 0);
	afb    = readOutput(outfile, iBin, Index, "F_afb", 0);  
	fh    = readOutput(outfile, iBin, Index, "F_fh", 0);  
////	afb    = readOutput(outfile, iBin, Index, "afb_migrad", 0);  
////	fh    = readOutput(outfile, iBin, Index, "fh_migrad", 0);  
	TCanvas *c = new TCanvas("c","c",800,600);
//	c->SetTitle("FCN distribution of A_{FB}");
	gr_afb->GetXaxis()->SetTitle("A_{FB}");
	gr_afb->GetYaxis()->SetTitle("NLL");
	gr_afb->Draw("AP");
	c->Print(TString::Format("./plots/%s_FCN_afb_bin%d.pdf",outfile,iBin));
	c->Clear();
//	c->SetTitle("FCN distribution of F_{H}");
	gr_fh->GetXaxis()->SetTitle("F_{H}");
	gr_fh->GetYaxis()->SetTitle("NLL");
	gr_fh->Draw("AP");
	c->Print(TString::Format("./plots/%s_FCN_fh_bin%d.pdf",outfile,iBin));
	c->Clear();
	cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
//	cout<<"Iafb  = "<<Iafb<<"   Ifh = "<<Ifh<<endl;
	cout<<"afb = "<<afb<<endl;
	cout<<"fh  = "<<fh<<endl;
	double val[3]={0,0,0};
	val[0] = Iafb; val[1] = readOutput(outfile, iBin, Index, "afb", 1);
	writeParam(iBin, TString::Format("Iafb_%s",outfile),val);
	val[1]=0; val[2]=0;
	val[0] = Ifh; val[1] = readOutput(outfile, iBin, Index, "fh", 1);
	writeParam(iBin, TString::Format("Ifh_%s",outfile),val);
	val[1]=0; val[2]=0;
	val[0] = afb;
	writeParam(iBin, TString::Format("afb_%s",outfile),val);
	val[1]=0; val[2]=0;
	val[0] = fh;
	writeParam(iBin, TString::Format("fh_%s",outfile),val);
}

void Plot2DFCN( int iBin, const char outfile[] = "FCN")
{
	setTDRStyle();
	TH2D *gr = new TH2D();
	TH2D *gr_afb_fh = new TH2D();
	TH2D *gr_afb_fh_f = new TH2D();
	double FCN = 999, fcn = 999;
	int Index = 0, index = 0, n = -1, NIndex = 1000;
	double afb, fh;
	do {
		index+=1;
		while ( ! (fcn = readOutput(outfile, iBin, index, "FCN", 0)) && index < NIndex) index+=1;
		afb = readOutput(outfile, iBin, index, "F_afb", 0);
		fh  = readOutput(outfile, iBin, index, "F_fh", 0);
		if (fcn < 0 ) {
			n+=1;
			gr_afb_fh->SetBinContent(afb, fh, fcn);
		}
		if ( FCN > fcn ) { FCN = fcn; Index = index; }
	} while ( index < NIndex);
	for (float nn = 0, x =-1.5, y = 0; x <= 1.5; nn+=1, x+=0.005) {
		y = 2 * fabs(x);
		if (x > -0.6 && x < 0.6) gr->SetBinContent(x, y, FCN);
	}
	
	TCanvas *c = new TCanvas("c","c",800,600);
	gr_afb_fh->GetXaxis()->SetTitle("A_{FB}");
	gr_afb_fh->GetYaxis()->SetTitle("F_{H}");
	gr->GetYaxis()->SetTitleOffset(1.15);
	
	gStyle->SetPalette(1);
	gr_afb_fh->Draw(" COL ");
	gr->Draw("SAME COL ");
	
	afb = readOutput(outfile, iBin, Index, "F_afb", 0);   // best fitted values
	fh  = readOutput(outfile, iBin, Index, "F_fh", 0);
	FCN  = readOutput(outfile, iBin, Index, "FCN", 0);
	cout<<Index<<"  "<<FCN<<endl;
	gr_afb_fh_f->SetBinContent(afb, fh, FCN);
	gr_afb_fh_f->Draw(" SAME p ");
	TLatex *t1 = new TLatex();
	t1->SetTextColor(4);
	t1->DrawLatex( fh-0.001, FCN+1.,  TString::Format("%d", iBin));
	
	c->Print(TString::Format("./plots/%s_FCN_afb_fh_bin%d.pdf",outfile,iBin));
	c->Clear();

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotAfbFh_i( int iBin, int Index, const char outfile[] = "AfbFh2D")
{
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
	
	TGraphAsymmErrors *gr_afb_fh = new TGraphAsymmErrors();
	double fcn = 999;
	int index = 0, n = -1, NIndex = 1000;
	double afb[9], fh[9];
	double aa, bb;
	if (Index == -999) {
		do {
			index+=1;
			while ( ! (fcn = readOutput(outfile, iBin, index, "FCN", 0)) && index < NIndex) index+=1;
			n+=1;
		//	aa = readOutput(outfile, iBin, index, "afb", 0);  // initial values Unconstrained
		//	bb  = readOutput(outfile, iBin, index, "fh", 0);
			double ifh, iafb;
			iafb = readOutput(outfile, iBin, index, "afb", 0);  // initial values Constrained
			ifh  = readOutput(outfile, iBin, index, "fh", 0);
			aa = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
			bb = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
			if (fcn < 0 ) {
				gr_afb_fh->SetPoint(n, aa, bb);
			}
		} while ( index < NIndex);
	} else if (Index == -1) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
		//	afb[i] = readParam(ibin,TString::Format("Iafb_%s",outfile), 0);  // best initial values Unconstrained
		//	fh[i]  = readParam(ibin,TString::Format("Ifh_%s",outfile), 0);
			double iafb, ifh;
			iafb = readParam(ibin,TString::Format("Iafb_%s",outfile), 0);  // best initial values Constrained
			ifh  = readParam(ibin,TString::Format("Ifh_%s",outfile), 0);
			afb[i] = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
			fh[i]  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
			cout<<"ibin = "<<ibin<<" afb = "<<afb[i]<<" fh = "<<fh[i]<<endl;
			gr_afb_fh->SetPoint(i, afb[i], fh[i]);
		}
	}
	gr->SetMarkerColor(2);
	gr->SetMarkerSize(0.1);
	gr->SetMarkerStyle(7);
	gr->GetYaxis()->SetTitleOffset(1.15);
	gr_0->SetLineColor(4);
	gr_1->SetLineColor(4);
	gr_0->GetXaxis()->SetTitle("A_{FB}");
	gr_0->GetYaxis()->SetTitle("F_{H}");
	gr_0->GetXaxis()->SetRangeUser(-1.0, 1.0);
	gr_0->GetYaxis()->SetRangeUser(-0.5, 3.);
	gr_0->Draw("AL");	
	gr_1->Draw("L");	
	gr->Draw("P");	
	gr_afb_fh->Draw("P");
	gr_afb_fh->SetMarkerStyle(20);
	
	if (Index == -999) {
		TGraphAsymmErrors *gr_afb_fh_f = new TGraphAsymmErrors();
	//	aa = readParam(iBin,TString::Format("Iafb_%s",outfile), 0); // best initial values  Unconstrained
	//	bb = readParam(iBin,TString::Format("Ifh_%s",outfile), 0);
		double ifh, iafb;
		iafb = readParam(iBin,TString::Format("Iafb_%s",outfile), 0); // best initial values Constrained
		ifh  = readParam(iBin,TString::Format("Ifh_%s",outfile), 0);
		aa = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
		bb = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
		gr_afb_fh_f->SetPoint(0, aa, bb);
		gr_afb_fh_f->SetMarkerColor(3);
		gr_afb_fh_f->SetMarkerStyle(20);
		gr_afb_fh_f->Draw(" P ");
		t1->SetTextColor(4);
		t1->DrawLatex( aa+0.002, bb-0.001,  TString::Format("%d", iBin));
		t1->DrawLatex( -0.85, 0.25,  TString::Format("Afb = %f", aa));    // initial values
		t1->DrawLatex(  0.25, 0.25,  TString::Format("Fh = %f", bb));
		
		TLegend *leg =new TLegend(0.16,0.75,0.38,0.88,NULL,"brNDC");
		leg->AddEntry(gr_afb_fh,  " Succeeded initial values ","P");
		leg->AddEntry(gr_afb_fh_f," Best initial value ",    "P");
		leg->SetLineColor(1);
		leg->SetFillColor(0);
		leg->SetTextSize(0.02);
		leg->Draw();					
		c->Print(TString::Format("./plots/%s_afb_fh_bin%d_i.pdf",outfile,iBin));  
	} else if (Index == -1 ) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
			t1->SetTextColor(4);
			t1->DrawLatex(afb[i]+0.001, fh[i]-0.001, TString::Format("%d", ibin));
		}
		TLegend *leg_1 =new TLegend(0.16,0.79,0.38,0.88,NULL,"brNDC");
		leg_1->AddEntry(gr_afb_fh," Best initial value ",    "P");
		leg_1->SetLineColor(1);
		leg_1->SetFillColor(0);
		leg_1->SetTextSize(0.02);
		leg_1->Draw();					
		c->Print(TString::Format("./plots/%s_afb_fh_bin%d_i.pdf",outfile,iBin));
	}
	
	c->Clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotAfbFh_f( int iBin, int Index, const char outfile[] = "AfbFh2D")
{
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
	
	TGraphAsymmErrors *gr_afb_fh = new TGraphAsymmErrors();
	double fcn = 999;
	int index = 0, n = -1, NIndex = 1000;
	double afb[9], afberrL[9], afberrH[9], fh[9], fherrL[9], fherrH[9];
	double aa, bb;
	if (Index == -999) {
		do {
			index+=1;
			while ( ! (fcn = readOutput(outfile, iBin, index, "FCN", 0)) && index < NIndex) index+=1;
			n+=1;
			aa = readOutput(outfile, iBin, index, "F_afb", 0);  // best fitted values
			bb  = readOutput(outfile, iBin, index, "F_fh", 0);
			if (fcn < 0 ) {
				gr_afb_fh->SetPoint(n, aa, bb);
			}
		} while ( index < NIndex);
	} else if (Index == -1) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
			afb[i] = readParam(ibin,TString::Format("afb_%s",outfile), 0);   // best fitted values
			fh[i]  = readParam(ibin,TString::Format("fh_%s",outfile), 0);
			cout<<"ibin = "<<ibin<<" afb = "<<afb[i]<<" fh = "<<fh[i]<<endl;
			gr_afb_fh->SetPoint(i, afb[i], fh[i]);
		}
	} else if (Index == -2) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
			TString outname = outfile;
			if (outname == "angular_reco")  {
				afb[i] = readParam(ibin,"recoafb", 0);
				afberrL[i] = fabs(readParam(ibin,"recoafb", 1));
				afberrH[i] = fabs(readParam(ibin,"recoafb", 2));
				fh[i]  = readParam(ibin,"recofh", 0);
				fherrL[i] = fabs(readParam(ibin,"recofh", 1));
				fherrH[i] = fabs(readParam(ibin,"recofh", 2));
				fcn = readParam(ibin,"recoFCN", 0);
				cout<<"ibin = "<<ibin<<" recoafb = "<<afb[i]<<" recofh = "<<fh[i]<<endl;
			} else if (outname == "angular_reco_R")  {
				afb[i] = readParam(ibin,"recoRafb", 0);
				afberrL[i] = fabs(readParam(ibin,"recoRafb", 1));
				afberrH[i] = fabs(readParam(ibin,"recoRafb", 2));
				fh[i]  = readParam(ibin,"recoRfh", 0);
				fherrL[i] = fabs(readParam(ibin,"recoRfh", 1));
				fherrH[i] = fabs(readParam(ibin,"recoRfh", 2));
				fcn = readParam(ibin,"recoRFCN", 0);
				cout<<"ibin = "<<ibin<<" reco_R_afb = "<<afb[i]<<" reco_R_fh = "<<fh[i]<<endl;
			} else if (outname == "angular_reco_D")  {
				afb[i] = readParam(ibin,"recoDafb", 0);
				afberrL[i] = fabs(readParam(ibin,"recoDafb", 1));
				afberrH[i] = fabs(readParam(ibin,"recoDafb", 2));
				fh[i]  = readParam(ibin,"recoDfh", 0);
				fherrL[i] = fabs(readParam(ibin,"recoDfh", 1));
				fherrH[i] = fabs(readParam(ibin,"recoDfh", 2));
				fcn = readParam(ibin,"recoDFCN", 0);
				cout<<"ibin = "<<ibin<<" reco_D_afb = "<<afb[i]<<" reco_D_fh = "<<fh[i]<<endl;
			} else if (outname == "angular_gen_R") {
				afb[i] = readParam(ibin,"genafb_R", 0);
				afberrL[i] = fabs(readParam(ibin,"genafb_R", 1));
				afberrH[i] = fabs(readParam(ibin,"genafb_R", 2));
				fh[i]  = readParam(ibin,"genfh_R", 0);
				fherrL[i] = fabs(readParam(ibin,"genfh_R", 1));
				fherrH[i] = fabs(readParam(ibin,"genfh_R", 2));
				fcn = readParam(ibin,"genFCN_R", 0);
				cout<<"ibin = "<<ibin<<" genafb_R = "<<afb[i]<<" genfh_R = "<<fh[i]<<endl;
			} else if (outname == "angular_gen") {
				afb[i] = readParam(ibin,"genafb", 0);
				afberrL[i] = fabs(readParam(ibin,"genafb", 1));
				afberrH[i] = fabs(readParam(ibin,"genafb", 2));
				fh[i]  = readParam(ibin,"genfh", 0);
				fherrL[i] = fabs(readParam(ibin,"genfh", 1));
				fherrH[i] = fabs(readParam(ibin,"genfh", 2));
				fcn = readParam(ibin,"genFCN", 0);
				cout<<"ibin = "<<ibin<<" genafb = "<<afb[i]<<" genfh = "<<fh[i]<<endl;
			} else if (outname == "angular2D") {
				afb[i] = readParam(ibin,"afb", 0);
				//afberrL[i] = fabs(readParam(ibin,"afb", 1));
				//afberrH[i] = fabs(readParam(ibin,"afb", 2));
				afberrL[i] = fabs(readParam(ibin,"FCErrAfb", 0));
				afberrH[i] = fabs(readParam(ibin,"FCErrAfb", 1));
				fh[i]  = readParam(ibin,"fh", 0);
				//fherrL[i] = fabs(readParam(ibin,"fh", 1));
				//fherrH[i] = fabs(readParam(ibin,"fh", 2));
				fherrL[i] = fabs(readParam(ibin,"FCErrFh", 0));
				fherrH[i] = fabs(readParam(ibin,"FCErrFh", 1));
				fcn = readParam(ibin,"FCN", 0);
				cout<<"ibin = "<<ibin<<" afb = "<<afb[i]<<" fh = "<<fh[i]<<endl;
			}
			gr_afb_fh->SetPoint(i, afb[i], fh[i]);
			gr_afb_fh->SetPointError(i, afberrL[i], afberrH[i], fherrL[i], fherrH[i]);
		}
	}
	gr->SetMarkerColor(2);
	gr->SetMarkerSize(0.1);
	gr->SetMarkerStyle(7);
	gr->GetYaxis()->SetTitleOffset(1.15);
	gr_0->SetLineColor(4);
	gr_1->SetLineColor(4);
	gr_0->GetXaxis()->SetTitle("A_{FB}");
	gr_0->GetYaxis()->SetTitle("F_{H}");
//	gr_0->GetXaxis()->SetRangeUser(-0.020,0.020);  // gen_R fitted
//	gr_0->GetYaxis()->SetRangeUser(-0.005,0.070);
//	gr_0->GetXaxis()->SetRangeUser(-0.035,0.035);  // reco fitted
//	gr_0->GetYaxis()->SetRangeUser(-0.005,0.070);
	gr_0->GetXaxis()->SetRangeUser(-0.35,0.35);  // data fitted
	gr_0->GetYaxis()->SetRangeUser(-0.05,1.00);
    if (iBin == 9 || iBin == 0) {
        gr_0->GetXaxis()->SetRangeUser(-0.30,0.30);  // data fitted
        gr_0->GetYaxis()->SetRangeUser(-0.05,0.60);
    } else if (iBin == 1) {
        gr_0->GetXaxis()->SetRangeUser(-0.50,0.50);  // data fitted
	    gr_0->GetYaxis()->SetRangeUser(-0.05,1.00);
    } else if (iBin == 2 || iBin == 4){
        gr_0->GetXaxis()->SetRangeUser(-0.03,0.03);  // data fitted
        gr_0->GetYaxis()->SetRangeUser(-0.005,0.06);
    } else if ( iBin == 6 || iBin == 7 || iBin == 8 || iBin == 10){
        gr_0->GetXaxis()->SetRangeUser(-0.075,0.075);  // data fitted
        gr_0->GetYaxis()->SetRangeUser(-0.01,0.15);
    }
	gr_0->Draw("AL");	
	gr_1->Draw("L");	
	gr->Draw("P");	
	gr_afb_fh->Draw("P");
	gr_afb_fh->SetMarkerStyle(20);
	
	if (Index == -999) {
		TGraphAsymmErrors *gr_afb_fh_f = new TGraphAsymmErrors();
		aa = readParam(iBin,TString::Format("afb_%s",outfile), 0);   // best fitted values
		bb = readParam(iBin,TString::Format("fh_%s",outfile), 0);
		gr_afb_fh_f->SetPoint(0, aa, bb);
		gr_afb_fh_f->SetMarkerColor(4);
		gr_afb_fh_f->SetMarkerStyle(20);
		gr_afb_fh_f->Draw(" P ");
		t1->SetTextColor(4);
		t1->DrawLatex( aa+0.002, bb-0.001,  TString::Format("%d", iBin));
        double a1, a2, b1, b2;
		TGraphAsymmErrors *gr_afb_fh_d = new TGraphAsymmErrors();
        aa = readParam(iBin,"afb", 0);
        a1 = fabs(readParam(iBin,"afb", 1));
        a2 = fabs(readParam(iBin,"afb", 2));
        bb = readParam(iBin,"fh", 0);
        b1 = fabs(readParam(iBin,"fh", 1));
        b2 = fabs(readParam(iBin,"fh", 2));
		gr_afb_fh_d->SetPoint(0, aa, bb);
        gr_afb_fh_d->SetPointError(0, a1, a2, b1, b2);
		gr_afb_fh_d->SetMarkerColor(3);
		gr_afb_fh_d->SetLineColor(3);
		gr_afb_fh_d->SetMarkerStyle(30);
		gr_afb_fh_d->Draw(" P ");
////        if (iBin == 9 || iBin == 0) {
////		    t1->DrawLatex( -0.25, 0.01,  TString::Format("Afb = %f", aa));  // data fitted
////		    t1->DrawLatex(  0.05, 0.01,  TString::Format("Fh = %f", bb));
////        } else if (iBin == 1) {
////		    t1->DrawLatex( -0.40, 0.01,  TString::Format("Afb = %f", aa));  // data fitted
////		    t1->DrawLatex(  0.10, 0.01,  TString::Format("Fh = %f", bb));
////        } else if (iBin == 2 || iBin == 4 || iBin == 6 || iBin == 7 || iBin == 8 || iBin == 10){
////		    t1->DrawLatex( -0.07, 0.005,  TString::Format("Afb = %f", aa));  // data fitted
////		    t1->DrawLatex(  0.01, 0.005,  TString::Format("Fh = %f", bb));
////        }
	//	t1->DrawLatex( -0.018, 0.003,  TString::Format("Afb = %f", aa));  // gen_R fitted
	//	t1->DrawLatex(  0.003, 0.003,  TString::Format("Fh = %f", bb));
	//	t1->DrawLatex( -0.018, 0.001,  TString::Format("Afb = %f", aa));  // reco fitted
	//	t1->DrawLatex(  0.003, 0.001,  TString::Format("Fh = %f", bb));
	//	t1->DrawLatex( -0.30, 0.01,  TString::Format("Afb = %f", aa));  // data fitted
	//	t1->DrawLatex(  0.08, 0.01,  TString::Format("Fh = %f", bb));
		
		TLegend *leg =new TLegend(0.16,0.75,0.38,0.88,NULL,"brNDC");
	//	TLegend *leg =new TLegend(0.69,0.75,0.91,0.88,NULL,"brNDC");
		leg->AddEntry(gr_afb_fh,  " Succeeded scanning results ","P");
		leg->AddEntry(gr_afb_fh_f," Best scanning results ",    "P");
	    leg->AddEntry(gr_afb_fh_d," Data fitting results ", "leP");
		leg->SetLineColor(1);
		leg->SetFillColor(0);
		leg->SetTextSize(0.02);
		leg->Draw();					
		c->Print(TString::Format("./plots/%s_afb_fh_bin%d_f.pdf",outfile,iBin));  // fitted
	} else if (Index == -1 || Index == -2) {
		for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
			if (i == 3) ibin++;
			if (i == 4) ibin++;
			t1->SetTextColor(4);
			t1->DrawLatex(afb[i]+0.001, fh[i]-0.001, TString::Format("%d", ibin));
		}
		if (Index == -1) {
			TLegend *leg_1 =new TLegend(0.16,0.79,0.38,0.88,NULL,"brNDC");
			leg_1->AddEntry(gr_afb_fh," Best fitted result ", "P");
			leg_1->SetLineColor(1);
			leg_1->SetFillColor(0);
			leg_1->SetTextSize(0.02);
			leg_1->Draw();					
			c->Print(TString::Format("./plots/%s_afb_fh_bin%d_f.pdf",outfile,iBin));
		} else if (Index == -2){
		//	TLegend *leg_2 =new TLegend(0.16,0.79,0.40,0.88,NULL,"brNDC");
			TLegend *leg_2 =new TLegend(0.67,0.79,0.91,0.88,NULL,"brNDC");
		//	leg_2->AddEntry(gr_afb_fh," GEN(unfiltered MC) results ", "leP");
		//	leg_2->AddEntry(gr_afb_fh," GEN fitting results ", "leP");
		//	leg_2->AddEntry(gr_afb_fh," RECO fitting results ", "leP");
			leg_2->AddEntry(gr_afb_fh," Data fitting results ", "leP");
			leg_2->SetLineColor(1);
			leg_2->SetFillColor(0);
			leg_2->SetTextSize(0.03);
			leg_2->Draw();					
			//c->Print(TString::Format("./plots/%s_afb_fh_bin%d.pdf",outfile,iBin));
			c->Print(TString::Format("./plots/%s_FC_afb_fh_bin%d.pdf",outfile,iBin));
		}
	}
	
	c->Clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

