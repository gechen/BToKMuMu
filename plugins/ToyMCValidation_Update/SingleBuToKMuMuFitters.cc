// vim: set sw=3 sts=3 filetype=cpp fdm=marker et: 
// vim: tw=60 ts=2: 
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
#include "tdrstyle.C"

#include <RooGlobalFunc.h>
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
#include <TLegend.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TPaveText.h>
#include <TLatex.h>
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

// Lumis.
//const double Lumi_JPSIX = 9.81, Lumi_Sig = 3296.81, Lumi_Jpsi = 18.58, Lumi_Psi = 212.50, Lumi_Data = 20.47;   // fb-1
//const double Lumi_Scale = 2.940;
double datasetLumi[6] = {20.47, 3296.81, 18.58, 212.50, 5951.1, 9.81};//data, signal MC, JpsiK MC, Psi'K MC, K*0mumu MC, JpsiX MC.
//Constants, Fit results for efficiency, etc.. //{{{
char genQ2range[13][200] = {"genQ2 <=  2.00 && genQ2 >  1.00",
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
									"(genQ2 > 1. && genQ2 < 8.68) || (genQ2 >= 14.18  && genQ2 < 23.04) || (genQ2 >= 10.09  && genQ2 < 12.86)", //};
									"genQ2 <= 22.00 && genQ2 >  1.00"};
char Q2range[13][200] = {"Q2 <=  2.00 && Q2 >  1.00",
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
								"(Q2 > 1.  && Q2 < 8.68) || (Q2 >= 14.18 && Q2 < 23.04) || (Q2 >= 10.09  && Q2 < 12.86)", //};
								"Q2 <= 22.00 && Q2 >  1.00"};
char CTL[2][100] = {"CosThetaL != -999", "0"};
double Q2rangedn[13] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 1.00 ,  1.00,  1.00, 1.00};
double Q2rangeup[13] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 18.00 , 22.00 , 6.00 , 22.00, 23.04, 23.1};
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
	" Bmass > 5.00 "};
int isCDFcut = 4; // -1 for off, 1 for cdf, 2 for LHCb . 3 for 16Aug reOptimization

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


///////////////////////////////////////////////////////////////////////////////////////////////////////
//_________________________________________________________________________________
void bmass(int iBin, const char outfile[] = "bmass")
{//{{{
	bool test = false; 
	
	RooRealVar     DimuonPt("DimuonPt","DimuonPt", 0., 400.);
	RooRealVar     MupEta("MupEta","MupEta", -10., 10.);
	RooRealVar     MumEta("MumEta","MumEta", -10., 10.);
	
   RooRealVar     Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.10, 5.60);
	RooRealVar     Q2("Q2","q^{2}",1.0,22.);
	RooRealVar     Mumumass("Mumumass","M^{#mu#mu}",1.,10.);
	RooRealVar     Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
	RooDataSet     *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr),
	TString::Format("(%s)",Q2range[iBin]),0);	// data_v1

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
	delete c;

}//}}}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular_reco_bin(int iBin, double Iafb, double Ifh, int Index, const char outfile[] = "angular_reco")
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
	//	f_sigA_argset.add(RooArgSet(accP0, accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(accP1, accP2, accP3));
		f_sigA_argset.add(RooArgSet(recoP0, recoP1, recoP2, recoP3, recoP4, recoP5, recoP6));
		f_rec_format = "( accP1 *exp(-0.5*(((CosThetaL-accP2)/accP3)**2)) ) * ( recoP0 + recoP1 * CosThetaL + recoP2 * CosThetaL**2 + recoP3 * CosThetaL**3 + recoP4 * CosThetaL**4 + recoP5 * CosThetaL**5 + recoP6 * CosThetaL**6  )";
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
	} else {  // ScanFit
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 12-10-2014 N.A.
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
//	c->Print(TString::Format("./plots/%s_bin%d.png",outfile,iBin));
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
//	RooRealVar bkgCombL_c4("bkgCombL_c4","c4",0.,-3.,3.);
//	RooRealVar bkgCombL_c5("bkgCombL_c5","c5",0.,-3.,3.);
	RooArgSet f_bkgCombL_argset;
	switch (iBin) {
	//	case 1:
		case 2:
		case 7:
		case 8:
		case 9:
	//	case 10:
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
			break;
	   default:
		   f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c0));
	   // f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5,bkgCombL_c4));
		   break;
	}
	RooPolynomial f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
//	RooChebychev f_bkgCombL_P("f_bkgCombL_P","f_bkgCombL_P",CosThetaL,f_bkgCombL_argset);
	// Create peak background distribution
//	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}",  0.15,  0.03,  0.25);
	RooRealVar bkgGauss_mean1("bkgGauss_mean1","cos#theta_{l}",  0.2,  0.08,  0.25);
//	RooRealVar bkgGauss_mean2("bkgGauss_mean2","cos#theta_{l}",  0.2, -0.2,  0.5);
	RooRealVar bkgGauss_mean2("bkgGauss_mean2","cos#theta_{l}",  0.2,  0.05,  0.5);
	RooRealVar bkgGauss_mean3("bkgGauss_mean3","cos#theta_{l}",  0.3,  0.1,  0.35);
//	RooRealVar bkgGauss_mean3("bkgGauss_mean3","cos#theta_{l}",  0.3,  0.2,  0.5);
//	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.2, -0.9,  0.5);
	RooRealVar bkgGauss_mean4("bkgGauss_mean4","cos#theta_{l}", -0.6, -0.9, -0.5);
//	RooRealVar bkgGauss_mean5("bkgGauss_mean5","cos#theta_{l}", -0.2, -0.3,  0.1);
	RooRealVar bkgGauss_mean5("bkgGauss_mean5","cos#theta_{l}", -0.2, -0.3, -0.1);
//	RooRealVar bkgGauss_mean6("bkgGauss_mean6","cos#theta_{l}",  0.1, -0.3,  0.3);
	RooRealVar bkgGauss_mean6("bkgGauss_mean6","cos#theta_{l}",  0.01, -0.2,  0.2);
	RooRealVar bkgGauss_mean7("bkgGauss_mean7","cos#theta_{l}",  0.2, 0.01,  0.3);
//	RooRealVar bkgGauss_mean7("bkgGauss_mean7","cos#theta_{l}",  0.2, 0.09,  0.25);
//	RooRealVar bkgGauss_mean8("bkgGauss_mean8","cos#theta_{l}", -0.7, -0.99, -0.45);
	RooRealVar bkgGauss_mean8("bkgGauss_mean8","cos#theta_{l}", -0.7, -0.9, -0.5);
	RooRealVar bkgGauss_mean9("bkgGauss_mean9","cos#theta_{l}",  0.1,  0.0,  0.2);
	
//	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",  .20,  .01,  0.5);
	RooRealVar bkgGauss_sigma1("bkgGauss_sigma1","#sigma_{1}",  .10,  .06,  0.2);
//	RooRealVar bkgGauss_sigma2("bkgGauss_sigma2","#sigma_{2}",  .20,  .01,  2);
	RooRealVar bkgGauss_sigma2("bkgGauss_sigma2","#sigma_{2}",  .20,  .01,  2);
//	RooRealVar bkgGauss_sigma3("bkgGauss_sigma3","#sigma_{3}",  .10,  .05, 0.3);
	RooRealVar bkgGauss_sigma3("bkgGauss_sigma3","#sigma_{3}",  .10,  .06, 0.3);
//	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  1.0,  .01,  5);
	RooRealVar bkgGauss_sigma4("bkgGauss_sigma4","#sigma_{4}",  0.6,  .1,   1);
//	RooRealVar bkgGauss_sigma5("bkgGauss_sigma5","#sigma_{5}",  .20,  .01,  1);
	RooRealVar bkgGauss_sigma5("bkgGauss_sigma5","#sigma_{5}",  .20,  .01,0.4);
	RooRealVar bkgGauss_sigma6("bkgGauss_sigma6","#sigma_{6}",  .12,  .10,0.3);
//	RooRealVar bkgGauss_sigma7("bkgGauss_sigma7","#sigma_{7}",  .20,  .10,0.5);
	RooRealVar bkgGauss_sigma7("bkgGauss_sigma7","#sigma_{7}",  .20,  .01,0.4);
//	RooRealVar bkgGauss_sigma8("bkgGauss_sigma8","#sigma_{8}",  .50,  .01,  1);
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

	// Get data and apply unbinned fit
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass,CosThetaL),TString::Format("(%s) && ( (Bmass > 5.349 && Bmass < 5.458) || (Bmass < 5.209))",Q2range[iBin]),0);
	// RooDataSet is an unbinned dataset (a collection of points in N-dimensional space)
	RooDataSet *d = new RooDataSet("d","d",RooArgSet(CosThetaL));
	if ( iBin == 1 || iBin == 0 ) {
	for (int i = 0; i < 1; i++) {
		CosThetaL = 1 + 1 * .03;
		d->add(RooArgSet(CosThetaL));
		CosThetaL = -1 - 1 * .03;
		d->add(RooArgSet(CosThetaL));
	} }
	if (iBin == 10 || iBin == 8 || iBin == 7 || iBin == 2 || iBin == 4 || iBin == 9) {
	for (int i = 0; i < 4; i++) {
		CosThetaL = 1 + 1 * .03;
		d->add(RooArgSet(CosThetaL));
		CosThetaL = -1 - 1 * .03;
		d->add(RooArgSet(CosThetaL));
	} }
	RooDataSet* d1 = (RooDataSet*) data->reduce(RooArgSet(CosThetaL));
	d1->append(*d);	
	RooFitResult *f_fitresult = f_bkgCombL->fitTo(*d1,Save(kTRUE),Strategy(2),Minimizer("Minuit"),Warnings(-1),PrintEvalErrors(-1));
	f_fitresult->Print();
	if (f_fitresult->status() != 0) return;
   // comb.bkg. map 
   RooWorkspace *wspace1 = new RooWorkspace("wspace","wspace");
   wspace1->import(*f_bkgCombL);
   wspace1->writeToFile(TString::Format("%s/wspace_prior_bin%d.root",owspacepath.Data(),iBin),true);
	
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
	
	// clear
	delete t1;
	delete c;
	delete data;
	
	// Prepare datacard
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
		if ( iBin == 9) {
			val[0] = bkgGauss_mean1.getVal();val[1] = bkgGauss_mean1.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma1.getVal();val[1] = bkgGauss_sigma1.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 0 ) {
			val[0] = bkgGauss_mean2.getVal();val[1] = bkgGauss_mean2.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma2.getVal();val[1] = bkgGauss_sigma2.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 4) {
			val[0] = bkgGauss_mean5.getVal();val[1] = bkgGauss_mean5.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma5.getVal();val[1] = bkgGauss_sigma5.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		} else if (iBin == 1 ) {
			val[0] = bkgGauss_mean7.getVal();val[1] = bkgGauss_mean7.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma7.getVal();val[1] = bkgGauss_sigma7.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		} else if (iBin == 10 ) {
			val[0] = bkgGauss_mean9.getVal();val[1] = bkgGauss_mean9.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma9.getVal();val[1] = bkgGauss_sigma9.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 2 ) {
			val[0] = bkgGauss_mean6.getVal();val[1] = bkgGauss_mean6.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma6.getVal();val[1] = bkgGauss_sigma6.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 6 ) {
			val[0] = bkgGauss_mean4.getVal();val[1] = bkgGauss_mean4.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma4.getVal();val[1] = bkgGauss_sigma4.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 7 ) {
			val[0] = bkgGauss_mean8.getVal();val[1] = bkgGauss_mean8.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma8.getVal();val[1] = bkgGauss_sigma8.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
	   } else if (iBin == 8) {
			val[0] = bkgGauss_mean3.getVal();val[1] = bkgGauss_mean3.getError();
			writeParam(iBin, "bkgGauss_mean", val);
			val[0] = bkgGauss_sigma3.getVal();val[1] = bkgGauss_sigma3.getError();
			writeParam(iBin, "bkgGauss_sigma", val);
		}
		val[0] = bkg_frac.getVal();val[1] = bkg_frac.getError();
		writeParam(iBin, "bkg_frac", val);
	}
	
}//}}}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> angular2D_bin(int iBin, double Iafb, double Ifh, int Index, const char outfile[] = "angular2D")
{//{{{
	setTDRStyle();
	// Remark: You must use RooFit!! It's better in unbinned fit.
	//         Extended ML fit is adopted by Mauro, just follow!!
	//         Need some modification for accXrecoEff.
	cout<<endl<<"iBin = "<<iBin<<endl<<endl; 
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
	if (Index == -1) {
		f_ang_format = "( 0.75*(1-fh)*(1-CosThetaL*CosThetaL) + 0.5*fh + afb*CosThetaL )";
	} else {
		f_ang_format = "( 0.75*(1-( 3./2. + 3. * atan(fh) / TMath::Pi() ))*(1-CosThetaL*CosThetaL) + 0.5* ( 3./2. + 3. * atan(fh) / TMath::Pi() ) + (( 1. * atan(afb) / TMath::Pi()) * ( 3./2. + 3. * atan(fh) / TMath::Pi() )  )*CosThetaL )";
	}	
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
	// Create peak background distribution
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
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
//	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	
   f_fitresult->Print();
	if (f_fitresult->status() != 0 || f_fitresult->covQual() !=3) {
	std::vector<double> output;
	output.push_back(fh.getVal());
	output.push_back(fh.getError());
	output.push_back(afb.getVal());
	output.push_back(afb.getError());
	return output;
	}

	// Draw the frame on the canvas
	TCanvas *c = new TCanvas();
	TPad *tp1 = new TPad("tp1","",0.0,0.15,1.0,1.00);
	TPad *tp2 = new TPad("tp2","",0.0,0.00,1.0,0.15);
	tp1->Draw();
	tp2->Draw();
	tp1->cd();
	RooPlot* framemass = Bmass.frame();
	data->plotOn(framemass,RooFit::Name("data"), Binning(20)); 
	f.plotOn(framemass,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullmass = framemass->pullHist();
	f.plotOn(framemass, Components(f_sig),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framemass,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(2),LineWidth(2));
	f.plotOn(framemass,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(5),LineWidth(4),LineStyle(4));
	
	framemass->SetTitle("");
	framemass->SetTitleOffset(1.1,"Y");
	framemass->SetMinimum(0);
	framemass->SetMaximum(framemass->GetMaximum() * 1.2);
	framemass->Draw();
    
	TLegend *leg =new TLegend(0.17,0.69,0.35,0.86,NULL,"brNDC");
	leg->AddEntry("data"," Data "," PE ");
	leg->AddEntry("pdf"," Total P.d.f. "," L ");
	leg->AddEntry("sig"," Signal "," L ");
	leg->AddEntry("bkgComb"," Comb. bkg. "," L");
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->Draw();

  double FH, FHu, FHd;
  FH   = toBoundedFh( fh.getVal() );
  FHu  = toFhUpErr( fh.getVal(), fh.getError() ); 
  FHd  = toFhDnErr( fh.getVal(), fh.getError() ); 
  double AFB, AFBu, AFBd;
  AFB  = toBoundedAfb( afb.getVal(), fh.getVal() );
  AFBu = toAfbUpErr( afb.getVal(), afb.getError(), fh.getVal(), fh.getError(), f_fitresult->correlation(fh, afb) );
  AFBd = toAfbDnErr( afb.getVal(), afb.getError(), fh.getVal(), fh.getError(), f_fitresult->correlation(fh, afb) );

//    val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
//		val[1]= 3./2. + 3. * atan( fh.getVal() + fh.getError() ) / TMath::Pi() - val[0]; 
//		val[2]= 3./2. + 3. * atan( fh.getVal() - fh.getError() ) / TMath::Pi() - val[0]; 
//		double FH, FHu, FHd;
//		FH = val[0]; FHu = val[1]; FHd = val[2];
//		val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );
//		val[1]=   sqrt( pow( ( 1. * FH * atan( afb.getVal() + afb.getError() ) / TMath::Pi() - val[0] ), 2) + pow( ( 1. * atan( afb.getVal() ) * FHu / TMath::Pi() ), 2)  + 2. * ( 1. * FH * atan( afb.getVal() + afb.getError() ) / TMath::Pi() - val[0] ) / afb.getError() * ( 1. * atan( afb.getVal() ) * FHu / TMath::Pi() ) / fh.getError() * f_fitresult->correlation(afb, fh) );
//		val[2]= - sqrt( pow( ( 1. * FH * atan( afb.getVal() - afb.getError() ) / TMath::Pi() - val[0] ), 2) + pow( ( 1. * atan( afb.getVal() ) * FHd / TMath::Pi() ), 2)  + 2. * ( 1. * FH * atan( afb.getVal() - afb.getError() ) / TMath::Pi() - val[0] ) / afb.getError() * ( 1. * atan( afb.getVal() ) * FHd / TMath::Pi() ) / fh.getError() * f_fitresult->correlation(afb, fh) );
//		double AFB, AFBu, AFBd;
//		AFB = val[0]; AFBu = val[1]; AFBd = val[2];

	TPaveText* paveText = new TPaveText( 0.67, 0.65, 0.88, 0.88, "NDC" ); 
	paveText->SetBorderSize(0);
	paveText->SetFillColor(kWhite);
////////////////////////////////////// refit  //////////////////////////////////.........................
	paveText->AddText(Form("F_{H}  =%6.4f + %6.4f  %6.4f", FH, FHu, FHd)); 
	paveText->AddText(Form("A_{FB} =%6.4f + %6.4f  %6.4f", AFB, AFBu, AFBd)); 
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
	TPaveText* paveText_l = new TPaveText( 0.17, 0.52, 0.25, 0.62, "NDC" );
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
	c->Print(TString::Format("./plots/%s_bin%d_Index_%d.pdf",outfile,iBin,Index));
	
	// Draw projection to CosThetaL
	tp1->cd();
	RooPlot* framecosl = CosThetaL.frame(); 
	data->plotOn(framecosl,RooFit::Name("data"), Binning(20)); 
	f.plotOn(framecosl,RooFit::Name("pdf"), LineColor(1)); 
	RooHist *pullcosl = framecosl->pullHist();
	f.plotOn(framecosl, Components(f_sig),FillColor(2),VLines(), DrawOption("F"));
	f.plotOn(framecosl,RooFit::Name("sig"), Components(f_sig),LineStyle(2),LineColor(2),LineWidth(2));
	f.plotOn(framecosl,RooFit::Name("bkgComb"), Components(f_bkgComb),LineColor(5),LineWidth(4),LineStyle(4));
	
	framecosl->SetTitle("");
	framecosl->SetTitleOffset(1.1,"Y");
	framecosl->SetMinimum(0);
	framecosl->SetMaximum(framecosl->GetMaximum() * 1.2);
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
	c->Print(TString::Format("./plots/%s_cosl_bin%d_Index_%d.pdf",outfile,iBin,Index));

	delete c;
	delete t1;
	delete data;


// write output
	double val[4]={0,0,0,0};
	if (Index == -1) {
		val[0] = fh.getVal();val[1] = fh.getAsymErrorLo(); val[2] = fh.getAsymErrorHi();
//		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeParam(iBin, "fh", val, 3);
		val[0] = afb.getVal();val[1] = afb.getAsymErrorLo(); val[2] = afb.getAsymErrorHi();
//		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeParam(iBin, "afb",val, 3);
		val[1]=0; val[2]=0;
		val[0] = f_fitresult->minNll();
		writeParam(iBin, "FCN", val);
//	cout<<"AFB = "<<AFB<<"  + "<<AFBu<<" "<<AFBd<<endl;
//	cout<<"FH  = "<<FH<<"  + "<<FHu<<" "<<FHd<<endl;
		
//		val[0] = FH; val[1] = FHu; val[2] = FHd;
//		writeParam(iBin, "fh", val, 3);
//		val[0] = AFB; val[1] = AFBu; val[2] = AFBd;
//		writeParam(iBin, "afb",val, 3);
//		val[1]=0; val[2]=0;
//		val[0] = f_fitresult->minNll();
//		writeParam(iBin, "FCN", val);
	} else {
		val[0] = Iafb; val[1] = afb.getVal(); val[2] = afb.getError();
		writeOutput(outfile,iBin, Index, "afb", val);
		val[0] = Ifh;  val[1] = fh.getVal();  val[2] = fh.getError();
		writeOutput(outfile,iBin, Index, "fh", val);
    
		val[0] = FH; val[1] = FHu; val[2] = FHd;
		writeOutput(outfile,iBin, Index, "F_fh", val); 
		val[0] = AFB; val[1] = AFBu; val[2] = AFBd;
		writeOutput(outfile,iBin, Index, "F_afb", val);
		
//    val[0]=(1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() );  // Constrain
//		writeOutput(outfile,iBin, Index, "F_afb", val);
//		val[0]= 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();  // Constrain
//		writeOutput(outfile,iBin, Index, "F_fh", val); 
		val[1]=0; val[2]=0;
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		if (fcn < 0 ) {
			n+=1;
			gr_afb->SetPoint(n, afb, fcn);
			gr_fh->SetPoint(n, fh, fcn);
		}
		if ( FCN > fcn ) { FCN = fcn; Index = index; }
	} while ( index < NIndex);
	Iafb   = readOutput(outfile, iBin, Index, "afb", 0);
	Ifh   = readOutput(outfile, iBin, Index, "fh", 0);
	afb    = readOutput(outfile, iBin, Index, "F_afb", 0);  
	fh    = readOutput(outfile, iBin, Index, "F_fh", 0);  
	TCanvas *c = new TCanvas("c","c",800,600);
	gr_afb->GetXaxis()->SetTitle("A_{FB}");
	gr_afb->GetYaxis()->SetTitle("NLL");
	gr_afb->Draw("AP");
	c->Print(TString::Format("./plots/%s_FCN_afb_bin%d.pdf",outfile,iBin));
	c->Clear();
	gr_fh->GetXaxis()->SetTitle("F_{H}");
	gr_fh->GetYaxis()->SetTitle("NLL");
	gr_fh->Draw("AP");
	c->Print(TString::Format("./plots/%s_FCN_fh_bin%d.pdf",outfile,iBin));
	c->Clear();
	cout<<"Index = "<<Index<<"   FCN = "<<FCN<<endl;
//	cout<<"Iafb  = "<<Iafb<<"   Ifh = "<<Ifh<<endl;
	double afbU, afbD, fhU, fhD;
	afbU    = readOutput(outfile, iBin, Index, "F_afb", 1);  
	afbD    = readOutput(outfile, iBin, Index, "F_afb", 2);  
	fhU    = readOutput(outfile, iBin, Index, "F_fh", 1);  
	fhD    = readOutput(outfile, iBin, Index, "F_fh", 2);  
	cout<<"afb = "<<afb<<" + "<<afbU<<" "<<afbD<<endl;
	cout<<"fh  = "<<fh<<" + "<<fhU<<" "<<fhD<<endl;
	double val[3]={0,0,0};
	val[0] = Iafb; val[1] = readOutput(outfile, iBin, Index, "afb", 1);
	writeParam(iBin, TString::Format("Iafb_%s",outfile),val);
	val[1]=0; val[2]=0;
	val[0] = Ifh; val[1] = readOutput(outfile, iBin, Index, "fh", 1);
	writeParam(iBin, TString::Format("Ifh_%s",outfile),val);
	val[1]=afbU; val[2]=afbD;
	val[0] = afb;
	writeParam(iBin, TString::Format("afb_%s",outfile),val);
	val[1]=fhU; val[2]=fhD;
	val[0] = fh;
	writeParam(iBin, TString::Format("fh_%s",outfile),val);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotAfbFh_f( int iBin, int Index, const char outfile[] = "AfbFh2D")
{
	setTDRStyle();
	TGraphAsymmErrors *gr   = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_0 = new TGraphAsymmErrors();
	TGraphAsymmErrors *gr_1 = new TGraphAsymmErrors();
	TCanvas *c = new TCanvas("c","c",800,600);
	TLatex *t1 = new TLatex();
	for (double n = 0, x =-1., y = 0; x <= 1.; n+=1, x+=0.0005) {
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
			if (outname == "angular2D") {
				afb[i] = readParam(ibin,"afb", 0);
				afberrL[i] = fabs(readParam(ibin,"afb", 1));
				afberrH[i] = fabs(readParam(ibin,"afb", 2));
				fh[i]  = readParam(ibin,"fh", 0);
				fherrL[i] = fabs(readParam(ibin,"fh", 1));
				fherrH[i] = fabs(readParam(ibin,"fh", 2));
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
	gr_0->GetXaxis()->SetRangeUser(-0.35,0.35);  // data fitted
	gr_0->GetYaxis()->SetRangeUser(-0.05,1.00);
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
		gr_afb_fh_f->SetMarkerColor(3);
		gr_afb_fh_f->SetMarkerStyle(20);
		gr_afb_fh_f->Draw(" P ");
		t1->SetTextColor(4);
		t1->DrawLatex( aa+0.002, bb-0.001,  TString::Format("%d", iBin));
		t1->DrawLatex( -0.18, 0.01,  TString::Format("Afb = %f", aa));  // data fitted
		t1->DrawLatex(  0.03, 0.01,  TString::Format("Fh = %f", bb));
		
		TLegend *leg =new TLegend(0.16,0.75,0.38,0.88,NULL,"brNDC");
		leg->AddEntry(gr_afb_fh,  " Succeeded fitted results ","P");
		leg->AddEntry(gr_afb_fh_f," Best fitted results ",    "P");
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
			c->Print(TString::Format("./plots/%s_afb_fh_bin%d_f.png",outfile,iBin));
		} else if (Index == -2){
			TLegend *leg_2 =new TLegend(0.67,0.79,0.91,0.88,NULL,"brNDC");
			leg_2->AddEntry(gr_afb_fh," Data fitting results ", "leP");
			leg_2->SetLineColor(1);
			leg_2->SetFillColor(0);
			leg_2->SetTextSize(0.03);
			leg_2->Draw();					
			c->Print(TString::Format("./plots/%s_afb_fh_bin%d.pdf",outfile,iBin));
		}
	}
	
	c->Clear();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void angular2D( const char outfile[] = "angular2D")
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

	double SMFh[9]    = { 0.10,  0.07,    0.04,  0.00,   0.02,   0.02,   0.03, 0.05, 0.02};
	double SMFherr[7] = { 0.008, 0.001, 0.0005,  0.00, 0.0005, 0.0005, 0.0015};

	double x[9]   ={1.50, 3.15, 6.49,  11.475,  15.09, 17.0, 20.0, 3.5, 11.5};
	double xerr[9]={0.5, 1.15, 2.19,   1.385,   0.91,  1.0,  2.0, 2.5, 10.5};
	double yfh[9], yuerrfh[9],yderrfh[9], yafb[9], yuerrafb[9], yderrafb[9];
	double sysfh[9]  = { 0.0929, 0.0634, 0.0600, 0.0779, 0.0505, 0.0356, 0.0355, 0.0966, 0.0434 };
	double sysafb[9] = { 0.0666, 0.0202, 0.0298, 0.0365, 0.0151, 0.0152, 0.0134, 0.0281, 0.0102 };
	for (int i = 0; i < 9; i++) {
		yfh[i]  = 0; yuerrfh[i] = 0; yderrfh[i] = 0;
		yafb[i] = 0; yuerrafb[i]= 0; yderrafb[i]= 0;
	}
// Checkout input data
	for(int i = 0, ibin = 0; i < 9 && ibin < 11; i++, ibin++){
		if (i == 3) ibin++;
		if (i == 4) ibin++;
		cout<<"iBin = "<<ibin<<endl;
		yafb[i]      = readParam(ibin,"afb",0);
		yuerrafb[i]  = fabs(readParam(ibin,"afb",2));
		yderrafb[i]  = fabs(readParam(ibin,"afb",1));
		yfh[i]       = readParam(ibin,"fh",0);
		yuerrfh[i]   = fabs(readParam(ibin,"fh",2));
		yderrfh[i]   = fabs(readParam(ibin,"fh",1));
      if (sysfh[i] > yfh[i]) sysfh[i] = yfh[i];
		if (yuerrfh[i] > fabs(yfh[i])) { yderrfh[i] = fabs(yfh[i]);}
		else { yderrfh[i] = yuerrfh[i]; }
		printf("Afb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[i],yuerrafb[i],yderrafb[i]);
		printf("Fh [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfh[i], yuerrfh[i], yderrfh[i]);
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
	frame->SetAxisRange(-0.1,1.5,"Y");
	TGraphAsymmErrors *d_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,yderrfh,yuerrfh);
	d_fh->SetMarkerColor(1);
	d_fh->SetMarkerStyle(20);
	d_fh->Draw("P");
	TGraphAsymmErrors *sys_fh  = new TGraphAsymmErrors(7,x,yfh,xerr,xerr,sysfh,sysfh);
	sys_fh->SetFillColor(4);
	sys_fh->SetFillStyle(3244);
	sys_fh->Draw("2");
	TGraphAsymmErrors *s_fh  = new TGraphAsymmErrors(7,x,SMFh,xerr,xerr,SMFherr,SMFherr);
	s_fh->SetFillColor(2);
	s_fh->Draw("2");
	TGraphAsymmErrors *c_fh  = new TGraphAsymmErrors(2,conX,conYF,conXerr,conXerr,conYFerr,conYFerr);
	c_fh->SetFillColor(1);
	c_fh->SetFillStyle(3003);
	c_fh->Draw("2");
	TLine *tl2 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl2->SetLineColor(1);
	tl2->Draw();

	TLegend *leg =new TLegend(0.69,0.69,0.90,0.86);
	leg->AddEntry(d_fh," Data ","lep");
	leg->AddEntry(sys_fh," Systematic error ","f");
	leg->AddEntry(s_fh," SM prediction ","f");
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
	frame->SetAxisRange(-0.5,0.4,"Y");
	frame->Draw();
	TGraphAsymmErrors *d_afb = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,yderrafb,yuerrafb);
	d_afb->SetMarkerColor(1);
	d_afb->SetMarkerStyle(20);
	d_afb->Draw("P");
	TGraphAsymmErrors *sys_afb  = new TGraphAsymmErrors(7,x,yafb,xerr,xerr,sysafb,sysafb);
	sys_afb->SetFillColor(4);
	sys_afb->SetFillStyle(3244);
	sys_afb->Draw("2");
	TGraphAsymmErrors *c_afb  = new TGraphAsymmErrors(2,conX,conYA,conXerr,conXerr,conYAerr,conYAerr);
	c_afb->SetFillColor(1);
	c_afb->SetFillStyle(3003);
	c_afb->Draw("2");
   
	TLine *tl1 =new TLine(0.0, 0.0, 22.0, 0.0);
	tl1->SetLineColor(2);
	tl1->Draw();
	TLegend *leg_1 =new TLegend(0.69,0.69,0.90,0.86);
	leg_1->AddEntry(d_afb," Data ","lep");
	leg_1->AddEntry(sys_afb," Systematic error ","f");
	leg_1->SetLineColor(0);
	leg_1->SetFillColor(0);
	leg_1->SetTextSize(0.02);
	leg_1->Draw();
	
	t1->DrawLatex(.15,.90,TString::Format("CMS Preliminary"));
	t1->DrawLatex(.62,.90,TString::Format("Data: 20.47 fb^{-1}(8TeV)"));
	c->Print(TString::Format("./plots/%s_afb.pdf",outfile));
}//}}}
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
       yuerrafb2[i] = yderrafb2[i];
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
	
  RooRealVar nsig("nsig","nsig",2E4,0,4E6);
	RooRealVar nbkgComb("nbkgComb","nbkgComb",1E5,0,6E6);
	
	RooAddPdf f("kernel","kernel",RooArgList(f_bkgComb,f_sig),RooArgList(nbkgComb,nsig));// no penalty term
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>> INFO: f_penalty NOT prepared. <<<<<<<<<<<<<<<<<<<<"<<endl;
///////////////////////////////////////////////////////////// p.d.f. ///////////////////////////////////////////////////	

	// Get data and apply unbinned fit
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
	RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Bmass,CosThetaL,Q2),Q2range[iBin],0);
	//RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE), Minimizer("Minuit"), Minos(RooArgSet(afb, fh)), Warnings(1), PrintEvalErrors(3), Verbose(1));	
	//cout<<"BUGGGGGGGG"<<endl;
   f_fitresult->Print();
   delete data;
   odatacardpath=summarypath;
   double val[2]={0,0};
  val[0] = f_fitresult->status();
  writeParam(iBin,  "migrad",  val, 1);
  val[0] = f_fitresult->minNll();
  writeParam(iBin, "FC_FCN", val);
  val[0] = 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi();
  writeParam(iBin,  "FC_fh",  val, 2);
  val[0] = (1. * atan( afb.getVal() ) / TMath::Pi()) * ( 3./2. + 3. * atan( fh.getVal()  ) / TMath::Pi() ); 
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

    RooDataSet *outdata = f.generate(RooArgSet(*Bmass,*CosThetaL),nEvents);
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

// void rndPickMCSamples(int nSets = 1, double oLumis=datasetLumi[0]) // Pick distinct events from MC
// {//{{{
//     static char decmode[20];
//     while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s")*strcmp(decmode, "signal")*strcmp(decmode, "manual") != 0){
//         printf("Please insert type [ signal / jpsi / psi2s / manual ]:");
//         scanf("%19s",decmode);
//     }
//     
//     double iLumis = 0.;
//     if (strcmp(decmode, "signal") == 0){
//         iLumis=datasetLumi[1];
//     }else if (strcmp(decmode, "jpsi") == 0){
//         iLumis = datasetLumi[2];
//     }else if (strcmp(decmode, "psi2s") == 0){
//         iLumis = datasetLumi[3];
//     }else{
//         printf("Input luminosity: ");
//         scanf("%lf",&iLumis);
//     }
//     printf("INFO\t\t: Randomly pick %d sets of events of %f/fb \n", nSets, oLumis);
// 
//     if (oLumis > iLumis){
//         printf("ERROR\t\t: Luminosity of output should be larger than the reservoir. EXIT 1.");
//         return;
//     }
// 
//     // Initialize random generator
//     TRandom3 *random = new TRandom3(time(0));
// 
//     int nEntries = ch->GetEntries();
//     int nEvtsPerSet = floor(oLumis/iLumis*nEntries);
//     //int mumuMassWindowBin = 0;
//     for(int iSet=0; iSet< nSets; iSet++){
//         TFile *fout = new TFile(TString::Format("%s/rndMC_%s_%devts_set%d.root",owspacepath.Data(),decmode,nEvtsPerSet,iSet+1), "RECREATE");
//         TTree *tout = ch->CloneTree(0);
//         int candEvtId=0;
//         std::vector<bool> bits(nEntries,false); // all false;
//         for(int iEvt = nEvtsPerSet; iEvt > 0; iEvt--){
//             do {
//                 printf("DEBUG\t\t: iEvt=%d, candEvtId=%d, Picking distict.\n",iEvt,candEvtId);
//                 candEvtId = random->Rndm()*nEntries;
//             }while(bits[candEvtId]);
//             bits[candEvtId]=true;
//             // Don't access files here, random access really takes time!
//         }
//         for(int iBit = 0; iBit < nEntries; iBit++){
//             if (bits[iBit] == false) continue;
//             ch->GetEntry(candEvtId);
//             tout->Fill();
//         }
//         //tout = tout->CopyTree("mumuMassWindow[mumuMassWindowBin]");
//         tout->FlushBaskets();
//         tout->AutoSave();
//         fout->Write();
//         fout->Close();
//     }
// }//}}}

void createFCToys(int iBin, int nToy=500)// Create toys for contour scanning
{//{{{
    // Remark: keep extension possibility for 2-D contour.
    // Buffers for checking directory/FILE
    //bool contourMode=false;// In case contour is needed.
    struct stat fiBuff;
    //FILE *fi = 0;

    // Setting
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    //const double dataAfb[11]   = {-0.2104, -0.0443, 0.0024, 0.00, 0.0019, 0.00, 0.0285, 0.0347, 0.0522,-0.1444, 0.0102};
    //const double dataFh[11]    = { 0.4885,  0.8491, 0.0694, 0.00, 0.0080, 0.00, 0.0571, 0.0736, 0.1088, 0.3447, 0.0305};
    const double stepSizeAfb[11] = { 0.0040,  0.0060, 0.0005, 0.00, 0.0001, 0.00, 0.0005, 0.0005, 0.0010, 0.0020, 0.0005};
    const double stepSizeFh[11]  = { 0.0040,  0.0080, 0.0010, 0.00, 0.0002, 0.00, 0.0005, 0.0005, 0.0010, 0.0040, 0.0005};
    double dataAfbLower[11], dataFhLower[11], dataFhBand[11];

    // Create output directory
    if (stat(TString::Format("./limit"),&fiBuff) != 0){
        mkdir(TString::Format("./limit"),0755);
        if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
            mkdir(TString::Format("%s",otoyspath.Data()),0755);
        }
    }

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh   = wspace->var("fh")->getVal();
    double  afb  = wspace->var("afb")->getVal();
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

    // Loop over phase space
    dataAfbLower[iBin] = -0.5*fh;
    double thisAfb = afb;
    double thisFh = fh;
    for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
        thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
        }
        owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
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
    }// Afb loop
    dataFhLower[iBin]  = fabs(afb)*2.0;
    dataFhBand[iBin]   = fh*2.0;
    thisAfb = afb;
    thisFh = fh;
    //for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
    for(int iFh = 0; (((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < 1.) or ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin])); iFh++){
        if ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]) {
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        } else {
           if (iBin == 2 || iBin == 10 || iBin == 4 || iBin == 6 or iBin == 7) iFh = iFh + 49;
           else iFh = iFh + 9;
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        }
        if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
        }
        owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
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
//    // 2-D loop for errFh and contour
//    thisAfb = afb;
//    thisFh = fh;
//    for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
//        thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
//        for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
//            thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
//            //cout<<thisAfb<<"   "<<thisFh<<endl;
//            if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
//            if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
//                mkdir(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
//            }
//            owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
//
//            genToySignal(iBin, nsigInToys,thisAfb,thisFh);
//            splitToySamples(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
//            genToyCombBkg(iBin,nbkgCombInToys);
//            splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
//            for(int iToy = 0; iToy<nToy; iToy++){
//                if (stat(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),&fiBuff) != 0){
//                    mkdir(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),0755);
//                }
//                rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
//                       TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
//                rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,owspacepath.Data(),       iBin,iToy+1).Data(),
//                       TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
//            }
//        }// Afb loop
//        //if (!contourMode) break;
//    }// Fh loop
    return;
}//}}}

void harvestFCFitResults(int iBin, int nToy=500, bool contourMode=false)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    contourMode=false;// In case contour is needed.
    // Setting
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const int    nHistBins   = 2000;
    //const double dataAfb[11]   = {-0.2104, -0.0443, 0.0024, 0.00, 0.0019, 0.00, 0.0285, 0.0347, 0.0522,-0.1444, 0.0102};
    //const double dataFh[11]    = { 0.4885,  0.8491, 0.0694, 0.00, 0.0080, 0.00, 0.0571, 0.0736, 0.1088, 0.3447, 0.0305};
    const double stepSizeAfb[11] = { 0.0040,  0.0060, 0.0005, 0.00, 0.0001, 0.00, 0.0005, 0.0005, 0.0010, 0.0020, 0.0005};
    const double stepSizeFh[11]  = { 0.0040,  0.0080, 0.0010, 0.00, 0.0002, 0.00, 0.0005, 0.0005, 0.0010, 0.0040, 0.0005};
    double dataAfbLower[11], dataFhLower[11], dataFhBand[11];

    // Get parameters for the q2 Bin
    //iwspacepath="./RootFiles";
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    //double  fh  = readParam("fh",TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal();
    //double  afb = readParam("afb",TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal();
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
    printf("INFO\t\t: bounded fh=%+.3f, bounded afb=%+.3f\n",fh,afb);

    // Loop over phase space
    dataAfbLower[iBin] = -0.5*fh;
    double thisAfb = afb;
    double thisFh = fh;
    for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
        thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;

        TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
        TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
        TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);

        for(int iToy = 0; iToy<nToy; iToy++){
            //iwspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            idatacardpath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            if (readParam(iBin,"migrad", 0) == 0.){
                if (readParam(iBin,"FC_fh",0) > 3e-3) { 
                h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));
                h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0));
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
    dataFhLower[iBin]  = fabs(afb)*2.0;
    dataFhBand[iBin]   = fh*2.0;
    thisAfb = afb;
    thisFh = fh;
    //for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
    for(int iFh = 0; (((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < 1.) or ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin])); iFh++){
        if ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]) {
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        } else {
           if (iBin == 2 || iBin == 10 || iBin == 4 || iBin == 6 or iBin == 7) iFh = iFh + 49;
           else iFh = iFh + 9;
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        }
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;

        TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
        TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
        TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);

        for(int iToy = 0; iToy<nToy; iToy++){
            idatacardpath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
            if (readParam(iBin,"migrad", 0) == 0.){
                if (readParam(iBin,"FC_fh",0) > 3e-3) { 
                h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));
                h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0));
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
//    // 2-D loop for errFh and contour
//    thisAfb = afb;
//    thisFh = fh;
//    for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
//        thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
//        for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
//            thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
//            //if ( (thisAfb<-0.006 && thisFh==0.032) || (thisFh>0.032) ) continue;
//            if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
//            owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
//            if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;
//
//            TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
//            TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
//            TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);
//
//            for(int iToy = 0; iToy<nToy; iToy++){
//                //iwspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
//                idatacardpath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
//                if (readParam(iBin,"migrad", 0) == 0.){
//                    h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));
//                    h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0));
//                }
//            }
//            h_setSummaryAfb->Draw();
//            h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb.cc",owspacepath.Data()));
//            h_setSummaryFh->Draw();
//            h_setSummaryFh->SaveAs(TString::Format("%s/h_setSummaryFh.cc",owspacepath.Data()));
//            fout->Write();
//            fout->Close();
//        }// Afb loop
//        //if (!contourMode) break;
//    }// Fh loop
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
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Settings
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    //const double dataAfb[11]   = {-0.2104, -0.0443, 0.0024, 0.00, 0.0019, 0.00, 0.0285, 0.0347, 0.0522,-0.1444, 0.0102};
    //const double dataFh[11]    = { 0.4885,  0.8491, 0.0694, 0.00, 0.0080, 0.00, 0.0571, 0.0736, 0.1088, 0.3447, 0.0305};
    const double stepSizeAfb[11] = { 0.0040,  0.0060, 0.0005, 0.00, 0.0001, 0.00, 0.0005, 0.0005, 0.0010, 0.0020, 0.0005};
    const double stepSizeFh[11]  = { 0.0040,  0.0080, 0.0010, 0.00, 0.0002, 0.00, 0.0005, 0.0005, 0.0010, 0.0040, 0.0005};
    double dataAfbLower[11], dataFhLower[11], dataFhBand[11];

    // Get parameters for the q2 Bin
    //iwspacepath="./RootFiles";
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
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
    dataAfbLower[iBin] = -0.5*fh;
    dataFhLower[iBin]  = fabs(afb)*2.0;
    dataFhBand[iBin]   = fh*2.0;
    double thisAfb = afb;
    double thisFh = fh;
    TFile *fin = 0;
    TFile *fin1 = 0;
    TH1F  *h_setSummaryAfb = 0;
    TH1F  *h_setSummaryFh  = 0;
    // loop for errAfb
    for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
        thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        fin = new TFile(TString::Format("%s/afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
        h_setSummaryAfb = (TH1F*)fin->Get("h_setSummaryAfb");
        getCIFromTH1F(h_setSummaryAfb,afbMeasuredLo[nAfbPoints],afbMeasuredHi[nAfbPoints]);
        afbTruth[nAfbPoints] = thisAfb;
        nAfbPoints++;
        fin->Close();
    }
    thisAfb = afb;
    thisFh = fh;
    //for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
    for(int iFh = 0; (((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < 1.) or ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin])); iFh++){
        if ((dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]) {
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        } else {
           if (iBin == 2 || iBin == 10 || iBin == 4 || iBin == 6 or iBin == 7) iFh = iFh + 49;
           else iFh = iFh + 9;
           thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
        }
        if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
        fin1 = new TFile(TString::Format("%s/afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
        h_setSummaryFh  = (TH1F*)fin1->Get("h_setSummaryFh");
        getCIFromTH1F(h_setSummaryFh,fhMeasuredLo[nFhPoints],fhMeasuredHi[nFhPoints]);
        fhTruth[nFhPoints] = thisFh;
        nFhPoints++;
        fin1->Close();
    }

//    for(int iFh = 0; (dataFhLower[iBin]+iFh*stepSizeFh[iBin]) < dataFhBand[iBin]; iFh++){
//        thisFh = dataFhLower[iBin]+iFh*stepSizeFh[iBin];
//        for(int iAfb = 0; (dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin]) < (0.5*thisFh); iAfb++){
//            thisAfb = dataAfbLower[iBin]+iAfb*stepSizeAfb[iBin];
//            //if ( fabs(thisAfb)==(0.5*thisFh) ) continue;
//            if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
//            fin = new TFile(TString::Format("%s/afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
//            h_setSummaryAfb = (TH1F*)fin->Get("h_setSummaryAfb");
//            getCIFromTH1F(h_setSummaryAfb,afbMeasuredLo[nAfbPoints],afbMeasuredHi[nAfbPoints]);
//            afbTruth[nAfbPoints] = thisAfb;
//            nAfbPoints++;
//            fin->Close();
//            
//            if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
//            fin1 = new TFile(TString::Format("%s/afb%+05.0f_fh%+05.0f/setSummary.root",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data());
//            h_setSummaryFh  = (TH1F*)fin1->Get("h_setSummaryFh");
//            getCIFromTH1F(h_setSummaryFh,fhMeasuredLo[nFhPoints],fhMeasuredHi[nFhPoints]);
//            fhTruth[nFhPoints] = thisFh;
//            nFhPoints++;
//            fin1->Close();
//        }
//    }

    // Make plots, fit, and find confidence interval
    TFile *fout = new TFile(TString::Format("%s/wspace_FCConfInterval_bin%d.root",owspacepath.Data(),iBin),"RECREATE");
    TCanvas *canvas = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLine *line = new TLine();
    TGraph *g_fhIntervalHi = new TGraph(nFhPoints,fhTruth, fhMeasuredHi);
    TGraph *g_fhIntervalLo = new TGraph(nFhPoints,fhTruth, fhMeasuredLo);
    TGraph *g_afbIntervalHi = new TGraph(nAfbPoints,afbTruth, afbMeasuredHi);
    TGraph *g_afbIntervalLo = new TGraph(nAfbPoints,afbTruth, afbMeasuredLo);
//    TF1 *f1_fhIntervalHi  = new TF1("f1_fhIntervalHi",  "[0]+[1]*x+[2]*x**2",0,3);
//    TF1 *f1_fhIntervalLo  = new TF1("f1_fhIntervalLo",  "[0]+[1]*x+[2]*x**2",0,3);
//    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x+[2]*x**2",-1,1);
//    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x+[2]*x**2",-1,1);
//    TFitResultPtr r_fhIntervalHi = g_fhIntervalHi->Fit("f1_fhIntervalHi","S","",max(0.01,fh-0.1),min(0.99,fh+0.3));
//    TFitResultPtr r_fhIntervalLo = g_fhIntervalLo->Fit("f1_fhIntervalLo","S","",max(0.01,fh-0.3),min(0.99,fh+0.1));
//    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",max(0.01,afb-0.1),min(0.99,afb+0.3));
//    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",max(0.01,afb-0.3),min(0.99,afb+0.1));
    TF1 *f1_fhIntervalHi  = new TF1("f1_fhIntervalHi",  "[0]+[1]*x+[2]*x**2",dataFhLower[iBin]*0.9,max(1., dataFhBand[iBin]));
    TF1 *f1_fhIntervalLo  = new TF1("f1_fhIntervalLo",  "[0]+[1]*x+[2]*x**2",dataFhLower[iBin]*0.9,max(1., dataFhBand[iBin]));
    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x",-0.55*fh,0.55*fh);
    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x",-0.55*fh,0.55*fh);
    TFitResultPtr r_fhIntervalHi = g_fhIntervalHi->Fit("f1_fhIntervalHi","S","",dataFhLower[iBin]*0.9,max(1., dataFhBand[iBin]));
    TFitResultPtr r_fhIntervalLo = g_fhIntervalLo->Fit("f1_fhIntervalLo","S","",dataFhLower[iBin]*0.9,max(1., dataFhBand[iBin]));
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",-0.55*fh,0.55*fh);
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",-0.55*fh,0.55*fh);
    if (iBin == 1) {
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",-0.2,0.2);
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",-0.2,0.2);
    }
    // Calculate the error, you must think about the case in which f1 not defined!
    double XPFh[2], XPAfb[2];
    if (fabs(f1_afbIntervalLo->GetX(afb)) <= 0.5*fh) XPAfb[0]=f1_afbIntervalLo->GetX(afb);
    else XPAfb[0]=-0.5*fh;
    if (fabs(f1_afbIntervalHi->GetX(afb)) <= 0.5*fh) XPAfb[1]=f1_afbIntervalHi->GetX(afb);
    else XPAfb[1]= 0.5*fh;
    if (f1_fhIntervalHi->GetX(fh) > dataFhLower[iBin] && f1_fhIntervalHi->GetX(fh) < fh) XPFh[0]=f1_fhIntervalHi->GetX(fh);
    else XPFh[0]=dataFhLower[iBin];
    if (f1_fhIntervalLo->GetX(fh) > dataFhLower[iBin] && f1_fhIntervalLo->GetX(fh) < 3.) XPFh[1]=f1_fhIntervalLo->GetX(fh);
    else XPFh[1]=3.;
    outFCErrFh[0] = r_fhIntervalHi.Get()  != 0 ? XPFh[0]-fh : 0.-fh;
    outFCErrFh[1] = r_fhIntervalLo.Get()  != 0 ? XPFh[1]-fh : 3.-fh;
    outFCErrAfb[0]= r_afbIntervalLo.Get() != 0 ? XPAfb[0]-afb : -0.5*fh-afb;
    outFCErrAfb[1]= r_afbIntervalHi.Get() != 0 ? XPAfb[1]-afb : 0.5*fh-afb;
//    outFCErrFh[0] = r_fhIntervalHi.Get()  != 0 ? max(dataFhLower[iBin],f1_fhIntervalHi->GetX(fh))-fh : 0.-fh;
//    outFCErrFh[1] = r_fhIntervalLo.Get()  != 0 ? min(3.,f1_fhIntervalLo->GetX(fh))-fh : 3.-fh;
//    outFCErrAfb[0]= r_afbIntervalLo.Get() != 0 ? max(-0.5*fh,f1_afbIntervalLo->GetX(afb))-afb : -0.5*fh-afb;
//    outFCErrAfb[1]= r_afbIntervalHi.Get() != 0 ? min(0.5*fh,f1_afbIntervalHi->GetX(afb))-afb : 0.5*fh-afb;
//    if ( outFCErrFh [1] < 0 ){
//        outFCErrFh [1] = 0;
//    }
//    if ( outFCErrAfb[0] > 0 ){
//        outFCErrAfb[0] = 0;
//        outFCErrFh [1] = 0;
//    }
//    if ( outFCErrAfb[1] < 0 ){
//        outFCErrAfb[1] = 0;
//        outFCErrFh [1] = 0;
//    }

    const double upAfb[11] = { 0.40,  0.99, 0.10, 0.00, 0.10, 0.00, 0.10, 0.10, 0.10, 0.50, 0.10};
    const double upFh[11]  = { 1.00,  1.50, 0.20, 0.00, 0.20, 0.00, 0.40, 0.40, 0.50, 1.00, 0.40};
    // Fh plots
    g_fhIntervalLo->SetTitle("");
    g_fhIntervalLo->GetYaxis()->SetTitle("Measured F_{H}");
    //g_fhIntervalLo->GetYaxis()->SetLimits(0.,1.1*max(fh,fhTruth[nFhPoints-1]));
    g_fhIntervalLo->GetYaxis()->SetRangeUser(0.,upFh[iBin]);
    g_fhIntervalLo->GetXaxis()->SetTitle("True F_{H}");
    //g_fhIntervalLo->GetXaxis()->SetRangeUser(0.,1.1*fhTruth[nFhPoints-1]);
    //g_fhIntervalHi->GetXaxis()->SetRangeUser(0.,1.1*fhTruth[nFhPoints-1]);
    //g_fhIntervalHi->GetXaxis()->SetRangeUser(*min_element(fhTruth, fhTruth+maxPoints), *max_element(fhTruth, fhTruth+maxPoints));
    g_fhIntervalLo->Draw("AP");
    g_fhIntervalHi->Draw("P SAME");
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(0.,fh,fhTruth[nFhPoints-1]*1.1,fh);
    //line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(fabs(afb*2.),0,fabs(afb*2.),upFh[iBin]);
    //line->DrawLine(g_fhIntervalHi->GetXaxis()->GetXmin(),0.5*fhTruth[nFhPoints-1],g_fhIntervalHi->GetXaxis()->GetXmax(),0.5*fhTruth[nFhPoints-1]);
    latex->DrawLatexNDC(0.15,0.84,TString::Format("F_{H}=%.3f^{%+.3f}_{%+.3f}",fh,outFCErrFh[1],outFCErrFh[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_fh_bin%d.pdf",plotpath.Data(),iBin));
    // Afb plots
    g_afbIntervalHi->SetTitle("");
    g_afbIntervalHi->GetYaxis()->SetTitle("Measured A_{FB}");
    //g_afbIntervalHi->GetYaxis()->SetLimits(-1.1*min(afb,afbTruth[nAfbPoints-1]),1.1*max(afb,afbTruth[nAfbPoints-1]));
    //g_afbIntervalHi->GetYaxis()->SetLimits(-0.5,0.5);
    //g_afbIntervalHi->GetYaxis()->SetLimits(-upAfb[iBin],upAfb[iBin]);
    g_afbIntervalHi->GetYaxis()->SetRangeUser(-upAfb[iBin],upAfb[iBin]);
    g_afbIntervalHi->GetXaxis()->SetTitle("True A_{FB}");
    g_afbIntervalHi->GetXaxis()->SetRangeUser(-1.1*afbTruth[nAfbPoints-1],1.1*afbTruth[nAfbPoints-1]);
    g_afbIntervalHi->Draw("AP");
    g_afbIntervalLo->Draw("P SAME");
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(-1*afbTruth[nAfbPoints-1]*1.1,afb,afbTruth[nAfbPoints-1]*1.1,afb);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(afbTruth[nAfbPoints-1],-upAfb[iBin],afbTruth[nAfbPoints-1],upAfb[iBin]);
    line->DrawLine(-1*afbTruth[nAfbPoints-1],-upAfb[iBin],-1*afbTruth[nAfbPoints-1],upAfb[iBin]);
    latex->DrawLatexNDC(0.15,0.84,TString::Format("A_{FB}=%.3f^{%+.3f}_{%+.3f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
    latex->DrawLatexNDC(0.79,0.91, Form("bin %d", iBin));
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_afb_bin%d.pdf",plotpath.Data(),iBin));

    fout->Close();

    writeParam(iBin,"FCErrFh" ,outFCErrFh );
    writeParam(iBin,"FCErrAfb",outFCErrAfb);

    return;

}//}}}
void getFCInterval2D(int iBin)
{//{{{
}//}}}

    // Sys. Error determination tools
bool isEffMapNonNegative(RooGenericPdf *f_sigA)
{//{{{
    return true;
    //return false;
}//}}}
void rndEfficiencyMap(int iBin, RooGenericPdf *f_sigA, TMatrixDSym *errMtx)
{//{{{

    //const char varNames[][32]={};
    // Prepare multivariate Gaussian
    //RooArgSet *f_sigA_arg = f_sigA->getVariables();

    // randomly generate maps
    do{
    }while(!isEffMapNonNegative(f_sigA));

    printf("INFO: f_sigA transformed randomly.\n");
    return;
}//}}}
void trimErrMtx(TMatrixDSym *iErrMtx, TMatrixDSym *oErrMtx)
{//{{{
    return;
}//}}}


void createValiToys(int iBin, int nToy=100)// Create toys for contour scanning
{//{{{
    struct stat fiBuff;

    // Setting
    TString otoyspath = TString::Format("./Toy/bin%d",iBin);
    // Create output directory
    if (stat(TString::Format("./Toy"),&fiBuff) != 0){
        mkdir(TString::Format("./Toy"),0755);
        if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
            mkdir(TString::Format("%s",otoyspath.Data()),0755);
        }
    }

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh   = wspace->var("fh")->getVal();
    double  afb  = wspace->var("afb")->getVal();
    double  nsig = wspace->var("nsig")->getVal();
    double  nbkgComb = wspace->var("nbkgComb")->getVal();

    TRandom3 *rndGenerator = new TRandom3();
    int nsigInToy[nToy];
    int nbkgCombInToy[nToy];
    int nsigInToys = 0;
    int nbkgCombInToys = 0;
    for(int iToy=0; iToy<nToy; iToy++){
        nsigInToy[iToy]=rndGenerator->Poisson(nsig) * 50;
        nbkgCombInToy[iToy]=rndGenerator->Poisson(nbkgComb) * 50;
        nsigInToys+=nsigInToy[iToy];
        nbkgCombInToys+=nbkgCombInToy[iToy];
    }

    // Loop over phase space
    double thisAfb = afb;
    double thisFh = fh;
    //if (!scanAfbFhPositivePdf(thisAfb,thisFh,true)) continue;
    if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0){
        mkdir(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),0755);
    }
    owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
    //if (stat(TString::Format("%s/set0001",owspacepath.Data()),&fiBuff) == 0) { cout<<thisAfb<<"  "<<thisFh<<endl; continue;}

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
    return;
}//}}}

void harvestValiFitResults(int iBin, int nToy=500, bool contourMode=false)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    // Setting
    TString otoyspath = TString::Format("./Toy/bin%d",iBin);
    const int    nHistBins   = 2000;

    // Get parameters for the q2 Bin
    //iwspacepath="./RootFiles";
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fh  = wspace->var("fh")->getVal();
    double  afb = wspace->var("afb")->getVal(); 
    printf("INFO\t\t: bounded fh=%+.3f, bounded afb=%+.3f\n",fh,afb);

    // Loop over phase space
    double thisAfb = afb;
    double thisFh = fh;
    //if (stat(TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000).Data(),&fiBuff) != 0) continue;
    owspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f",otoyspath.Data(),thisAfb*10000,thisFh*10000);
    //if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;

    TFile *fout = new TFile(TString::Format("%s/setSummary.root",owspacepath.Data()),"RECREATE");
    TH1F *h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb%+05.0f",thisAfb*10000).Data(), nHistBins, -1.,1.);
    TH1F *h_setSummaryFh  = new TH1F("h_setSummaryFh", TString::Format("h_fh%+05.0f",thisFh*10000).Data(), nHistBins, 0.,3.);

    for(int iToy = 0; iToy<nToy; iToy++){
        //iwspacepath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
        idatacardpath=TString::Format("%s/afb%+05.0f_fh%+05.0f/set%04d",otoyspath.Data(),thisAfb*10000,thisFh*10000,iToy+1);
        if (readParam(iBin,"migrad", 0) == 0.){
            if (readParam(iBin,"FC_fh",0) > 3e-3) { 
            h_setSummaryFh  ->Fill(readParam(iBin,"FC_fh",0));
            h_setSummaryAfb ->Fill(readParam(iBin,"FC_afb",0));
            }
        }
    }
    h_setSummaryAfb->Draw();
    h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb.cc",owspacepath.Data()));
    h_setSummaryFh->Draw();
    h_setSummaryFh->SaveAs(TString::Format("%s/h_setSummaryFh.cc",owspacepath.Data()));
    fout->Write();
    fout->Close();
    return;

}//}}}



// Other tools

//_________________________________________________________________________________

void printListOfTChainElements(TChain *chain){
    TObjArray *fileElements=chain->GetListOfFiles();
    int nFiles = fileElements->GetEntries();
    //printf("DEBUG\t\t: %d files in the chain\n",nFiles);
    TIter next(fileElements);
    TChainElement *chEl=0;
    for( int entry=0; entry < nFiles; entry++ ) {
        chEl=(TChainElement*)next();
        printf("%s\n",chEl->GetTitle());
    }
    printf("DEBUG\t\t: %d files in the chain\n",nFiles);
}



std::string GetDirectory (const std::string& path)
{
   size_t found = path.find_last_of("/");
   return(path.substr(0, found));
}





//_________________________________________________________________________________
//_________________________________________________________________________________
int main(int argc, char** argv) {
//	Tags
	is7TeVCheck = false;   
//	Help message
	if (argc <= 2) {
		cout<<"Usage       : ./fit Function infile binID"<<endl;
		cout<<"Functions   :"<<endl;
		cout<<"    0. bmass               Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg."<<endl;
		cout<<"    5. angular2D_1a_Sm     Leading step1 to angular2D, determine signal shape from simulation."<<endl;
		cout<<"    8. angular2D_prior     Leading step4 to angular2D, fit to data sideband to get initial values of combinatorial bkg."<<endl;
		cout<<"    9. angular2D           Derive F_{H} and A_{FB} by fitting to mass and angular distribution."<<endl;
		cout<<"       For data fitting test:        ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
		cout<<"       For initial values scanning:          ./fit <function> <input.root> <iBin>  "<<endl;
		cout<<"   10. PlotFCN             Plot FCN distribution for each q^{2} bin, find the samllest FCN value"<<endl;
		cout<<"                              and it's initial value, save the fitting results."<<endl;
		cout<<"       For refit with final initial values:  ./fit <function> <input.root> <iBin> refit  "<<endl;
		cout<<"   11. PlotAfbFh           Plot Afb-Fh 2D distribution for each q^{2} bin, and the final results"<<endl;
		cout<<"Remark      :"<<endl;
		cout<<"    1. Outputs will be stored in ./plots, please keep the directory."<<endl;
		cout<<"    2. Fitted parameters will be stored in ./fitParameters/*.txt, please keep the directory."<<endl;
		cout<<"    3. Scaned fitted parameters will be stored in ./OutputValues/bin*/*.txt, please keep these directories."<<endl;
	//	cout<<"    3. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!"<<endl;
		return 0;
	}
//	main
	if (argc != 4 && argc != 7 && argc != 5 && argc != 17 && argc != 14 && argc != 15){
		cout<<"./fit func infile binID"<<endl;
		for (int i = 0; i < 11; i++) {
		//	if (i == 3 || i == 5) continue;
			printf("    Bin %d : %s\n",i,Q2range[i]);
		}
		return 0;
	}
  int nToy=500;
	TString func    = argv[1];
	TString infile  = argv[2];
	int iBin        = atoi(argv[3]);
	
	if (func == "angular2D_Toy") {
		ch->Add(infile.Data());
      summarypath=GetDirectory(argv[2]);
      cout<<summarypath<<endl;
		if (ch == NULL) gSystem->Exit(0);
		angular2D_Toy_bin(iBin); 
	} else if (func == "angular2D_Toy_unCons") {
		ch->Add(infile.Data());
      summarypath=GetDirectory(argv[2]);
      cout<<summarypath<<endl;
		if (ch == NULL) gSystem->Exit(0);
		angular2D_Toy_unCons_bin(iBin); 
	} else if (func == "bmass") {
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		const char outfile[]="bmass";
		bmass(iBin, outfile); 
	}else if (func == "angular2D_1a_Sm" || func == "angular2D_prior" ){
		void (*fx)(int, const char*, bool);
		if ( func == "angular2D_1a_Sm" ){
			fx = angular2D_1a_Sm;
		}else if( func == "angular2D_prior" ){
			fx = angular2D_prior;
		} else {
			return 0;
		}	
		ch->Add(infile.Data());
		if (ch == NULL) gSystem->Exit(0);
		fx(iBin,func,true);// By default overwrite exist parameters.
   }else if (func == "angular2D" || func == "angular_reco"){
		if (iBin >= 0 && iBin < 11) {
			ch->Add(infile.Data());
			if (ch == NULL) gSystem->Exit(0);
			std::vector<double> vbin;
			
			double Iafb, Ifh, iafb, ifh;
			int Index = -999;
			if (argc == 7) {   //////////////////////////////   scanning   ///////////////////////////////////////////
				Index  = atoi(argv[4]);
				Iafb   = atof(argv[5]);
				Ifh    = atof(argv[6]);
				cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>>>> Index = "<<Index<<" >>>>>> Iafb = "<<Iafb<<" >>>>>> Ifh = "<<Ifh<<endl;
				if (func == "angular2D") {
					vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
				} else if (func == "angular_reco") {
					vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
				} else {
				//	cout<<"A test of data fitting with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> <afb> <fh> test "<<endl;
					return 0;
				}
			} else if (argc == 5) {  //////////////////////////////   refit ////////////////////////////////////
				TString type  = argv[4];
				if ( type == "refit") {
					Index = -1;
					if (func == "angular2D") {
						const char outfile[] = "angular2D";
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);  // Selected initial values!
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular2D_bin(iBin, Iafb, Ifh, Index);
            //vbin = angular2D_bin(iBin, iafb, ifh, Index);   // unconstrained fitting
					} else if (func == "angular_reco") {
						const char outfile[] = "angular_reco";
						iafb = readParam(iBin, TString::Format("Iafb_%s",outfile), 0);
						ifh  = readParam(iBin, TString::Format("Ifh_%s",outfile), 0);
						Iafb = (1. * atan( iafb ) / TMath::Pi()) * ( 3./2. + 3. * atan( ifh  ) / TMath::Pi() );
						Ifh  = 3./2. + 3. * atan( ifh  ) / TMath::Pi();
						vbin = angular_reco_bin(iBin, Iafb, Ifh, Index);
				   } else {
					   return 0;
					}
					cout<<endl<<">>>>>>>>>> iBin = "<<iBin<<" >>>>>> Iafb_"<<func<<" = "<<Iafb<<" >>>>>> Ifh_"<<func<<" = "<<Ifh<<endl;
				} else {
					cout<<"Refit with initial values: "<<endl<<" ./fit <function> <input.root> <iBin> refit "<<endl;
					return 0;
				}
			}
		}else if (iBin == 999) {
			if (func == "angular2D") {
				const char outfile[] = "angular2D";
				angular2D(outfile);
			} else if (func == "angular_reco") {
				const char outfile[] = "angular_reco";
				angular_reco(outfile);
			}
		}else if (iBin == 998) {
			if (func == "angular2D") {
				const char outfile[] = "angular2D";
				angular2D_vali(outfile);
			}
		}else { 
			cout<<"Refit data, iBin counts from 0 to 10; or 999 to plot the results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
	} else if (func == "PlotFCN") {
		TString label = argv[4];
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			PlotFCN(iBin, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			PlotFCN(iBin, outfile);
		}
	//	const char outfile[]="FCN";
	} else if (func == "PlotAfbFh_f") {	
		TString label = argv[4];
		int Index;
		if (iBin >= 0 && iBin < 11) {
			Index = -999;
		}else if (iBin == 999) {
			Index = -1;
		}else if (iBin == 998) {
			Index = -2;
		}else { 
			cout<<"Plot Afb and Fh 2D distribution: for each q2 iBin counts from 0 to 10; or 999 to plot the refitted results!"<<endl;
			cout<<" Please check the Usage : ./fit"<<endl;
			return 0; 
		}
		if (label == "angular2D") {
			const char outfile[] = "angular2D";
			PlotAfbFh_f(iBin, Index, outfile);
		} else if (label == "angular_reco") {
			const char outfile[] = "angular_reco";
			PlotAfbFh_f(iBin, Index, outfile);
		}
   }else if (func == "createToys"){
      scaleFactor = 100;
      double expNCombBkg[11] = { 518.5 , 1125.1 , 1586.7 , 0 , 755.8 , 0 , 248.2 , 242.8 , 286.1 , 2355.5 , 4838.5 };
      double expNSig[11]     = { 173.5 ,  330.9 ,  783.3 , 0 , 360.2 , 0 , 214.8 , 262.1 , 225.9 ,  779.5 , 2275.4 };
      for (int iBin = 0; iBin < 11; iBin++) {
            genToySignal(iBin,expNSig[iBin]*scaleFactor);// (iBin, nEvts)
            genToyCombBkg(iBin,expNCombBkg[iBin]*scaleFactor,"combBkgToy");
      }
      ch->Add(infile.Data());
      if (ch == NULL) return 1;
      splitMCSamples();
      //rndPickMCSamples(scaleFactor);// (nSets,lumis)
   } else if (func == "ToyVali"){
        // Get parameters for the q2 Bin
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_pdf_bin%d.root",iwspacepath.Data(),iBin));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        if (!wspace) return 1;
        double  fh   = wspace->var("fh")->getVal();
        double  afb  = wspace->var("afb")->getVal();
        double  nsig = wspace->var("nsig")->getVal();
        double  nbkgComb = wspace->var("nbkgComb")->getVal();
        owspacepath=TString::Format("./Toy/bin%d",iBin);
        //genToySignal(iBin, nsig*100, 0.10,0.40); // v0
        //genToySignal(iBin, nsig*100, 0.01,0.05); // v1
        //genToySignal(iBin, nsig*100, -0.30,0.80); // v2
        genToySignal(iBin, nsig*100, afb,fh); // 
        genToyCombBkg(iBin,nbkgComb*100);
    } else if (func == "createValiToys"){
        nToy = 100;
        createValiToys(iBin,nToy);
    }else if (func == "ValiScan"){
        nToy = 100;
        int whichBin=-1;
        whichBin=iBin;
        harvestValiFitResults(whichBin,nToy);
    } else if (func == "createFCToys"){
      createFCToys(iBin,nToy);
    }else if (func == "FCScan"){
        int whichBin=-1;
        whichBin=iBin;
        harvestFCFitResults(whichBin,nToy);
        cout<<"----------------"<<endl;
        getFCInterval(whichBin,nToy);
    }
/////////////////////////////////////////////////////////////////////////////////////////////
	 else if (func == "test"){
        ch->Add(infile.Data());
        if (ch == NULL) gSystem->Exit(0);
        for (int iBin = 0; iBin < 8; iBin++) {
        }
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    printf("%lld entries processed.\n",ch->GetEntries());
    gSystem->Exit(0);

    return 0 ;
}
