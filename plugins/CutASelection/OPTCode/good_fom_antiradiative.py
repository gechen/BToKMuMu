#!/bin/env python

from ROOT import TFile, TCut, TGraphErrors, TTree, TCanvas, gDirectory
from math import *
from array import array


allcut = { "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.9950,0.9999, 0.0001], "TrkPt":["Tk_4vec.Pt() >", 1.4,0.6,3.0,0.1], "Trkdcasigbs":["Trkdcasigbs >",3.3, 1.0, 4.5, 0.1]}  # OPT

#allcut = { "Bvtxcl":["Bvtxcl >",0.12, 0.00, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.9950,0.9999, 0.0001], "Blxysig":["Blxysig >", 10.6, 3.0, 12.0, 0.1], "TrkPt":["Tk_4vec.Pt() >", 1.3,0.6,3.0,0.1], "Trkdcasigbs":["Trkdcasigbs >",3.3, 1.0, 4.5, 0.1]}  # OPT
#allcut = {  "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 10.7, 7.0, 12.0, 0.1], "TrkPt":["Tk_4vec.Pt() >", 1.5,0.6,3.0,0.1], "Trkdcasigbs":["Trkdcasigbs >",3.8, 2.5, 4.5, 0.1]}  # OPT
#allcut = {"Bvtxcl":["Bvtxcl >",0.08, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 7.0, 5.0, 10.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.3,0.6,2.0,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.7, 1.8, 3.5, 0.1]}  # OPT
#allcut = { "antirad_jpsi":["abs(Bmass - Mumumass - 2.182) >", 0.125, 0.05, 0.25,0.005],"antirad_psip":["abs(Bmass - Mumumass - 1.593) >", 0.07, 0.02, 0.12,0.005]  }  # veto 

#allcut = { "jpsimass_cdf":["(Bmass < 5.23) && abs(Bmass - Mumumass - 2.182) >", 0.140, 0.05, 0.25,0.005],"psimass_cdf":["(Bmass < 5.23) && abs(Bmass - Mumumass - 1.593) >", 0.06, 0.02, 0.12,0.005]  }  # veto 
#allcut = { "antirad_lhcbjpsi":["(Bmass < 5.23) && (Mumumass < 3.15) && (3.097 - Mumumass) <", 0.10, 0.01, 0.30,0.005]  }  # LHCb jpsi
#allcut = { "antirad_lhcbpsip":["(Bmass < 5.23) && (Mumumass < 3.75) && (3.686 - Mumumass) <", 0.12, 0.01, 0.30,0.005]  }  # LHCb psi



#allcut = { "Bvtxcl":["Bvtxcl >",0.12, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9995, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 7.7, 6.0, 14.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.2,0.8,2.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.8, 0.5, 3.5, 0.1], "jpsimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 2.182) >", 0.17, 0.04, 0.25,0.01],"psimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 1.593) >", 0.05, 0.01, 0.20,0.01]  } # veto
#allcut = { "jpsimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 2.182) >", 0.170, 0.05, 0.25,0.005],"psimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 1.593) >", 0.075, 0.02, 0.12,0.005]  }  # veto 

#allcut = { "jpsimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 2.182) >", 0.138, 0.05, 0.25,0.001],"psimass":["(Bmass < 5.23) &&  abs(Bmass - Mumumass - 1.593) >", 0.056, 0.02, 0.12,0.001]  }  # veto 



#allcut = { "jpsimass_H":["(Bmass > 5.33) &&  abs(Bmass - Mumumass - 2.182) >", 0.09, 0.04, 0.25,0.005],"psimass_H":["(Bmass > 5.33) &&  abs(Bmass - Mumumass - 1.593) >", 0.06, 0.01, 0.20,0.005]  }  # veto 


#allcut = { "Bvtxcl":["Bvtxcl >",0.09, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9996, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 7.0, 6.0, 12.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.3,0.8,2.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.7, 0.5, 3.5, 0.1], "lhcbjpsi":["(Bmass < 5.23) && (Mumumass < 3.15) && (3.097 - Mumumass) <", 0.13, 0.01, 0.30,0.01]  }  # LHCb jpsi
#allcut = { "Bvtxcl":["Bvtxcl >",0.09, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9996, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 7.0, 6.0, 12.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.3,0.8,2.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.7, 0.5, 3.5, 0.1], "lhcbpsip":["(Bmass < 5.23) && (Mumumass < 3.75) && (3.686 - Mumumass) <", 0.16, 0.01, 0.30,0.01]  }  # LHCb psi



#allcut = {"Bvtxcl":["Bvtxcl >",0.16, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9995, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 7.0, 6.0, 12.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.2,0.8,2.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.8, 0.5, 3.5, 0.1]}
#allcut = {"Bvtxcl":["Bvtxcl >",0.11, 0.02, 0.22, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9996, 0.9990,0.9999, 0.0001], "Blxysig":["Blxysig >", 8.1, 6.0, 12.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.2,0.8,2.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",2.8, 0.5, 3.5, 0.1]}

#allcut = {"Bvtxcl":["Bvtxcl >",0.09, 0.06, 0.18, 0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.9992 ,1.0, 0.0001], "Blxysig":["Blxysig >", 8.1, 6.0, 11.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.6,1.2,2.7,0.1], "Trkdcasigbs":["Trkdcasigbs >",1.5, 0.3, 2.0, 0.1], "psimassh":["Mumumass <", 3.20, 0.040, 0.09,0.002],"psimassl":["3.097-Mumumass >", 0.18, 0.040, 0.30,0.02]}

#allcut = {"Bvtxcl":["Bvtxcl >",0.09, 0.06, 0.2,0.01], "Bcosalphabs2D":["Bcosalphabs2D >", 0.9997, 0.999 ,1.0, 0.0001], "Blxysig":["Blxysig >", 8.1, 6.0, 12.0, 0.1], "Tk_4vec.Pt()":["Tk_4vec.Pt() >", 1.8,1.0,3.6,0.1], "Trkdcasigbs":["Trkdcasigbs >",1.0, 0.0, 1.2, 0.1], "cdfjpsi":["abs(Bmass-Mumumass-2.182) >", 0.18, 0.06, 0.20,0.01],"cdfpsip":["abs(Bmass-Mumumass-1.593) >", 0.09, 0.06, 0.20,0.01]}





def cut1(thiskey):
    ss =")"
    for k, v in allcut.items():
        if k == thiskey:
            continue
        else:
            ss += "&&("
            thisv = v[0] + str(v[1])
            ss += thisv
            ss += ")"

    print ss 
    return ss


def cut2(thiskey, cc):
    ss = "("
    ss += allcut[thiskey][0] + str(cc)
    print ss
    return ss


def m_fom(t1, t2, thiskey, thiscut):
    N_sig = t1.GetEntries(cut2(thiskey, thiscut) + cut1(thiskey))
    N_inc = t2.GetEntries(cut2(thiskey, thiscut) + cut1(thiskey))
    print "This cut is " + thiskey + str(thiscut) + ", N_sig is " + str(N_sig) + " and " + "N_inc is " + str(N_inc)
    if N_inc == 0:
        N_inc += 1
    f = 0.006213575
    N1p = f*N_sig
    N2 = N_inc
    N01 = 1962416
    N02 = 897075
    fomerr2=N1p*(f*N01-N1p)/(N2*N01)+N1p**2*(N02-N2)/(4*N02*N2**2)
    fomerr=sqrt(fomerr2)
    fom = f * N_sig / sqrt(1.0 * N_inc);
    print "The FOM is "+ str(fom) + "+/-"+str(fomerr)
    return fom, fomerr


def good_fom(thiskey):
    global tg, ftitle

    ## OPT
    #dtfile = "../RootFiles/bfOPT/Data_2012_bfOPT_cut1.root"
    #mcfile = "../RootFiles/bfOPT/BToKMuMu_SignalMC_8TeV_v4_bfOPT_cut1.root"
    dtfile = "../RootFiles/bfOPT/Data_2012_8TeV_v4_OPTonlyVTXLoose.root"
    mcfile = "../RootFiles/bfOPT/MC_Signal_8TeV_v4_OPTonlyVTXLoose_NoGEN.root"
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+resonance-1.root"
    #mcfile = "../RootFiles/Files/MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1_NoGEN.root"
    ## veto
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT.root"
    #mcfile = "../RootFiles/Files/MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT_NoGEN.root"
    ## LHCb jpsi
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root"
    #mcfile = "../RootFiles/Files/MC_JpsiK_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root"
    ## LHCb psi
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root"
    #mcfile = "../RootFiles/Files/MC_Psi2SK_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root"

    cutop = allcut[thiskey][0]
    start = allcut[thiskey][2]
    lim = allcut[thiskey][3]
    step = allcut[thiskey][4]    
    
    f1 =  TFile(mcfile)
    t1 =  f1.Get("tree")

    f2 = TFile(dtfile)
    t2 = f2.Get("tree")

    cc = start;
    x, y, st, yer=array("d"), array("d"), array("d"), array("d")

    while (cc < lim):
        x.append(cc);
        print  "cut is " + str(cc)
        dd, derr = m_fom(t1, t2, thiskey, cc)
        y.append(dd)
        st.append(0.5*step)
        yer.append(derr)
        cc += step

    tg = TGraphErrors(len(x), x, y, st, yer)
    tg.SetMarkerStyle(8)
    tg.SetMarkerColor(2)
    #tg.SetMarkerSize(2)    
    tg.SetFillColor(4)
    tg.SetFillStyle(3006)
    # tg.SetMaximum(34.0)
    # tg.SetMinimum(30.0)
    ftitle = thiskey + ".pdf"
    tg.GetXaxis().SetTitle(cutop)
    tg.GetYaxis().SetTitle("FOM")
    tg.SetTitle("")
    

def plot_same(thisrange):

    ## OPT
    #dtfile = "../RootFiles/bfOPT/Data_2012_bfOPT_cut1.root"
    #mcfile = "../RootFiles/bfOPT/BToKMuMu_SignalMC_8TeV_v4_bfOPT_cut1.root"
    dtfile = "../RootFiles/bfOPT/Data_2012_8TeV_v4_OPTonlyVTXLoose.root"
    mcfile = "../RootFiles/bfOPT/MC_Signal_8TeV_v4_OPTonlyVTXLoose_NoGEN.root"
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+resonance-1.root"
    #mcfile = "../RootFiles/Files/MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1_NoGEN.root"
    ## veto
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT.root"
    #mcfile = "../RootFiles/Files/MC_Signal_8TeV_v4_cut0+HLT+Q+B2+resonance-1+OPT_NoGEN.root"
    ## LHCb jpsi
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root"
    #mcfile = "../RootFiles/Files/MC_JpsiK_8TeV_v4_cut0+HLT+Q+B2+JpsiK+OPT.root"
    ## LHCb psi
    #dtfile = "../RootFiles/Files/Data_2012_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root"
    #mcfile = "../RootFiles/Files/MC_Psi2SK_8TeV_v4_cut0+HLT+Q+B2+Psi2SK+OPT.root"

    f1 =  TFile(mcfile)
    t1 =  f1.Get("tree")

    f2 = TFile(dtfile)
    t2 = f2.Get("tree")



    thecut = "("+thisrange + cut1("null")
    t1.Draw("(Bmass-Mumumass)"+">>h1", thecut)
    t2.Draw("(Bmass-Mumumass)"+">>h2", thecut)
    hsig = gDirectory.Get('h1');
    hinc = gDirectory.Get('h2');
    hsig.SetLineColor(2)
    hinc.SetLineColor(4)
    hsig.Scale(0.006213575) 
    hinc.Draw()
    hsig.Draw("SAME")
    ftitle = "r1" + ".pdf"
    c1.Print(ftitle)
    
def plot_same2(thisrange):

    dtfile = "../RootFiles/Input.root"

    f1 = TFile(dtfile)
    t1 = f1.Get("tree")

    thecut = "("+thisrange + cut1("null")
    t1.Draw("(Bmass-Mumumass):Mumumass**2"+">>h2", thecut)
    hsig = gDirectory.Get('h2');
    #hsig.SetLineColor(2)
    #hsig.SetTitle(";q^2(GeV^2/c^2);DeltaM(GeV)")
    #hsig.Scale(0.006213575) 
    hsig.Draw("BOX")
    ftitle = "r1" + ".pdf"
    c1.Print(ftitle)



if __name__ == "__main__":
    c1 = TCanvas()
    #c1.SetLogy(1)
    for k in allcut.keys():
    #for k in ["antirad_lhcbjpsi"]:
    #for k in ["antirad_lhcbpsip"]:
    #for k in ["antirad_jpsi", "antirad_psip"]:
    #for k in ["jpsimass_cdf", "psimass_cdf"]:
        good_fom(k)
        tg.Draw("AP3")
        c1.Print(ftitle)
    #plot_same2("Mumumass<3.097") #R1
    #plot_same("Mumumass>3.686+3*Mumumasserr") #R3
    #plot_same("(Mumumass>3.097+3*Mumumasserr)&&(Mumumass<3.686-3*Mumumasserr)") #R2
    #plot_same("(Mumumass>3.097-5*Mumumasserr)&&(Mumumass<3.097+3*Mumumasserr)")
    #plot_same("(Mumumass>3.686-3*Mumumasserr)&&(Mumumass<3.686+3*Mumumasserr)")
