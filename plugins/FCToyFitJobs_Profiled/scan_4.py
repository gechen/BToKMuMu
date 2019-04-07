#!/usr/bin/env python
"""
Main program for fitting 

"""

import sys
import linecache
from tls import *

## total number of scan seeds
N = 100
NN = 100
NNN = 1000
bandFh  = 0.10
bandAfb = 0.05
def main():
    args = sys.argv[1:]  ### "1:" after ./scan.py 
    if args[0] == 'angular_gen':
        ## this is input file
        file_in = "../RootFiles/MC_GENOnly/BToKMuMu_GENOnly_8TeV_genonly_v5+7.root"
        ibin = args[1]
        return call_fit_gen(file_in, ibin)
    elif args[0] == 'angular_gen_R':
        file_in = "../CutASelection/RootFiles/MC_Signal_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_gen_R(file_in, ibin)
    elif args[0] == 'angular_reco_D':
        file_in = "../CutASelection/RootFiles/MC_Signal_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_reco_D(file_in, ibin)
    elif args[0] == 'angular_reco_R':
        file_in = "../CutASelection/RootFiles/MC_Signal_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_reco_R(file_in, ibin)
    elif args[0] == 'angular_reco':
        file_in = "../CutASelection/RootFiles/MC_Signal_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_reco(file_in, ibin)
    elif args[0] == 'angular2D':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'angular2D_Limited':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D_L(file_in, ibin)
    elif args[0] == 'angular2D_NLL':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D_NLL(file_in, ibin)
    elif args[0] == 'angular2D_Profiled_Afb':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D_Afb(file_in, ibin)
    elif args[0] == 'angular2D_Profiled_Afb_refit':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D_Afb_refit(file_in, ibin)
    elif args[0] == 'angular2D_Profiled_Fh':
        file_in = "../CutASelection/RootFiles/Data_2012_8TeV_v4_AllCut.root"
        ibin = args[1]
        return call_fit_2D_Fh(file_in, ibin)
    elif args[0] == 'accXrecoEffLimited':
        file_in = "fitParameters/eff_par_7.txt"
        ibin = args[1]
        return call_fit_Eff_M(file_in, ibin)
    else:
        raise NameError(args)

def call_fit_gen(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit.py angular_gen mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))
        #os.system("./fit angular_gen %s %s %s %s %s >& AllCuts_v3/gen_%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), ibin, str(idx) ))

def call_fit_gen_R(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit.py angular_gen_R mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))

def call_fit_reco_D(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit.py angular_reco_D mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))

def call_fit_reco_R(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit.py angular_reco_R mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))

def call_fit_reco(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit.py angular_reco mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))

def call_fit_2D(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        #os.system("./fit.py angular2D data %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))
        os.system("./fit angular2D %s %s %s %s %s >& AllCuts_v3/2D_%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), ibin, str(idx) ))

def call_fit_2D_NLL(infile, ibin):
    vini = genNLL(NNN)
    print "xxxxxxxx"
    for idx in xrange(456,NNN):
        afb_i, fh_i = vini[idx]
        if (abs(afb_i)>0.5*fh_i): print afb_i, fh_i
        print "xxx > "+str(idx)+" < xxx"
#        if (abs(afb_i)<0.5*fh_i): print afb_i, fh_i
        os.system("./fit angular2D_NLL %s %s %s %s %s >& AllCuts_v3/2D_bin%s_%s.log" %(infile,ibin, str(idx), str(afb_i), str(fh_i), ibin, str(idx) ))

def call_fit_2D_Afb(infile, ibin):
    print str(bandFh),str(ibin)
    stepFh = bandFh/NN
    for ifh in xrange(NN):
        fh_ii = ifh * stepFh + stepFh/2.
        print str(ifh)+" >>> fix fh = "+str(fh_ii)
        vini = genafb(N, fh_ii)
        for idx in xrange(N):
            afb_i, fh_i = vini[idx]
            fh_hat_i = fh_to_hat(fh_i)
            afb_hat_i = afb_to_hat(afb_i, fh_i)
            #os.system("./fit angular2D_Profiled_Afb %s %s %s %s %s %s >& AllCuts_v3/2D_bin%s_Fh%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), str(ifh), ibin, str(ifh), str(idx) ))
            os.system("./fit angular2D_Profiled_Afb %s %s %s %s %s %s >& AllCuts_v3/2D_bin%s_Fh%s_%s.log" %(infile,ibin, str(idx), str(afb_i), str(fh_i), str(ifh), ibin, str(ifh), str(idx) ))

def call_fit_2D_Afb_refit(infile, ibin):
    NN = int(bandFh/stepFh)
    for ifh in xrange(NN):
        fh_ii = ifh * stepFh + stepFh/2.
        print str(ifh)+" >>> fix fh = "+str(fh_ii)
        os.system("./fit angular2D_Profiled_Afb %s %s refit %s  >& AllCuts_v3/2D_bin%s_Fh%s.log" %(infile,ibin, str(ifh), ibin, str(ifh) ))

def call_fit_2D_Fh(infile, ibin):
    stepAfb = 2.*bandAfb/NN
    for iafb in xrange(NN):
        afb_ii = iafb * stepAfb - bandAfb + stepAfb/2.
        print str(iafb)+" >>> fix afb = "+str(afb_ii)
        vini = genfh(N, afb_ii)
        for idx in xrange(N):
            afb_i, fh_i = vini[idx]
            fh_hat_i = fh_to_hat(fh_i)
            afb_hat_i = afb_to_hat(afb_i, fh_i)
            #if afb_i != hat_to_afb(afb_hat_i, fh_hat_i): 
            #    continue
            #print str(iafb)+">>> fix afb = "+str(afb_ii)+" afb_i "+str(afb_i)+"; fh_i "+str(fh_i)+"; afb_hat_i "+str(afb_hat_i)+"; fh_hat_i "+str(fh_hat_i) +" afb_i "+str(hat_to_afb(afb_hat_i, fh_hat_i))
            #os.system("./fit angular2D_Profiled_Fh %s %s %s %s %s %s >& AllCuts_v3/2D_bin%s_Afb%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), str(iafb), ibin, str(iafb), str(idx) ))
            #os.system("./fit angular2D_Profiled_Fh %s %s %s %s %s %s >& AllCuts_v3/2D_bin%s_Afb%s_%s.log" %(infile,ibin, str(idx), str(afb_i), str(fh_hat_i), str(iafb), ibin, str(iafb), str(idx) ))
            os.system("./fit angular2D_Profiled_Fh %s %s %s %s %s %s >& AllCuts_v3/2D_bin%s_Afb%s_%s.log" %(infile,ibin, str(idx), str(afb_i), str(fh_i), str(iafb), ibin, str(iafb), str(idx) ))







def call_fit_2D_L(infile, ibin):
    for idx in xrange(N):
        vini = genseed(N)
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit angular2D_Limited %s %s %s %s %s >& AllCuts_v3/2D_%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), ibin, str(idx) ))
        #os.system("./fit angular2D_Limited %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s >& AllCuts_v3/2D_%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), str(EffP1), str(EffP2), str(EffP3), str(EffP4), str(EffP5), str(EffP6), str(EffP7), str(EffP8), str(EffP9), str(EffP10), ibin, str(idx) ))

def call_fit_Eff_L(infile, ibin):
    for idx in xrange(1):
        effP = genseed_E8(1)
        EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10 = effP[idx]
        print EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10
        os.system("./fit accXrecoEffLimited %s %s %s %s %s %s %s %s %s %s %s %s >& AllCuts_v3/eff_%s_%s.log" %(infile,ibin, str(EffP1), str(EffP2), str(EffP3), str(EffP4), str(EffP5), str(EffP6), str(EffP7), str(EffP8), str(EffP9), str(EffP10), ibin, str(idx) ))


def call_fit_Eff_M(infile, ibin):
    #NLines = len(open('fitParameters/eff_par_2.txt','rU').readlines())
    NLines = len(open(infile,'rU').readlines())
    #print NLines
    for n in range(NLines):
        #count = linecache.getline('fitParameters/eff_par_2.txt',n+1)
        count = linecache.getline(infile,n+1)
        #print count
        with open('fitParameters/temp.txt', 'w') as f:
            f.write(count)
        file = open('fitParameters/temp.txt', 'r')
        lines = [[float(digit) for digit in line.split()] for line in file]
        #print lines[0][2]
        EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7 = lines[0]
        EffP8, EffP9, EffP10 = [0., 0., 0.]
        os.system("./fit accXrecoEffLimited %s %s %s %s %s %s %s %s %s %s %s %s %s >& AllCuts_v3/eff_%s_%s.log" %(infile,ibin, str(EffP1), str(EffP2), str(EffP3), str(EffP4), str(EffP5), str(EffP6), str(EffP7), str(EffP8), str(EffP9), str(EffP10), str(n+1), ibin, str(n+1) ))

if __name__ == '__main__':    ### main function
    main()
















