'''
#!/usr/bin/env python
'''
import random, os
from subprocess import call
from random import gauss, uniform
from math import *

'''
## total number of scan seeds
N = 1
'''

def genseed(N):
## generate the random seed for initial values
    sigma = 0.4

    ini_afb = []
    ##ini_fh = [random.random() for i in xrange(N)]
    ini_fh = [random.uniform(0,3) for i in xrange(N)]  ## 2015-09-15 ++
    
    ## reweighting to increase scan efficiency
    j = 0
    while len(ini_afb) < N:
        value = gauss(0.0, sigma)
        if -1.0 < value < 1.0 and abs(value) < 0.5*ini_fh[j]:
            ini_afb.append(value)
            j += 1
    ini_v = zip(ini_afb, ini_fh)
    return ini_v

def fh_to_hat(fh_in):
    vhat = tan ( (fh_in - 3./2.) * pi / 3. );
    return vhat

def afb_to_hat(afb_in, fh_in):
    vhat = tan (  afb_in * pi / fh_to_hat(fh_in) );
    return vhat

def hat_to_fh(fh_in):
    vhat = 3./2. + 3. * atan( fh_in  ) / pi;  ## hat to fh
    return vhat

def hat_to_afb(afb_in, fh_in):
    vhat = (1. * atan( afb_in ) / pi) * ( 3./2. + 3. * atan( fh_to_hat(fh_in)  ) / pi );  ## hat to afb
    return vhat



def genfh(N, afb_in):
## generate the random seed for initial values
    ini_afb = []
    ini_fh = [random.uniform(abs(2*afb_in),abs(5*afb_in)) for i in xrange(N)]  ## 2015-09-15 ++
    while len(ini_afb) < N:
        ini_afb.append(afb_in)
    ini_v = zip(ini_afb, ini_fh)
    return ini_v

def genafb(N, fh_in):
## generate the random seed for initial values
    sigma = 0.4
    ini_afb = []
    ini_fh  = []
    ## reweighting to increase scan efficiency
    while len(ini_afb) < N:
        value = gauss(0.0, sigma)
        if -1.0 < value < 1.0 and abs(value) < 0.5*fh_in:
            ini_afb.append(value)
            ini_fh.append(fh_in)
    ini_v = zip(ini_afb, ini_fh)
    return ini_v

def genNLL(N):
## generate the random seed for initial values
    sigma = 0.4

    ini_afb = []
    ##ini_fh = [random.random() for i in xrange(N)]
#    ini_fh = [random.uniform(0,0.4) for i in xrange(N)]  ## bin6, 7, 8, 
    ini_fh = [random.uniform(0,0.12) for i in xrange(N)]  ## bin2, 4, 10, 
#    ini_fh = [random.uniform(0,1.2) for i in xrange(N)]  ## bin0, 1,  
#    ini_fh = [random.uniform(0,0.8) for i in xrange(N)]  ## bin9,  
    
    ## reweighting to increase scan efficiency
    j = 0
    while len(ini_afb) < N:
        value = gauss(0.0, sigma)
        if abs(value) < 0.5*ini_fh[j]:
            ini_afb.append(value)
            j += 1
    ini_v = zip(ini_afb, ini_fh)
    for ii in xrange(N): 
        if (abs(ini_afb[ii])>0.5*ini_fh[ii]): print ini_afb[ii], ini_fh[ii]
    return ini_v

'''
def call_fit(infile, ibin):
    for idx in xrange(N):
        vini = genseed()
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)

#        os.system(./fit angular2D file_in  1  afb_i fh_i Index  test)
#        os.system("./fit angular2D   %s %s %s %s %s test >  ./test_%s.log" %(infile,ibin, str(afb_hat_i), str(fh_hat_i), str(idx),ibin))
        #os.system("./fit.py angular2D data %s %s %s %s test -b" %(infile,ibin, str(afb_hat_i), str(fh_hat_i)))
        #os.system("./fit.py angular_gen_R mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))
        os.system("./fit.py angular_reco mc %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))
        #os.system("./fit.py angular2D data %s %s %s %s %s -b" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i) ))


## the main program now
## this is input file
#file_in = "../RootFiles/AllCut/Data_KMuMu_AllCuts_v3.root"
file_in = "../RootFiles/AllCut/MC_Signal_AllCuts_v3.root"

for ibin in ["0"]:
#for ibin in ["0","1","2","4","6","7","8","9","10"]:
    call_fit(file_in, ibin)
'''


