#!/usr/bin/env python
"""
Main program for fitting 

"""

import sys
import linecache
from tls import *

## total number of scan seeds
N = 25
def main():
    args = sys.argv[1:]  ### "1:" after ./scan.py 
    if args[0] == 'ToyVali':
        for ibin in range(0,11):
            if ibin is (3 or 5): continue
            #print ibin
            file_in = ("./Toy/bin%s/Summary.root"%ibin)
            #print ("./Toy/bin%s/Summary.root"%ibin)
            return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_0':
        ibin = 0
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_1':
        ibin = 1
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_2':
        ibin = 2
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_4':
        ibin = 4
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_6':
        ibin = 6
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_7':
        ibin = 7
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_8':
        ibin = 8
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_9':
        ibin = 9
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    elif args[0] == 'ToyVali_10':
        ibin = 10
        file_in = ("./Toy/bin%s/Summary.root"%ibin)
        return call_fit_2D(file_in, ibin)
    else:
        raise NameError(args)

def call_fit_2D(infile, ibin):
    vini = genseed(N)
    for idx in xrange(N):
        afb_i, fh_i = vini[idx]
        fh_hat_i = fh_to_hat(fh_i)
        afb_hat_i = afb_to_hat(afb_i, fh_i)
        os.system("./fit angular2D %s %s %s %s %s >& AllCuts_v3/2D_%s_%s.log" %(infile,ibin, str(idx), str(afb_hat_i), str(fh_hat_i), ibin, str(idx) ))
        #os.system("./fit angular2D %s %s %s %s %s >& AllCuts_v3/2D_%s_%s.log" %(infile,ibin, str(idx), str(afb_i), str(fh_i), ibin, str(idx) ))
if __name__ == '__main__':    ### main function
    main()
















