#!/usr/bin/env python
"""
Main program for fitting 

"""

__author__ = "Geng CHEN <geng.chen@cern.ch>"
__copyright__ = "Copyright (c) Geng CHEN"


import sys 
from tls import * 

def main():
    #dirname = 'AllCuts_v3'
    args = sys.argv[1:]  ### "1:" after ./fit.py 
    if args[0] == 'angular_reco':
        return angular_reco(args[1:])
    elif args[0] == 'angular2D_1a_Sm':
        return angular2D_1a_Sm(args[1:])
    elif args[0] == 'angular2D':
        return angular2D(args[1:])   ### after function
    elif args[0] == 'angular2D_Toy':
        return angular2D(args[1:])   ### after function
    else:
        raise NameError(args)

def angular_reco(args):
    dirname = 'AllCuts_v3'
    datatype = args[0] ### data type
    label = args[1]    ### label
    ibin = args[2]    ### iBin
    index = args[3]    ### index 2015-05-04
    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')   ## b sub
    figname = 'angular_reco' + '_' + ibin + '_' + index

    if batch:
        jobname = 'fit'
        queue = '1nd'   ## 1 day
        cmd = create_batch_cmd(main='./fit')   ## 
        bashfile = set_file('/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/FCToyFitJobs/', dirname, figname, '.sh', test=test)
        pwd = os.getcwd()  ## current path
        pre = 'cd %s' % pwd 
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file('./', dirname, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:  ### 10 jobs
            logfile += '.%I'

        #test = True 
        test = False 
        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
        
def angular2D_1a_Sm(args):
    dirname = 'AllCuts_v3'
    datatype = args[0] ### data type
    label = args[1]    ### label
    ibin = args[2]    ### iBin
    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')   ## b sub
    figname = 'angular2D_1a_Sm' + '_' + ibin

    if batch:
        jobname = 'fit'
        queue = '1nd'   ## 1 day
        cmd = create_batch_cmd(main='./fit')   ## 
        bashfile = set_file('/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/FCToyFitJobs/', dirname, figname, '.sh', test=test)
        pwd = os.getcwd()  ## current path
        pre = 'cd %s' % pwd 
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file('./', dirname, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:  ### 10 jobs
            logfile += '.%I'

        #test = True 
        test = False 
        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
        
def angular2D(args):
    dirname = 'AllCuts_v3'
    datatype = args[0] ### data type
    label = args[1]    ### label
    ibin = args[2]    ### iBin
    index = args[3]    ### index 2015-05-04
    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')   ## b sub
    figname = 'angular2D' + '_' + ibin + '_' + index

    if batch:
        jobname = 'fit'
        queue = '1nd'   ## 1 day
        cmd = create_batch_cmd(main='./fit')   ## 
        bashfile = set_file('/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/FCToyFitJobs/', dirname, figname, '.sh', test=test)
        pwd = os.getcwd()  ## current path
        pre = 'cd %s' % pwd 
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file('./', dirname, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:  ### 10 jobs
            logfile += '.%I'

        #test = True 
        test = False 
        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
        
def angular2D_Toy(args):
    dirname = 'AllCuts_v3'
    datatype = args[0] ### data type
    label = args[1]    ### label
    ibin = args[2]    ### iBin
    index = args[3]    ### index 2015-05-04
    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')   ## b sub
    figname = 'angular2D' + '_' + ibin + '_' + index

    if batch:
        jobname = 'fit'
        queue = '2nm'   ## 1 day
        cmd = create_batch_cmd(main='./fit')   ## 
        bashfile = set_file('/afs/cern.ch/user/g/gechen/gechen/work/BToKMuMu/Ana/FCToyFitJobs/', dirname, figname, '.sh', test=test)
        pwd = os.getcwd()  ## current path
        pre = 'cd %s' % pwd 
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file('./', dirname, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:  ### 10 jobs
            logfile += '.%I'

        #test = True 
        test = False 
        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
        
if __name__ == '__main__':    ### main function
    main()

    
