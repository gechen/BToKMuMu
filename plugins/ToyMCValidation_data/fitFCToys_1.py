#!/usr/bin/env python

import os
import subprocess
import time
import sys
import re
import shutil
from threading import Thread

number_of_threads = 11

class transfer_thread(Thread):
    def __init__ (self,headfolder,subfolder):
        Thread.__init__(self)
        self.headfolder = headfolder
        self.subfolder  = subfolder
        self.fullpath   = headfolder+"/"+subfolder
        self.status=-1

    def run(self):
        print time.ctime()
        print "Processing folder",self.headfolder,self.subfolder
        sys.stdout.flush()
        cmd=["./fit", "angular2D_Toy"  ,self.fullpath+"/*Toy*.root", theBin]
        #cmd=["./fit", "angular2D_Toy"  ,self.fullpath+"/*Toy*.root", "--iallpath="+self.headfolder, "--oallpath="+self.fullpath, "--iCombBkgWspacepath="+self.fullpath, "--keeplog", "--keepparam"]
        # cmd=["echo", "Dry run in "+self.fullpath]
        # if "UnderProcessing" in os.listdir(self.fullpath):
            # cmd=["echo", "Skip UnderProcessing tag in "+self.fullpath]
        # print cmd
        time.sleep(3)
        operation = subprocess.Popen(cmd)
        operation.wait()
        self.status=0

proclist = []

binfolder="./Toy/bin"
if len(sys.argv) == 2 :
    theBin = sys.argv[1]
    binfolder = binfolder+sys.argv[1]
    print "Running with "+binfolder
else:
    exit(1)

headfolders = []
for idx in os.listdir(binfolder):
    if 'afb' not in idx: continue
    headfolders.append(idx)


# Start loop
#subfolderPattern = re.compile("^set[0-9]...$")
subfolderPattern = re.compile("^set[0-9]...$")
#fitResultPattern = "wspace_FC_bin"+str(theBin)+".root"
#fitResultPattern = "setSummary.root"
fitResultPattern1 = "h_setSummaryAfb.cc"
fitResultPattern2 = "h_setSummaryFh.cc"
#fitResultPattern3 = "wspace_FC_bin8.root"
fitResultPattern4 = "signalToy_Bin"+str(theBin)+".root"
fitResultPattern5 = "combBkgToy_Bin"+str(theBin)+".root"
for headfolder in headfolders:
    print "Processing headfolder {0}".format(headfolder)
    for idx in os.listdir(binfolder+'/'+headfolder):
        #print idx
        fullpath = binfolder+'/'+headfolder+'/'+idx
        if not subfolderPattern.match(idx):
            continue

        fullpath1 = binfolder+'/'+headfolder
#        if fitResultPattern in os.listdir(fullpath1):
#            os.remove(fullpath1+"/setSummary.root")
        if fitResultPattern1 in os.listdir(fullpath1):
            os.remove(fullpath1+"/h_setSummaryAfb.cc")
        if fitResultPattern2 in os.listdir(fullpath1):
            os.remove(fullpath1+"/h_setSummaryFh.cc")
        #if fitResultPattern3 in os.listdir(fullpath):
        #    os.remove(fullpath+"/wspace_FC_bin8.root")
        if fitResultPattern4 in os.listdir(fullpath1):
            os.remove(fullpath1+"/"+fitResultPattern4)
        if fitResultPattern5 in os.listdir(fullpath1):
            os.remove(fullpath1+"/"+fitResultPattern5)
        fitResultPattern6 = "signalToy_Bin"+str(theBin)+"_"+str(idx)+".root"
        fitResultPattern7 = "combBkgToy_Bin"+str(theBin)+"_"+str(idx)+".root"
        fitResultPattern8 = "fitParameters"+str(theBin)+".txt"
        if fitResultPattern8 in os.listdir(fullpath):
            os.remove(fullpath+"/"+"signalToy_Bin"+str(theBin)+"_"+str(idx)+".root")
            os.remove(fullpath+"/"+"combBkgToy_Bin"+str(theBin)+"_"+str(idx)+".root")

#        if fitResultPattern in os.listdir(fullpath):
#            print "Old fitting result "+fullpath+'/'+fitResultPattern+" is found. Skip!"
#            #os.remove(fullpath+"/setSummary.root")
#            #os.remove(fullpath+"/"+fitResultPattern)
#            if "UnderProcessing" in os.listdir(fullpath):
#                os.remove(fullpath+"/UnderProcessing")
#
#        elif "UnderProcessing" in os.listdir(fullpath):
#            # Clean the tags during development
#            #if fitResultPattern not in os.listdir(fullpath):
#                # os.remove(fullpath+"/UnderProcessing")
#            #if "UnderProcessing" in os.listdir(fullpath):
#            #    os.remove(fullpath+"/UnderProcessing")
#
#            print "Tag for running job is found in "+fullpath+". Skip!"
        else:
            print "Processing subfolder {0}".format(idx)
            #file=open(fullpath+"/UnderProcessing",'w')
            #file.close()

            #n_active_proc = number_of_threads
            #while n_active_proc>=number_of_threads:
            #    n_active_proc = 0
            #    for proc in proclist:
            #        if proc.status<0: n_active_proc = n_active_proc + 1
            #    time.sleep(1)

#            current = transfer_thread(binfolder+'/'+headfolder,idx)
#            proclist.append(current)
#            current.start()
#            time.sleep(1)

