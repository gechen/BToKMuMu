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

binfolder="./NLLMC/bin"
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
subfolderPattern = re.compile("^set[0-9]...$")
#subfolderPattern = re.compile("^set0001$")
#fitResultPattern = "wspace_FC_bin"+str(theBin)+".root"
#fitResultPattern = "setSummary.root"
#fitResultPattern1 = "h_setSummaryAfb.cc"
#fitResultPattern2 = "h_setSummaryFh.cc"
#fitResultPattern3 = "fitParameters"+str(theBin)+".txt"
#headfolder = "afb-2950_fh+8500"
for headfolder in headfolders:
    print "\nProcessing headfolder {0}".format(headfolder)
    #print "No fitting result is found. Skip!"
    for idx in os.listdir(binfolder+'/'+headfolder):
        #print idx
        fullpath = binfolder+'/'+headfolder+'/'+idx
        if not subfolderPattern.match(idx):
            continue
        for ii in os.listdir(fullpath):
            print ii
