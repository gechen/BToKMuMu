'''
#!/usr/bin/env python
'''
import random, os
from subprocess import call
from random import gauss, uniform
from math import *
#from diff import *
from Parameters4 import *

#file = open('fitParameters/fitParameters4.txt', 'r')
#print file.read()

def genseed_E4(N):
## generate the random seed for initial values
    sigma1 = 0.1*accXrecoEffErr[1]
    mean1  = accXrecoEff[1]
    sigma2 = 0.1*accXrecoEffErr[2]
    mean2  = accXrecoEff[2]
    sigma3 = 0.1*accXrecoEffErr[3]
    mean3  = accXrecoEff[3]
    sigma4 = 0.1*accXrecoEffErr[4]
    mean4  = accXrecoEff[4]
    sigma5 = 0.1*accXrecoEffErr[5]
    mean5  = accXrecoEff[5]
    sigma6 = 0.1*accXrecoEffErr[6]
    mean6  = accXrecoEff[6]
    sigma7 = 0.1*accXrecoEffErr[7]
    mean7  = accXrecoEff[7]
    sigma8 = 0.1*accXrecoEffErr[8]
    mean8  = accXrecoEff[8]
    sigma9 = 0.1*accXrecoEffErr[9]
    mean9  = accXrecoEff[9]
    sigma10 = 0.1*accXrecoEffErr[10]
    mean10  = accXrecoEff[10]
    EffP1 = []
    EffP2 = []
    EffP3 = []
    EffP4 = []
    EffP5 = []
    EffP6 = []
    EffP7 = []
    EffP8 = []
    EffP9 = []
    EffP10 = []
    ## reweighting to increase scan efficiency
    while len(EffP1) < 1: 
        j = 0
        #effp1 = [random.gauss(mean1,sigma1) for i in xrange(N)]
        effp1 = random.gauss(mean1,sigma1)
        effp2 = random.gauss(mean2,sigma2)
        effp3 = random.gauss(mean3,sigma3)
        effp4 = random.gauss(mean4,sigma4)
        effp5 = random.gauss(mean5,sigma5)
        effp6 = random.gauss(mean6,sigma6)
        effp7 = random.gauss(mean7,sigma7)
        effp8 = random.gauss(mean8,sigma8)
        effp9 = random.gauss(mean9,sigma9)
        effp10= random.gauss(mean10,sigma10)
        for y in xrange(-10, 11):
            x = 0.1 * y 
            #print x
            #value1 =  effp1 * exp(-0.5*(  pow (((x-effp2)/effp3),2) ))  
            #print value1
            #value2 =  effp4 + effp5 * x + effp6 * pow(x,2) + effp7 * pow(x,3) + effp8 * pow(x,4) + effp9 * pow(x,5) + effp10 * pow(x,6)  
            #print value2
            #value = value1 * value2
            value =  effp1 + effp2 * x + effp3 * pow(x,2) + effp4 * pow(x,3) + effp5 * pow(x,4) + effp6 * pow(x,5) + effp7 * pow(x,6)  
            #print value
            if value < 0.0:
                y += 1
            elif value < 0.0 and x == 1:
                return 0;
            elif value > 0.0 and x == 1:
            #elif value > 0.0 :
                EffP1.append(effp1)
                EffP2.append(effp2)
                EffP3.append(effp3)
                EffP4.append(effp4)
                EffP5.append(effp5)
                EffP6.append(effp6)
                EffP7.append(effp7)
                EffP8.append(effp8)
                EffP9.append(effp9)
                EffP10.append(effp10)
                j += 1
    #print EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10
    #ini_effP = [random.gauss(accXrecoEff[i],0.1*accXrecoEffErr[i]) for i in xrange(N)]
    ini_eff = zip(EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10 )
    return ini_eff

#for idx in xrange(3):
#    effP1 = genseed_E(3)
#    EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10 = effP1[idx]
#    print EffP1, EffP2, EffP3, EffP4, EffP5, EffP6, EffP7, EffP8, EffP9, EffP10



