import ROOT
import numpy as np
from array import array
from time import time
NPOINTS = 10000
x = ROOT.RooRealVar('x', 'x', -5, 5)
mean = ROOT.RooRealVar('mean', 'mean', 0, -1, 1)
sigma = ROOT.RooRealVar('sigma', 'sigma', 1, 0.1, 2)
gaus = ROOT.RooGaussian('gaus', 'gaus', x, mean, sigma)
data = gaus.generate(ROOT.RooArgSet(x), 100)
nll = gaus.createNLL(data)
start = time()
nll_profile = nll.createProfile(ROOT.RooArgSet(mean))
frame_mean = mean.frame(ROOT.RooFit.Bins(NPOINTS))
nll_profile.plotOn(frame_mean, ROOT.RooFit.Precision(1))
time1 = time() - start
start = time()
m = ROOT.RooMinuit(nll)
m.setVerbose(False)
m.setPrintLevel(-1)
m.migrad()
nll_minimum = nll.getVal()
xs, ys = np.linspace(-1, 1, NPOINTS + 1), []
for x in xs:
	mean.setVal(x)
	mean.setConstant(True)
	m.migrad()
	ys.append(nll.getVal() - nll_minimum)
time2 = time() - start
gr = ROOT.TGraph(len(xs), array('f', xs), array('f', ys))
gr.SetMarkerStyle(20)
canvas = ROOT.TCanvas("canvas","trst",600,800)
frame_mean.addObject(gr, 'P')
frame_mean.Draw()
print "time automatic: ", time1
print "time manual:    ", time2




