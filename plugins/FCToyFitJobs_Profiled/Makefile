GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

all: fit nll the1 test1 nll2 nll3 cov
	@echo "Full build successful."


#sel: SingleBuToKMuMuSelector.cc
#	g++ $(SELECTORFLAGS) $< -o $@

#fig: SingleBuToKMuMuFigures.cc
#	g++ $(FIGURESFLAGS) $< -o $@

## Toys
fit: SingleBuToKMuMuFitters.cc
	g++ $(FITFLAGS) $< -o $@	

## bin7_new
nll: SingleBuToKMuMuFitters_NLL.cc
	g++ $(FITFLAGS) $< -o $@	
#### bin2,4,10
nll2: SingleBuToKMuMuFitters_NLL2.cc
	g++ $(FITFLAGS) $< -o $@	
nll3: SingleBuToKMuMuFitters_NLL3.cc
	g++ $(FITFLAGS) $< -o $@	
#### bin2,4,10
## Coverage
cov: SingleBuToKMuMuFitters_Cov.cc
	g++ $(FITFLAGS) $< -o $@	


## bin7_old, bin7_points
the1: SingleBuToKMuMuFitters_1.cc
	g++ $(FITFLAGS) $< -o $@	

## bin7_test
test1: SingleBuToKMuMuFitters_backFor1.cc
	g++ $(FITFLAGS) $< -o $@	

## 2D_NLL
#the2: SingleBuToKMuMuFitters_2.cc
#	g++ $(FITFLAGS) $< -o $@	


clean:
	rm -rf fit nll the1 test1 nll2 nll3 cov
