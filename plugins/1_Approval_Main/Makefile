GCCFLAGS=-Wall -g
ROOTFLAGS=`root-config --cflags --libs`

TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c
SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof
FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

#all: sel fig fit
all: fit sel sel_noSB sel_OtherB sel_OBbfvtx
	@echo "Full build successful."

sel_OBbfvtx: SingleBuToKMuMuSelector_OBbfvtx.cc
	g++ $(SELECTORFLAGS) $< -o $@

sel_OtherB: SingleBuToKMuMuSelector_OtherB.cc
	g++ $(SELECTORFLAGS) $< -o $@

sel_noSB: SingleBuToKMuMuSelector_noSB.cc
	g++ $(SELECTORFLAGS) $< -o $@


sel: SingleBuToKMuMuSelector.cc
	g++ $(SELECTORFLAGS) $< -o $@

#fig: SingleBuToKMuMuFigures.cc
#	g++ $(FIGURESFLAGS) $< -o $@

fit: SingleBuToKMuMuFitters.cc
	g++ $(FITFLAGS) $< -o $@	

clean:
	rm -rf fit sel sel_noSBi sel_OtherB sel_OBbfvtx
#	rm -rf sel fig fit 
