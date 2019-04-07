GCCFLAGS=-Wall -g -lm
#ROOTFLAGS=`root-config --cflags --libs`

#TOOLSFLAGS=$(GCCFLAGS) -I$(ROOTSYS)/include -c
#SELECTORFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lProof 
#FIGURESFLAGS=$(GCCFLAGS) $(ROOTFLAGS)
#FITFLAGS=$(GCCFLAGS) $(ROOTFLAGS) -lRooFit -lRooFitCore -lMinuit

#all: optcut
#	@echo "Full build successful."

#optcut: cut_opt_I.cxx
#	g++ $(SELECTORFLAGS) $< -o $@



Target  = optcut
Objects = cut_optimize.o

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

#all:$(Target)
all: optcut anti

optcut: cut_opt_I.cxx cut_optimize.h
	g++   $(ROOTFLAGS) $(ROOTLIBS) $(GCCFLAGS) $< -o $@

anti: cut_anti.cxx cut_optimize.h
	g++   $(ROOTFLAGS) $(ROOTLIBS) $(GCCFLAGS) $< -o $@

clean:
	rm -rf  optcut anti
