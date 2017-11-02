.EXPORT_ALL_VARIABLES:

.PHONY: clean all

BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/lib
TARTSYS = /home/daq/tmp/anaroot/sources/Core/
COMMON_DIR = $(HOME)/TINAanalysis/common
ELOSS_DIR = $(HOME)/TINAanalysis/ELoss
KIN_DIR = $(HOME)/TINAanalysis/Reaction

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

CPP             = g++
CFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC

INCLUDES        = -I./inc -I$(COMMON_DIR) -I$(ELOSS_DIR) -I$(KIN_DIR) -I$(TARTSYS)/include
BASELIBS 	= -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR) -L$(TARTSYS)/.libs
ALLIBS  	=  $(BASELIBS) -lCommandLineInterface -lanacore -lEnergyLoss -lKinematics
LIBS 		= $(ALLIBS)
LFLAGS		= -g -fPIC -shared
CFLAGS += -Wl,--no-as-needed
LFLAGS += -Wl,--no-as-needed
CFLAGS += -Wno-unused-variable -Wno-write-strings

CFLAGS += -UKYUSHU

all: Turner CsIke

Turner:	Turner.cc
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) -o $(BIN_DIR)/$@ 

CsIke:	CsIke.cc
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) -o $(BIN_DIR)/$@ 

clean:
	@echo "Cleaning up"
	@rm -rf build doc
	@rm -f inc/*~ src/*~ scripts/*~ *~
