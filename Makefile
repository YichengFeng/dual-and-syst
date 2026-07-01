###########################################################################
# Makefile — Build TDualGraph_cxx.so and TSystGraph_cxx.so with
#            hand-written Streamers
#
# Usage:
#   make            # builds both shared libraries
#   make clean
#
# Author: Yicheng Feng
# Email: fengyich@outlook.com
###########################################################################

CXX       := g++

ROOT_BINDIR := /home/fengyich/Programs/root_v6.34.10/root/bin
ROOT_INCDIR := $(shell $(ROOT_BINDIR)/root-config --incdir)
ROOT_LIBDIR := $(shell $(ROOT_BINDIR)/root-config --libdir)

CXXFLAGS  := -std=c++17 -O2 -fPIC -m64 -pthread -I. -I$(ROOT_INCDIR)
LDFLAGS   := -shared -L$(ROOT_LIBDIR) -lCore -lRIO -lHist -lGraf \
             -Wl,-rpath,$(ROOT_LIBDIR) -pthread -lm -ldl

ROOTCLING := $(ROOT_BINDIR)/rootcling

# rootcling's -c flag is deprecated/ignored.  We must NOT pass -std=...
# because rootcling interprets -s as its own option.
ROOTCLING_FLAGS := -I. -I$(ROOT_INCDIR)

# =========================================================================
#  Targets
# =========================================================================

all: TDualGraph_cxx.so TSystGraph_cxx.so TGraphErrorsOperator_cxx.so

clean:
	@rm -f TDualGraph_cxx.so TSystGraph_cxx.so TGraphErrorsOperator_cxx.so
	@rm -f TDualGraph.o TSystGraph.o TGraphErrorsOperator.o
	@rm -f TDualGraph_dict.o TDualGraph_dict.cxx TDualGraph_dict_rdict.pcm
	@rm -f TSystGraph_dict.o   TSystGraph_dict.cxx   TSystGraph_dict_rdict.pcm

# =========================================================================
#  Library 1 — TDualGraph (self-contained)
# =========================================================================

TDualGraph_cxx.so: TDualGraph.o TDualGraph_dict.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# =========================================================================
#  Library 2 — TSystGraph (self-contained; includes TDualGraph's
#  dictionary so Cling sees it without loading TDualGraph_cxx.so first)
# =========================================================================

TSystGraph_cxx.so: TSystGraph.o TDualGraph.o TDualGraph_dict.o TSystGraph_dict.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# =========================================================================
#  Library 3 — TGraphErrorsOperator (free functions, no dictionary)
# =========================================================================

TGraphErrorsOperator_cxx.so: TGraphErrorsOperator.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# =========================================================================
#  Compilation units
# =========================================================================

TDualGraph.o: TDualGraph.cxx TDualGraph.h DualGraph.h DualNumber.h DualMultiv.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

TSystGraph.o: TSystGraph.cxx TSystGraph.h SystGraph.h TDualGraph.h DualGraph.h \
              DualNumber.h DualMultiv.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

TGraphErrorsOperator.o: TGraphErrorsOperator.cxx TGraphErrorsOperator.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

TDualGraph_dict.o: TDualGraph_dict.cxx
	$(CXX) -c $(CXXFLAGS) -o $@ $<

TSystGraph_dict.o: TSystGraph_dict.cxx
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# =========================================================================
#  ROOT dictionaries
# =========================================================================

TDualGraph_dict.cxx: TDualGraph.h LinkDef_TDualGraph.h
	$(ROOTCLING) -f $@ -c $(ROOTCLING_FLAGS) TDualGraph.h LinkDef_TDualGraph.h

TSystGraph_dict.cxx: TSystGraph.h LinkDef_TSystGraph.h
	$(ROOTCLING) -f $@ -c $(ROOTCLING_FLAGS) TSystGraph.h LinkDef_TSystGraph.h

.PHONY: all clean
