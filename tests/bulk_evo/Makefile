
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

GSLLIBS      := $(shell gsl-config --libs)

#EXTRA_FLAGS = -D SIMPLE # EoS p=e/3
EXTRA_FLAGS  = -D TABLE # Laine EoS, tabulated

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) $(GSLLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp


OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
TARGET	   = 3Dhydro
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

cfiles : 
	rm -rf hydro_output/*
	rm -rf output_after/* 

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
