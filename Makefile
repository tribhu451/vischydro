
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

SRC        = main.cpp master.cpp grid.cpp cell.cpp cnvrt.cpp opt_glau.cpp fluid_info.cpp read_ic.cpp hydro.cpp \
             trancoeff.cpp hrr.cpp eos1.cpp eos0.cpp eos2.cpp mc_glau.cpp surface.cpp cornelius.cpp freeze.cpp Vector3D.cpp \
              sVector3D.cpp 
        #eos0.cpp  eos1.cpp cnvrt.cpp cell.cpp fluid.cpp hydro.cpp opt_glau.cpp mc_glau.cpp trancoeff.cpp \
	#fluid_info.cpp read_ic_from_file.cpp freeze.cpp Vector3D.cpp read_gubser_from_file.cpp\
        #sVector3D.cpp cornelius.cpp hrr.cpp
             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = vischydro
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
