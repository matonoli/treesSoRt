#
# Modified PYTHIA Makefile
#
# 1) how to compile PYTHIA
# 
# Download and unpack the latest PYTHIA from:
# http://home.thep.lu.se/Pythia/
# 
# In the directory where you have unpacked it:
# ./configure --enable-shared
# make 
#
# Currently you have by hand to remove the pedantic options in:
# $(MYPYTHIA)/Makefile.inc
# to make charged_particle_tree compile (step 3)
#
# 2) Variables to define in .bashrc before you can compile charged_particle_tree
#
# export MYPYTHIA=/home/pchristi/work/pythia/pythia8210
# (This of course has to point to the PYTHIA version you download)
# export PYTHIA8DATA=${MYPYTHIA}/share/Pythia8/xmldoc
# if [ -z "$LD_LIBRARY_PATH" ]
# then                        
#     export LD_LIBRARY_PATH=${MYPYTHIA}/lib
# else                                     
#     export LD_LIBRARY_PATH=${MYPYTHIA}/lib:${LD_LIBRARY_PATH}
# fi
#
# 3) how to compile the program charged_particle_tree
#
# Note that to do this, ROOT has to be loaded
#
# make charged_particle_tree

SHELL = /bin/sh

-include $(MYPYTHIA)/Makefile.inc

# PYTHIA variables
PYTHIA_INCDIR=$(MYPYTHIA)/include
PYTHIA_LIBDIR=$(MYPYTHIA)/lib

# ROOT variables (ROOTCFLAGS also includes include path)
ROOTCFLAGS	= $(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --glibs)

# There is no default behaviour, so remind user.
all:
	@echo "Usage: make XXX, where XXX.cc is your program"
	@echo $(PYTHIA_LIBDIR)
	@echo $(CXX)
	@echo $(CXXFLAGS)
	@echo $(ROOTCFLAGS)
	@echo $(PYTHIA_INCDIR)
	@echo $(ROOTLIBS)

# Create an executable for one of the normal test programs
%:	%.cc $(PYTHIA_LIBDIR)/libpythia8.so #dependencies
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -I$(PYTHIA_INCDIR) \
	$@.cc -o $@.exe \
	-L$(PYTHIA_LIBDIR) -lpythia8 \
	$(ROOTLIBS) -lEG \
	TransverseSpherocity/TransverseSpherocity_cxx.so #-g


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -f *.exe
	rm -f *~; rm -f \#*; rm -f core*
