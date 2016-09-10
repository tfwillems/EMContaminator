##
## Makefile for all executables
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX -pthread #-pedantic -Wunreachable-code -Weverything

## To build static executables, run:
##   rm -f HipSTR BamSieve
##   make STATIC=1
## verify:
##   ldd HipSTR BamSieve
##
## To Create a static distribution file, run:
##   make static-dist
ifeq ($(STATIC),1)
LDFLAGS=-static
else
LDFLAGS=
endif

## Source code files, add new files to this list
SRC_COMMON = main.cpp em_analyzer.cpp error.cpp mathops.cpp snp.cpp stringops.cpp vcf_reader.cpp version.cpp

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)

BAMTOOLS_ROOT=bamtools
HTSLIB_ROOT=htslib

LIBS              = -L./ -lm -L$(HTSLIB_ROOT)/ -L$(BAMTOOLS_ROOT)/lib -lz
INCLUDE           = -I$(BAMTOOLS_ROOT)/src
BAMTOOLS_LIB      = $(BAMTOOLS_ROOT)/lib/libbamtools.a
HTSLIB_LIB        = $(HTSLIB_ROOT)/libhts.a

.PHONY: all
all: version EMContaminator
	rm version.cpp
	touch version.cpp

# Create a tarball with static binaries
.PHONY: static-dist
static-dist:
	rm -f HipSTR BamSieve
	$(MAKE) STATIC=1
	( VER="$$(git describe --abbrev=7 --dirty --always --tags)" ;\
	  DST="EMContaminator-$${VER}-static-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
            cp EMContaminator && \
            tar -czvf "$${DST}.tar.gz" "$${DST}" && \
            rm -r "$${DST}/" \
        )

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > version.cpp

# Clean the generated files of the main project only (leave Bamtools/htslib alone)
.PHONY: clean
clean:
	rm -f *.o *.d EMContaminator

# Clean all compiled files, including bamtools/hstlib
.PHONY: clean-all
clean-all: clean
	if test -d bamtools/build ; then \
		$(MAKE) -C bamtools/build clean ; \
		rm -rf bamtools/build ; \
	fi

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable
EMContaminator: $(OBJ_COMMON) $(BAMTOOLS_LIB) $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp $(BAMTOOLS_LIB)
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

# Rebuild BAMTools if needed
$(BAMTOOLS_LIB):
	git submodule update --init --recursive bamtools
	git submodule update --recursive bamtools
	( cd bamtools && mkdir build && cd build && cmake .. && $(MAKE) )

# Rebuild htslib library if needed
$(HTSLIB_LIB):
	git submodule update --init --recursive htslib
	git submodule update --recursive htslib
	cd htslib && $(MAKE)
