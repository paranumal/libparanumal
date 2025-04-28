#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

define STRESS_HELP_MSG

Stress solver makefile targets:

   make stressMain (default)
   make lib
   make clean
   make clean-libs
   make clean-kernels
   make realclean
   make info
   make help
   make test

Usage:

make stressMain
   Build stressMain executable.
make lib
   Build libstress.a solver library.
make clean
   Clean the stressMain executable, library, and object files.
make clean-libs
   In addition to "make clean", also clean needed libraries.
make clean-kernels
   In addition to "make clean-libs", also cleans the cached OCCA kernels.
make realclean
   In addition to "make clean-kernels", also clean 3rd party libraries.
make info
   List directories and compiler flags in use.
make help
   Display this help message.
make test
   Run tests.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter stressMain lib clean clean-libs clean-kernels \
                realclean info help test, $(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${STRESS_HELP_MSG})
endif
endif

ifndef LIBP_MAKETOP_LOADED
ifeq (,$(wildcard ../../make.top))
$(error cannot locate ${PWD}/../../make.top)
else
include ../../make.top
endif
endif

#libraries
STRESS_LIBP_LIBS=parAlmond linearSolver mesh parAdogs ogs linAlg prim core

#includes
INCLUDES=${LIBP_INCLUDES} \
				 -I.
#defines
DEFINES =${LIBP_DEFINES} \
         -DLIBP_DIR='"${LIBP_DIR}"'

#.cpp compilation flags
STRESS_CXXFLAGS=${LIBP_CXXFLAGS} ${DEFINES} ${INCLUDES}
#.cu compilation flags
STRESS_NVCCFLAGS=${LIBP_NVCCFLAGS} ${DEFINES} ${INCLUDES}

#link libraries
LIBS=-L${LIBP_LIBS_DIR} $(addprefix -l,$(STRESS_LIBP_LIBS)) \
     ${LIBP_LIBS}

#link flags
LFLAGS=${STRESS_CXXFLAGS} ${LIBS}

#object dependancies
DEPS=$(wildcard *.hpp) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.hpp)

CXXSRC =$(wildcard src/*.cpp) 
NVCCSRC =$(wildcard src/*.cu) 

CXXOBJS=$(CXXSRC:.cpp=.o)
NVCCOBJS=$(NVCCSRC:.cu=.o)

.PHONY: all libp_libs clean clean-libs \
		clean-kernels realclean help info

all: stressMain

lib: libstress.a

libp_libs:
ifneq (,${verbose})
	${MAKE} -C ${LIBP_LIBS_DIR} $(STRESS_LIBP_LIBS) verbose=${verbose}
else
	@${MAKE} -C ${LIBP_LIBS_DIR} $(STRESS_LIBP_LIBS) --no-print-directory
endif

stressMain:$(CXXOBJS) $(NVCCOBJS) stressMain.o libp_libs
ifneq (,${verbose})
	$(LIBP_LD) -o stressMain stressMain.o $(CXXOBJS) $(NVCCOBJS) $(MESH_OBJS) $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LIBP_LD) -o stressMain stressMain.o $(CXXOBJS) $(NVCCOBJS) $(MESH_OBJS) $(LFLAGS)
endif

libstress.a: $(CXXOBJS) $(NVCCOBJS)
ifneq (,${verbose})
	ar -cr libstress.a $(CXXOBJS) $(NVCCOBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libstress.a $(CXXOBJS) $(NVCCOBJS)
endif

# rule for .cu files
%.o: %.cu $(DEPS) 
ifneq (,${verbose})
	$(LIBP_NVCC) -o $*.o -c $*.cu $(STRESS_NVCCFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_NVCC) -o $*.o -c $*.cu $(STRESS_NVCCFLAGS)
endif

# rule for .cpp files
%.o: %.cpp $(DEPS) 
ifneq (,${verbose})
	$(LIBP_CXX) -o $*.o -c $*.cpp $(STRESS_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $*.o -c $*.cpp $(STRESS_CXXFLAGS)
endif

#cleanup
clean:
	rm -f src/*.o *.o stressMain libstress.a

clean-libs: clean
	${MAKE} -C ${LIBP_LIBS_DIR} clean

clean-kernels: clean-libs
	rm -rf ${LIBP_DIR}/.occa/

realclean: clean
	${MAKE} -C ${LIBP_LIBS_DIR} realclean

help:
	$(info $(value STRESS_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(STRESS_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true

test: stressMain
	@${MAKE} -C $(LIBP_TEST_DIR) --no-print-directory  test-stress
