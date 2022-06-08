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

define INS_HELP_MSG

Incompressible Navier-Stokes solver makefile targets:

   make insMain (default)
   make lib
   make clean
   make clean-libs
   make clean-kernels
   make realclean
   make info
   make help
   make test

Usage:

make insMain
   Build insMain executable.
make lib
   Build libins.a solver library.
make clean
   Clean the insMain executable, library, and object files.
make clean-libs
   In addition to "make clean", also clean the mesh and ogs libraries.
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

ifeq (,$(filter insMain lib clean clean-libs clean-kernels \
                realclean info help test,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${INS_HELP_MSG})
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
ELLIPTIC_DIR =${LIBP_DIR}/solvers/elliptic
FPE_LIBP_LIBS=timeStepper linearSolver parAlmond mesh parAdogs ogs linAlg core

#includes
INCLUDES=-I${ELLIPTIC_DIR} \
			${LIBP_INCLUDES} \
			-I.

#defines
DEFINES =${LIBP_DEFINES} \
         -DLIBP_DIR='"${LIBP_DIR}"'

#.cpp compilation flags
INS_CXXFLAGS=${LIBP_CXXFLAGS} ${DEFINES} ${INCLUDES}

#link libraries
LIBS=-L${ELLIPTIC_DIR} -lelliptic \
	  -L${LIBP_LIBS_DIR} $(addprefix -l,$(FPE_LIBP_LIBS)) \
     ${LIBP_LIBS}

#link flags
LFLAGS=${INS_CXXFLAGS} ${LIBS}

#object dependancies
DEPS=$(wildcard *.hpp) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.hpp) \
     $(wildcard $(ELLIPTIC_DIR)/*.hpp)

SRC =$(wildcard src/*.cpp)

OBJS=$(SRC:.cpp=.o)

.PHONY: all lib libp_libs libelliptic clean clean-libs \
		clean-kernels realclean help info

all: insMain

lib: libins.a

libp_libs:
ifneq (,${verbose})
	${MAKE} -C ${LIBP_LIBS_DIR} $(FPE_LIBP_LIBS) verbose=${verbose}
else
	@${MAKE} -C ${LIBP_LIBS_DIR} $(FPE_LIBP_LIBS) --no-print-directory
endif

libelliptic: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${ELLIPTIC_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${ELLIPTIC_DIR} lib --no-print-directory
endif

insMain:$(OBJS) insMain.o libelliptic
ifneq (,${verbose})
	$(LIBP_LD) -o insMain insMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LIBP_LD) -o insMain insMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
endif

libins.a: $(OBJS)
ifneq (,${verbose})
	ar -cr libins.a $(OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libins.a $(OBJS)
endif

# rule for .cpp files
%.o: %.cpp $(DEPS) | libelliptic
ifneq (,${verbose})
	$(LIBP_CXX) -o $*.o -c $*.cpp $(INS_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $*.o -c $*.cpp $(INS_CXXFLAGS)
endif

#cleanup
clean:
	rm -f src/*.o *.o insMain libins.a

clean-libs: clean
	${MAKE} -C ${ELLIPTIC_DIR} clean
	${MAKE} -C ${LIBP_LIBS_DIR} clean

clean-kernels: clean-libs
	rm -rf ${LIBP_DIR}/.occa/

realclean: clean
	${MAKE} -C ${ELLIPTIC_DIR} clean
	${MAKE} -C ${LIBP_LIBS_DIR} realclean

help:
	$(info $(value INS_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(INS_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true

test: insMain
	@${MAKE} -C $(LIBP_TEST_DIR) --no-print-directory  test-ins
