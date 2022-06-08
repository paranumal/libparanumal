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

define CNS_HELP_MSG

Compressible Navier-Stokes flow solver makefile targets:

   make cnsMain (default)
   make lib
   make clean
   make clean-libs
   make clean-kernels
   make realclean
   make info
   make help
   make test

Usage:

make cnsMain
   Build cnsMain executable.
make lib
   Build libcns.a solver library.
make clean
   Clean the cnsMain executable, library, and object files.
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

ifeq (,$(filter cnsMain lib clean clean-libs clean-kernels \
                realclean info help test,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${CNS_HELP_MSG})
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
CNS_LIBP_LIBS=timeStepper mesh parAdogs ogs linAlg core

#includes
INCLUDES=${LIBP_INCLUDES} \
				 -I.

#defines
DEFINES =${LIBP_DEFINES} \
         -DLIBP_DIR='"${LIBP_DIR}"'

#.cpp compilation flags
CNS_CXXFLAGS=${LIBP_CXXFLAGS} ${DEFINES} ${INCLUDES}

#link libraries
LIBS=-L${LIBP_LIBS_DIR} $(addprefix -l,$(CNS_LIBP_LIBS)) \
     ${LIBP_LIBS}

#link flags
LFLAGS=${CNS_CXXFLAGS} ${LIBS}

#object dependancies
DEPS=$(wildcard *.hpp) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.hpp)

SRC =$(wildcard src/*.cpp)

OBJS=$(SRC:.cpp=.o)

.PHONY: all lib libp_libs clean clean-libs \
		clean-kernels realclean help info

all: cnsMain

lib: libcns.a

libp_libs:
ifneq (,${verbose})
	${MAKE} -C ${LIBP_LIBS_DIR} $(CNS_LIBP_LIBS) verbose=${verbose}
else
	@${MAKE} -C ${LIBP_LIBS_DIR} $(CNS_LIBP_LIBS) --no-print-directory
endif

cnsMain:$(OBJS) cnsMain.o libp_libs
ifneq (,${verbose})
	$(LIBP_LD) -o cnsMain cnsMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LIBP_LD) -o cnsMain cnsMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
endif

libcns.a: $(OBJS)
ifneq (,${verbose})
	ar -cr libcns.a $(OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libcns.a $(OBJS)
endif

# rule for .cpp files
%.o: %.cpp $(DEPS) | libp_libs
ifneq (,${verbose})
	$(LIBP_CXX) -o $*.o -c $*.cpp $(CNS_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $*.o -c $*.cpp $(CNS_CXXFLAGS)
endif

#cleanup
clean:
	rm -f src/*.o *.o cnsMain libcns.a

clean-libs: clean
	${MAKE} -C ${LIBP_LIBS_DIR} clean

clean-kernels: clean-libs
	rm -rf ${LIBP_DIR}/.occa/

realclean: clean
	${MAKE} -C ${LIBP_LIBS_DIR} realclean

help:
	$(info $(value CNS_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(CNS_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true

test: cnsMain
	@${MAKE} -C $(LIBP_TEST_DIR) --no-print-directory  test-cns
