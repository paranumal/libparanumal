#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
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

define ADVECTION_HELP_MSG

Advection solver makefile targets:

   make advectionMain (default)
   make lib
   make clean
   make clean-libs
   make clean-kernels
   make realclean
   make info
   make help

Usage:

make advectionMain
   Build advectionMain executable.
make lib
   Build libadvection.a solver library.
make clean
   Clean the advectionMain executable, library, and object files.
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

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter advectionMain lib clean clean-libs clean-kernels \
                realclean info help,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${ADVECTION_HELP_MSG})
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
GS_DIR       =${LIBP_TPL_DIR}/gslib
BLAS_DIR     =${LIBP_TPL_DIR}/BlasLapack
OGS_DIR      =${LIBP_LIBS_DIR}/ogs
MESH_DIR     =${LIBP_LIBS_DIR}/mesh
CORE_DIR     =${LIBP_DIR}/core

#includes
INCLUDES=${LIBP_INCLUDES} \
		 -I${LIBP_LIBS_DIR}/include \
         -I.


#defines
DEFINES =${LIBP_DEFINES} \
         -DLIBP_DIR='"${LIBP_DIR}"'

#.cpp compilation flags
ADVECTION_CXXFLAGS=${CXXFLAGS} ${DEFINES} ${INCLUDES}

#link libraries
LIBS=-L${LIBP_LIBS_DIR} -lmesh -logs\
     -L$(GS_DIR)/lib -lgs \
     -L${CORE_DIR} -lcore \
     -L${BLAS_DIR} -lBlasLapack \
     -lgfortran \
     ${LIBP_LIBS}

#link flags
LFLAGS=${ADVECTION_CXXFLAGS} ${LIBS}

#object dependancies
DEPS=$(wildcard *.hpp) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.hpp)

SRC =$(wildcard src/*.cpp)

OBJS=$(SRC:.cpp=.o)

.PHONY: all lib libcore libmesh libogs libblas clean clean-libs \
		clean-kernels realclean help info advectionMain

all: advectionMain

lib: libadvection.a

advectionMain:$(OBJS) advectionMain.o | libmesh
ifneq (,${verbose})
	$(LD) -o advectionMain advectionMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LD) -o advectionMain advectionMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
endif

libadvection.a: $(OBJS)
ifneq (,${verbose})
	ar -cr libadvection.a $(OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libadvection.a $(OBJS)
endif

libmesh: libogs libgs libblas libcore
ifneq (,${verbose})
	${MAKE} -C ${MESH_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${MESH_DIR} lib --no-print-directory
endif

libogs: libcore
ifneq (,${verbose})
	${MAKE} -C ${OGS_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${OGS_DIR} lib --no-print-directory
endif

libcore: libgs
ifneq (,${verbose})
	${MAKE} -C ${CORE_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${CORE_DIR} lib --no-print-directory
endif

libgs: libblas
ifneq (,${verbose})
	${MAKE} -C $(GS_DIR) install verbose=${verbose}
else
	@${MAKE} -C $(GS_DIR) install --no-print-directory
endif

libblas:
ifneq (,${verbose})
	${MAKE} -C ${BLAS_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${BLAS_DIR} lib --no-print-directory
endif

# rule for .cpp files
%.o: %.cpp $(DEPS) | libmesh
ifneq (,${verbose})
	$(CXX) -o $*.o -c $*.cpp $(ADVECTION_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(CXX) -o $*.o -c $*.cpp $(ADVECTION_CXXFLAGS)
endif

#cleanup
clean:
	rm -f src/*.o *.o advectionMain libadvection.a

clean-libs: clean
	${MAKE} -C ${OGS_DIR} clean
	${MAKE} -C ${MESH_DIR} clean
	${MAKE} -C ${CORE_DIR} clean

clean-kernels: clean-libs
# 	$(shell ${OCCA_DIR}/bin/occa clear all -y)
	rm -rf ~/.occa/

realclean: clean-libs
	${MAKE} -C ${GS_DIR} clean
	${MAKE} -C ${BLAS_DIR} clean

help:
	$(info $(value ADVECTION_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(ADVECTION_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true
