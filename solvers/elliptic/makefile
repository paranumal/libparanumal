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

define ELLIPTIC_HELP_MSG

Elliptic solver makefile targets:

   make ellipticMain (default)
   make lib
   make clean
   make clean-libs
   make realclean
   make info
   make help

Usage:

make ellipticMain
   Build ellipticMain executable.
make lib
   Build libelliptic.a solver library.
make clean
   Clean the ellipticMain executable, library, and object files.
make clean-libs
   In addition to "make clean", also clean the mesh, ogs, and parAlmond libraries.
make realclean
   In addition to "make clean-libs", also clean 3rd party libraries.
make info
   List directories and compiler flags in use.
make help
   Display this help message.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter ellipticMain lib clean clean-libs \
                realclean info help,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${ELLIPTIC_HELP_MSG})
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
OGS_DIR      =${LIBP_LIBS_DIR}/gatherScatter
PARALMOND_DIR=${LIBP_LIBS_DIR}/parAlmond
MESH_DIR     =${LIBP_DIR}/src

#includes
ELLIPTIC_INCLUDES=-I. \
         -I$(OGS_DIR) \
         -I$(PARALMOND_DIR) \
         ${LIBP_INCLUDES}

#defines
ELLIPTIC_DEFINES =${LIBP_DEFINES} \
         -DDHOLMES='"${LIBP_DIR}"' \
         -DDELLIPTIC='"${CURDIR}"'

#.cpp compilation flags
ELLIPTIC_CXXFLAGS=${CXXFLAGS} ${ELLIPTIC_DEFINES} ${ELLIPTIC_INCLUDES}

#link libraries
LIBS=-L$(MESH_DIR) -lmesh \
     -L$(PARALMOND_DIR) -lparAlmond  \
     -L$(OGS_DIR) -logs \
     -L$(GS_DIR)/lib -lgs \
     -L${BLAS_DIR} -lBlasLapack \
     -lgfortran \
     ${LIBP_LIBS}

#link flags
LFLAGS=${ELLIPTIC_CXXFLAGS} ${LIBS}

#object dependancies
DEPS=$(wildcard *.h) \
     $(wildcard $(LIBP_INCLUDE_DIR)/*.h) \
     $(OGS_DIR)/ogs.hpp \
     $(PARALMOND_DIR)/parAlmond.hpp

SRC =$(wildcard src/*.c)
OBJS=$(SRC:.c=.o)

.PHONY: all lib libmesh libogs libblas libparAlmond clean clean-libs realclean help info

all: ellipticMain

lib: libelliptic.a silentUpdate

ellipticMain:$(OBJS) $(MESH_OBJS) ellipticMain.o | libparAlmond
ifneq (,${verbose})
	$(LD) -o ellipticMain ellipticMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
else
	@printf "%b" "$(EXE_COLOR)Linking $(@F)$(NO_COLOR)\n";
	@$(LD) -o ellipticMain ellipticMain.o $(OBJS) $(MESH_OBJS) $(LFLAGS)
endif

libelliptic.a: $(OBJS)
ifneq (,${verbose})
	ar -cr libelliptic.a $(OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libelliptic.a $(OBJS)
endif

libmesh: libogs libgs libblas
ifneq (,${verbose})
	${MAKE} -C ${MESH_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${MESH_DIR} lib --no-print-directory
endif

libogs: libgs
ifneq (,${verbose})
	${MAKE} -C ${OGS_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${OGS_DIR} lib --no-print-directory
endif

libgs:
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

libparAlmond: | libmesh
ifneq (,${verbose})
	${MAKE} -C ${PARALMOND_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${PARALMOND_DIR} lib --no-print-directory
endif

# rule for .c files
%.o: %.c $(DEPS) | libparAlmond
ifneq (,${verbose})
	$(CXX) -o $*.o -c $*.c $(ELLIPTIC_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(CXX) -o $*.o -c $*.c $(ELLIPTIC_CXXFLAGS)
endif

#cleanup
clean:
	rm -f src/*.o ellipticMain libelliptic.a

clean-libs: clean
	${MAKE} -C ${OGS_DIR} clean
	${MAKE} -C ${PARALMOND_DIR} clean
	${MAKE} -C ${MESH_DIR} clean

realclean: clean-libs
	${MAKE} -C ${GS_DIR} clean
	${MAKE} -C ${BLAS_DIR} clean

help:
	$(info $(value ELLIPTIC_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(ELLIPTIC_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true

silentUpdate:
	@true
