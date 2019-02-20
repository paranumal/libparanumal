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

define LIBP_HELP_MSG

LibParanumal makefile targets:

	 make solvers (default)
	 make {solver}
	 make clean
	 make clean-libs
	 make realclean
	 make info
	 make help

Usage:

make solvers
	 Builds each solver executable.
make {solver}
	 Builds a solver executable,
	 solver can be acoustics/advection/bns/cns/elliptic/gradient/ins.
make clean
	 Cleans all solver executables, libraries, and object files.
make clean-{solver}
	 Cleans a solve executable, library, and object files.
make clean-libs
	 In addition to "make clean", also clean the ogs and parAlmond libraries.
make realclean
	 In addition to "make clean-libs", also clean 3rd party libraries.
make info
	 List directories and compiler flags in use.
make help
	 Display this help message.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter solvers \
								acoustics advection bns cns elliptic gradient ins \
								lib clean clean-libs \
								realclean info help,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${LIBP_HELP_MSG})
endif
endif

ifndef LIBP_MAKETOP_LOADED
ifeq (,$(wildcard make.top))
$(error cannot locate ${PWD}/make.top)
else
include make.top
endif
endif

#libraries
GS_DIR       =${LIBP_TPL_DIR}/gslib
BLAS_DIR     =${LIBP_TPL_DIR}/BlasLapack
OGS_DIR      =${LIBP_LIBS_DIR}/gatherScatter
PARALMOND_DIR=${LIBP_LIBS_DIR}/parAlmond
MESH_DIR     =${LIBP_DIR}/src
SOLVER_DIR   =${LIBP_DIR}/solvers

.PHONY: all solvers \
				acoustics advection bns cns elliptic gradient ins \
				libparAlmond libmesh libogs ligs libblas \
				clean clean-libs realclean help info

all: solvers

solvers: acoustics advection bns cns elliptic gradient ins

acoustics: libmesh
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

advection: libmesh
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

bns: libmesh
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

cns: libmesh
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

elliptic: libparAlmond
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

gradient: libmesh
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

ins: elliptic
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

libmesh: libogs libgs libblas
ifneq (,${verbose})
	${MAKE} -C ${MESH_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${MESH_DIR} lib --no-print-directory
endif

libparAlmond: libmesh
ifneq (,${verbose})
	${MAKE} -C ${PARALMOND_DIR} lib verbose=${verbose}
else
	@${MAKE} -C ${PARALMOND_DIR} lib --no-print-directory
endif

libogs: libgs libblas
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

#cleanup
clean: clean-acoustics clean-advection clean-bns clean-cns \
			 clean-elliptic clean-gradient clean-ins

clean-acoustics:
	${MAKE} -C ${SOLVER_DIR}/acoustics clean

clean-advection:
	${MAKE} -C ${SOLVER_DIR}/advection clean

clean-bns:
	${MAKE} -C ${SOLVER_DIR}/bns clean

clean-cns:
	${MAKE} -C ${SOLVER_DIR}/cns clean

clean-elliptic:
	${MAKE} -C ${SOLVER_DIR}/elliptic clean

clean-gradient:
	${MAKE} -C ${SOLVER_DIR}/gradient clean

clean-ins:
	${MAKE} -C ${SOLVER_DIR}/ins clean

clean-libs: clean
	${MAKE} -C ${MESH_DIR} clean
	${MAKE} -C ${OGS_DIR} clean
	${MAKE} -C ${PARALMOND_DIR} clean

realclean: clean-libs
	${MAKE} -C ${GS_DIR} clean
	${MAKE} -C ${BLAS_DIR} clean

help:
	$(info $(value LIBP_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(CXXFLAGS))
	@true
