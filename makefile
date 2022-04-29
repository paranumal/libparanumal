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

define LIBP_HELP_MSG

LibParanumal makefile targets:

	 make solvers (default)
	 make {solver}
	 make clean
	 make clean-kernels
	 make realclean
	 make info
	 make help
	 make test

Usage:

make solvers
	 Builds each solver executable.
make {solver}
	 Builds a solver executable,
	 solver can be acoustics/advection/bns/cns/elliptic/fokkerPlanck/gradient/ins.
make clean
	 Cleans all solver executables, libraries, and object files.
make clean-{solver}
	 Cleans a solve executable, library, and object files.
make clean-kernels
	 In addition to "make clean", also cleans the cached OCCA kernels.
make realclean
	 In addition to "make clean", also clean 3rd party libraries.
make info
	 List directories and compiler flags in use.
make help
	 Display this help message.
make test
	 Run the included solver examples.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter solvers \
				acoustics advection bns cns elliptic fokkerPlanck gradient ins \
				lib clean clean-kernels \
				realclean info help test,$(MAKECMDGOALS)))
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

#gslib
GS_DIR=${LIBP_TPL_DIR}/gslib

#libraries
LIBP_CORE_LIBS=timeStepper linearSolver parAlmond mesh ogs linAlg core
SOLVER_DIR   =${LIBP_DIR}/solvers

.PHONY: all solvers libp_libs \
			acoustics advection bns lbs cns elliptic fokkerPlanck gradient ins \
			clean clean-libs realclean help info

all: solvers

solvers: acoustics advection bns lbs cns elliptic fokkerPlanck gradient ins

libp_libs:
ifneq (,${verbose})
	${MAKE} -C ${LIBP_LIBS_DIR} $(LIBP_CORE_LIBS) verbose=${verbose}
else
	@${MAKE} -C ${LIBP_LIBS_DIR} $(LIBP_CORE_LIBS) --no-print-directory
endif

acoustics: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

advection: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

bns: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

lbs: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

cns: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

elliptic: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

fokkerPlanck: libp_libs | elliptic
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

gradient: libp_libs
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

ins: libp_libs | elliptic
ifneq (,${verbose})
	${MAKE} -C ${SOLVER_DIR}/$(@F) verbose=${verbose}
else
	@printf "%b" "$(SOL_COLOR)Building $(@F) solver$(NO_COLOR)\n";
	@${MAKE} -C ${SOLVER_DIR}/$(@F) --no-print-directory
endif

#cleanup
clean: clean-acoustics clean-advection clean-bns clean-lbs clean-cns \
	   clean-elliptic clean-fokkerPlanck clean-gradient clean-ins \
	   clean-libs

clean-acoustics:
	${MAKE} -C ${SOLVER_DIR}/acoustics clean

clean-advection:
	${MAKE} -C ${SOLVER_DIR}/advection clean

clean-bns:
	${MAKE} -C ${SOLVER_DIR}/bns clean

clean-lbs:
	${MAKE} -C ${SOLVER_DIR}/lbs clean

clean-cns:
	${MAKE} -C ${SOLVER_DIR}/cns clean

clean-elliptic:
	${MAKE} -C ${SOLVER_DIR}/elliptic clean

clean-fokkerPlanck:
	${MAKE} -C ${SOLVER_DIR}/fokkerPlanck clean

clean-gradient:
	${MAKE} -C ${SOLVER_DIR}/gradient clean

clean-ins:
	${MAKE} -C ${SOLVER_DIR}/ins clean

clean-libs:
	${MAKE} -C ${LIBP_LIBS_DIR} clean

clean-kernels: clean
# 	$(shell ${OCCA_DIR}/bin/occa clear all -y)
	rm -rf ~/.occa/

realclean: clean
	${MAKE} -C ${LIBP_LIBS_DIR} realclean

help:
	$(info $(value LIBP_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(CXXFLAGS))
	@true

test: all
	@${MAKE} -C $(LIBP_TEST_DIR) --no-print-directory  test
