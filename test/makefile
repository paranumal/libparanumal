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

define TEST_HELP_MSG

Testing makefile targets:

	 make test (default)
	 make info
	 make help

Usage:

make test
	 Run tests.
make info
	 List directories and compiler flags in use.
make help
	 Display this help message.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter info help test,$(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${TEST_HELP_MSG})
endif
endif

ifndef LIBP_MAKETOP_LOADED
ifeq (,$(wildcard ../make.top))
$(error cannot locate ${PWD}/../make.top)
else
include ../make.top
endif
endif

#libraries
GS_DIR       =${LIBP_TPL_DIR}/gslib
BLAS_DIR     =${LIBP_TPL_DIR}/BlasLapack
OGS_DIR      =${LIBP_LIBS_DIR}/ogs
MESH_DIR     =${LIBP_LIBS_DIR}/mesh
CORE_DIR     =${LIBP_DIR}/core
TEST_DIR     =${LIBP_DIR}/test

.PHONY: all help info test

all: test

help:
	$(info $(value TEST_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info LIBS      = $(LIBS))
	@true

test:
	@./testMesh.py
	@./testGradient.py
	@./testAdvection.py
	@./testAcoustics.py
	@./testElliptic.py
	@./testFokkerPlanck.py
	@./testCns.py
	@./testBns.py
	@./testIns.py
	@./testTimeStepper.py
	@./testLinearSolver.py
	@./testParAlmond.py
