
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib.github/src
ALMONDDIR = ../../examples/parALMOND
ELLIPTICDIR = ../../examples/ellipticTri2D

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -I$(ELLIPTICDIR) -g  -D LIBP_DIR='"${CURDIR}/../.."'

# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND -L$(ELLIPTICDIR) -lellipticTri2D

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

INCLUDES = ins2D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.hpp \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/ogs_t.h \
$(HDRDIR)/hgs_t.h \
$(ALMONDDIR)/parALMOND.h \
$(ELLIPTICDIR)/ellipticTri2D.h \
$(ELLIPTICDIR)/ellipticPreconTri2D.h

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
insSetup2D.o \
insRunBenchmark2D.o

# library objects
LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderTri2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/meshSetupTri2D.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshParallelGather.o \
../../src/meshParallelScatter.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherSetup.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup2D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
../../src/timer.o\

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

insBenchmarkTri2D:$(AOBJS) $(LOBJS) insBenchmarkTri2D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../examples/parALMOND; make -j lib; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../examples/ellipticTri2D; make -j lib; cd ../../benchmarks/insBenchmarkTri2D
	$(LD)  $(LDFLAGS)  -o insBenchmarkTri2D insBenchmarkTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../benchmarks/insBenchmarkTri2D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib insBenchmarkTri2D

# what to do if user types "make clean"
clean:
	cd ../../examples/parALMOND; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../examples/ellipticTri2D; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../src; rm *.o; cd ../benchmarks/insBenchmarkTri2D
	rm *.o insBenchmarkTri2D libinsTri2D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../examples/parALMOND; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../examples/ellipticTri2D; make clean; cd ../../benchmarks/insBenchmarkTri2D
	cd ../../src; rm *.o; cd ../benchmarks/insBenchmarkTri2D
	rm *.o insBenchmarkTri2D