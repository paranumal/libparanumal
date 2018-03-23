
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib.github/src
ALMONDDIR = ../parALMOND
ELLIPTICDIR = ../ellipticTri2D

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc.openmpi
CC	= mpic++.openmpi
LD	= mpic++.openmpi

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -I$(ELLIPTICDIR) -g  -D DHOLMES='"${CURDIR}/../.."' 
		
# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ELLIPTICDIR) -lellipticTri2D -L$(ALMONDDIR) -lparALMOND 

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

INCLUDES = insTri2D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/ogs_t.h \
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
insSetupTri2D.o \
insPlotVTUTri2D.o \
insErrorTri2D.o \
insErrorNormsTri2D.o \
insReportTri2D.o \
insRunTri2D.o \
insAdvectionStepTri2D.o \
insAdvectionSubCycleStepTri2D.o \
insHelmholtzStepTri2D.o \
insPoissonStepTri2D.o \
insUpdateStepTri2D.o 

# library objects
LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshConnectFaceModes2D.o \
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
../../src/meshParallelGatherScatter.o \
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

insMainTri2D:$(AOBJS) $(LOBJS) insMainTri2D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/insTri2D
	cd ../parALMOND; make -j lib; cd ../insTri2D
	cd ../ellipticTri2D; make -j lib; cd ../insTri2D
	$(LD)  $(LDFLAGS)  -o insMainTri2D insMainTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libinsTri2D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/insTri2D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib insMainTri2D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../insTri2D
	cd ../ellipticTri2D; make clean; cd ../insTri2D
	cd ../../src; rm *.o; cd ../examples/insTri2D
	rm *.o insMainTri2D libinsTri2D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/insTri2D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/insTri2D
	cd ../parALMOND; make clean; cd ../insTri2D
	cd ../ellipticTri2D; make clean; cd ../insTri2D
	cd ../../src; rm *.o; cd ../examples/insTri2D
	rm *.o insMainTri2D libinsTri2D.a