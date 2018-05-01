
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
ELLIPTICDIR = ../ellipticHex3D

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc  = mpicc.openmpi
CC  = mpic++.openmpi
LD  = mpic++.openmpi

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I$(GSDIR) -I$(ELLIPTICDIR) -D DHOLMES='"${CURDIR}/../.."'  

# link flags to be used 
LDFLAGS = $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ELLIPTICDIR) -lellipticHex3D -L$(ALMONDDIR) -lparALMOND

# libraries to be linked in
LIBS  =  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran


INCLUDES = insHex3D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h \
$(HDRDIR)/ogs_t.h \
$(ALMONDDIR)/parALMOND.h \
$(ELLIPTICDIR)/ellipticHex3D.h \
$(ELLIPTICDIR)/ellipticPreconHex3D.h  

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
insSetupHex3D.o \
insRunHex3D.o \
insErrorHex3D.o \
insReportHex3D.o \
insPlotVTUHex3D.o \
insPlotSliceHex3D.o \
insPlotContourHex3D.o \
insAdvectionStepHex3D.o \
insAdvectionSubCycleStepHex3D.o \
insHelmholtzStepHex3D.o \
insPoissonStepHex3D.o \
insUpdateStepHex3D.o 

LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU3D.o \
../../src/meshPlotContour3D.o \
../../src/meshPlotAdaptiveContour3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/meshSetupHex3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup3D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o \
../../src/timer.o

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

insMainHex3D:$(AOBJS) $(LOBJS) insMainHex3D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/insHex3D
	cd ../parALMOND; make -j lib; cd ../insHex3D
	cd ../ellipticHex3D; make -j lib; cd ../insHex3D
	$(LD)  $(LDFLAGS)  -o insMainHex3D insMainHex3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libinsHex3D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/insHex3D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib insMainHex3D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../insHex3D
	cd ../ellipticHex3D; make clean; cd ../insHex3D
	cd ../../src; rm *.o; cd ../examples/insHex3D
	rm *.o insMainHex3D libinsHex3D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/insHex3D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/insHex3D
	cd ../parALMOND; make clean; cd ../insHex3D
	cd ../ellipticHex3D; make clean; cd ../insHex3D
	cd ../../src; rm *.o; cd ../examples/insHex3D
	rm *.o insMainHex3D libinsHex3D.a
