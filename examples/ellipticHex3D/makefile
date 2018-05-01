ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR = ../../include
GSDIR  = ../../3rdParty/gslib.github/src
ALMONDDIR = ../parALMOND

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc.openmpi
CC	= mpic++.openmpi
LD	= mpic++.openmpi

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -g  -D DHOLMES='"${CURDIR}/../.."'


# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

#-llapack -lblas

INCLUDES = ellipticHex3D.h ellipticPreconHex3D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h \
$(HDRDIR)/ogs_t.h \
$(ALMONDDIR)/parALMOND.h \

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
PCG.o \
ellipticSetupHex3D.o \
ellipticSolveSetupHex3D.o \
ellipticSolveHex3D.o \
ellipticHaloExchange3D.o \
ellipticParallelGatherScatterHex3D.o \
ellipticParallelGatherScatterSetup.o \
ellipticPreconditionerHex3D.o\
ellipticPreconditionerSetupHex3D.o\
ellipticBuildIpdgHex3D.o \
ellipticBuildContinuousHex3D.o \
ellipticBuildJacobiHex3D.o \
ellipticBuildLocalPatchesHex3D.o \
ellipticBuildMultigridLevelHex3D.o \
ellipticMultiGridSetupHex3D.o \
ellipticSEMFEMSetupHex3D.o\
ellipticSetupSmoothersHex3D.o \
ellipticSmoothersHex3D.o


# library objects
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
../../src/meshParallelConsecutiveGlobalNumbering.o\
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/meshSetupHex3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshOccaSetup3D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
../../src/timer.o\
../../src/readArray.o


COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o 

ellipticMainHex3D:$(AOBJS) $(LOBJS) ellipticMainHex3D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/ellipticHex3D
	cd ../parALMOND; make -j lib; cd ../ellipticHex3D
	$(LD)  $(LDFLAGS)  -o ellipticMainHex3D ellipticMainHex3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libellipticHex3D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/ellipticHex3D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib ellipticMainHex3D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../ellipticHex3D
	cd ../../src; rm *.o; cd ../examples/ellipticHex3D
	rm *.o ellipticMainHex3D libellipticHex3D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/ellipticHex3D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/ellipticHex3D
	cd ../parALMOND; make clean; cd ../ellipticHex3D
	cd ../../src; rm *.o; cd ../examples/ellipticHex3D
	rm *.o ellipticMainHex3D libellipticHex3D.a
