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

INCLUDES = ellipticQuad2D.h ellipticPreconQuad2D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
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
ellipticSetupQuad2D.o \
ellipticSolveSetupQuad2D.o \
ellipticSolveQuad2D.o \
ellipticHaloExchange2D.o \
ellipticParallelGatherScatterQuad2D.o \
ellipticParallelGatherScatterSetup.o \
ellipticPreconditionerQuad2D.o\
ellipticPreconditionerSetupQuad2D.o\
ellipticBuildIpdgQuad2D.o \
ellipticBuildContinuousQuad2D.o \
ellipticBuildFullPatchesQuad2D.o \
ellipticBuildFacePatchesQuad2D.o \
ellipticBuildJacobiQuad2D.o \
ellipticBuildLocalPatchesQuad2D.o \
ellipticBuildMultigridLevelQuad2D.o \
ellipticMultiGridSetupQuad2D.o \
ellipticSEMFEMSetupQuad2D.o\
ellipticSetupSmoothersQuad2D.o \
ellipticSmoothersQuad2D.o


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
../../src/meshParallelReaderQuad2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConsecutiveGlobalNumbering.o\
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshOccaSetup2D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
../../src/timer.o\
../../src/readArray.o\
./diagnostic.o\


COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o 

ellipticMainQuad2D:$(AOBJS) $(LOBJS) ellipticMainQuad2D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/ellipticQuad2D
	cd ../parALMOND; make -j lib; cd ../ellipticQuad2D
	$(LD)  $(LDFLAGS)  -o ellipticMainQuad2D ellipticMainQuad2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libellipticQuad2D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/ellipticQuad2D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib ellipticMainQuad2D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../ellipticQuad2D
	cd ../../src; rm *.o; cd ../examples/ellipticQuad2D
	rm *.o ellipticMainQuad2D libellipticQuad2D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/ellipticQuad2D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/ellipticQuad2D
	cd ../parALMOND; make clean; cd ../ellipticQuad2D
	cd ../../src; rm *.o; cd ../examples/ellipticQuad2D
	rm *.o ellipticMainQuad2D libellipticQuad2D.a

