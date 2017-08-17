ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include
ALMONDDIR = ../parALMOND

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I../../3rdParty/gslib.github/src -D DHOLMES='"${CURDIR}/../.."'  \
		 	-I$(ALMONDDIR)/include

# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../utilities/BlasLapack -lBlasLapack -lgfortran

#-llapack -lblas

INCLUDES = ellipticTri2D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/ogs_t.h \
$(HDRDIR)/hgs_t.h \
$(ALMONDDIR)/parALMOND.h \

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
ellipticSetupTri2D.o \
ellipticPreconditionerSetupTri2D.o\
ellipticCoarsePreconditionerSetupTri2D.o\
ellipticPreconditionerTri2D.o\
ellipticSolveTri2D.o\
ellipticSolveSetupTri2D.o\
ellipticHaloExchange2D.o\
ellipticParallelGatherScatterTri2D.o \
ellipticParallelGatherScatterSetup.o \
ellipticBuildIpdgTri2D.o \
ellipticBuildContinuousTri2D.o \
ellipticBuildFullPatchesIpdgTri2D.o \
ellipticBuildFacePatchesIpdgTri2D.o \
ellipticBuildJacobiIpdgTri2D.o \
ellipticBuildLocalPatchesIpdgTri2D.o \
ellipticBuildMultigridLevelTri2D.o \
ellipticMultiGridSetupTri2D.o \
ellipticSetupSmoothersTri2D.o \
ellipticSmoothersTri2D.o

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
./diagnostic.o\

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

ellipticMainTri2D:$(AOBJS) $(LOBJS) ellipticMainTri2D.o gslibInterface
	cd ../parALMOND; make -j; cd ../ellipticTri2D
	$(LD)  $(LDFLAGS)  -o ellipticMainTri2D ellipticMainTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j; cd ../../examples/ellipticTri2D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

# what to do if user types "make clean"
clean :
	cd ../parALMOND; make clean; cd ../ellipticTri2D
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/ellipticTri2D
	cd ../../src; rm *.o; cd ../examples/ellipticTri2D
	rm *.o ellipticMainTri2D


