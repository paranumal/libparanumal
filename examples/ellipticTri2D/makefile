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
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -g  -D DHOLMES='"${CURDIR}/../.."'


# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

#-llapack -lblas

INCLUDES = ellipticTri2D.h ellipticPreconTri2D.h
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
PCG.o \
PGMRES.o \
PBiCGStab.o \
ellipticSetupTri2D.o \
ellipticPreconditionerSetupTri2D.o\
ellipticPreconditionerTri2D.o\
ellipticSolveTri2D.o\
ellipticSolveSetupTri2D.o\
ellipticHaloExchange2D.o\
ellipticParallelGatherScatterTri2D.o \
ellipticParallelGatherScatterSetup.o \
ellipticBuildIpdgTri2D.o \
ellipticBuildBRdgTri2D.o \
ellipticBuildContinuousTri2D.o \
ellipticBuildFullPatchesTri2D.o \
ellipticBuildFacePatchesTri2D.o \
ellipticBuildJacobiTri2D.o \
ellipticBuildLocalPatchesTri2D.o \
ellipticBuildMultigridLevelTri2D.o \
ellipticMultiGridSetupTri2D.o \
ellipticSEMFEMSetupTri2D.o\
ellipticSetupSmoothersTri2D.o \
ellipticSmoothersTri2D.o

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
../../src/meshParallelConsecutiveGlobalNumbering.o\
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
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/ellipticTri2D
	cd ../parALMOND; make -j lib; cd ../ellipticTri2D
	$(LD)  $(LDFLAGS)  -o ellipticMainTri2D ellipticMainTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libellipticTri2D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/ellipticTri2D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib ellipticMainTri2D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../ellipticTri2D
	cd ../../src; rm *.o; cd ../examples/ellipticTri2D
	rm *.o ellipticMainTri2D libellipticTri2D.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/ellipticTri2D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/ellipticTri2D
	cd ../parALMOND; make clean; cd ../ellipticTri2D
	cd ../../src; rm *.o; cd ../examples/ellipticTri2D
	rm *.o ellipticMainTri2D libellipticTri2D.a

