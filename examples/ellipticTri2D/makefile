ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I../../3rdParty/gslib.github/src -D DHOLMES='"${CURDIR}/../.."'  \
		 	-I../../utilities/parALMOND/include

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L../../utilities/parALMOND -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../utilities/BlasLapack -lBlasLapack -lgfortran

#-llapack -lblas


# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
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
ellipticBuildIpdgTri2D.o \
ellipticBuildContinuousTri2D.o \
ellipticBuildPatchesIpdgTri2D.o 

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
	$(LD)  $(LDFLAGS)  -o ellipticMainTri2D ellipticMainTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

ellipticMainSplitPml2D:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o ellipticMainSplitPml2D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

test2D: $(LOBJS) $(AOBJS) testMainTri2D.o
	$(LD)  $(LDFLAGS)  -o testMain2D testMainTri2D.o $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ellipticMainTri2D.o


