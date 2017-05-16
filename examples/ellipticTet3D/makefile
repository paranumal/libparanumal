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
LIBS	=  $(links) -llapack -lblas

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
ellipticSetupTet3D.o \
ellipticPreconditionerSetupTet3D.o\
ellipticCoarsePreconditionerSetupTet3D.o\
ellipticSolveTet3D.o\
ellipticSolveSetupTet3D.o\
ellipticHaloExchange3D.o\
ellipticParallelGatherScatterTet3D.o

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
../../src/meshParallelReaderTet3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/meshSetupTet3D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshLoadReferenceNodesTet3D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshParallelGather.o \
../../src/meshParallelScatter.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherSetup.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup3D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
./diagnostic.o

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

ellipticMainTet3D:$(AOBJS) $(LOBJS) ellipticMainTet3D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o ellipticMainTet3D ellipticMainTet3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

ellipticMainSplitPml3D:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o ellipticMainSplitPml3D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

test3D: $(LOBJS) $(AOBJS) testMainTet3D.o
	$(LD)  $(LDFLAGS)  -o testMain3D testMainTet3D.o $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ellipticMainTet3D.o


