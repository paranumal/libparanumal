ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linkingcc	
cc 	= mpicc
CC	= mpic++
LD	= mpic++


# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g  -D DHOLMES='"${CURDIR}/../.."'

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g

# libraries to be linked in
LIBS	=  $(links)

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
AOBJS    = \
boltzmannMain2D.o \
boltzmannComputeVorticity2D.o \
boltzmannError2D.o \
boltzmannCouetteError2D.o \
boltzmannPeriodic2D.o\
boltzmannRun2D.o \
boltzmannSetup2D.o \
boltzmannSRABStep2D.o \
boltzmannLSERKStep2D.o \
boltzmannLSIMEXStep2D.o \
boltzmannMRABStep2D.o \
boltzmannMRSAABStep2D.o \
boltzmannSAABStep2D.o \
boltzmannSARKStep2D.o \
boltzmannRampFunction2D.o \
boltzmannReport2D.o \
boltzmannPlotVTU2D.o \
boltzmannPlotTEC2D.o \
boltzmannMRABSetup2D.o \
boltzmannMRABPmlSetup2D.o \
boltzmannPmlSetup2D.o \
boltzmannMRRun2D.o \
boltzmannSACoefs2D.o \

# library objects
LOBJS = \
../../src/meshParallelReaderTri2D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshSetupTri2D.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelPrint2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshReader2D.o \
../../src/meshMRABSetup2D.o \
../../src/meshBuildMRABClusters2D.o \
../../src/meshClusteredGeometricPartition2D.o \
../../src/meshMRABWeightedPartitionTri2D.o \
../../src/meshOccaSetup2D.o\
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o

boltzmannMain2D:$(AOBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o boltzmannMain2D $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) boltzmannMain2D
 

