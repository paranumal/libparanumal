ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
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
OBJS    = \
boltzmannSplitPmlMain2D.o \
boltzmannComputeVorticity2D.o \
boltzmannError2D.o \
boltzmannCouetteError2D.o \
boltzmannPeriodic2D.o\
boltzmannSplitPmlRun2D.o \
boltzmannSplitPmlSetup2D.o \
boltzmannSplitPmlLserkStep2D.o \
boltzmannSplitPmlMrabStep2D.o \
boltzmannSplitPmlSaabStep2D.o \
boltzmannSplitPmlLsimexStep2D.o \
boltzmannSplitPmlSarkStep2D.o \
boltzmannSplitPmlSark3Step2D.o \
boltzmannSplitPmlLserk3Step2D.o \
boltzmannRampFunction2D.o \
boltzmannReport2D.o \
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
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshReader2D.o \
../../src/meshOccaSetup2D.o\
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o

boltzmannSplitPmlMain2D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o boltzmannSplitPmlMain2D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS) boltzmannSplitPmlMain2D
 

