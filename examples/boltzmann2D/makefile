ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++-mpich-mp
LD	= mpic++-mpich-mp

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g

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
boltzmannSplitPmlRun2D.o \
boltzmannSplitPmlSetup2D.o \
boltzmannRampFunction2D.o \
../../src/meshConnect2D.o \
../../src/meshConnectBoundary2D.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshGeometricFactors2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange2D.o \
../../src/meshHaloExtract2D.o \
../../src/meshHaloSetup2D.o \
../../src/meshLoadReferenceNodes2D.o \
../../src/meshParallelConnect2D.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReader2D.o \
../../src/meshPartitionStatistics2D.o \
../../src/meshPhysicalNodes2D.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshReader2D.o \
../../src/meshSetup2D.o \
../../src/meshSurfaceGeometricFactors2D.o \
../../src/meshOccaSetup2D.o \
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o

boltzmannSplitPmlMain2D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o boltzmannSplitPmlMain2D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)


