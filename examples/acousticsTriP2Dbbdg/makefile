ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ./.

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++
LD	= mpic++
#CC	= mpic++-mpich-mp
#LD	= mpic++-mpich-mp

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -D DHOLMES='"${CURDIR}/../.."'

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
acousticsGaussianPulse2D.o \
acousticsError2D.o \
acousticsRun2D.o \
acousticsSetup2D.o \
acousticsSurface2D.o \
acousticsUpdate2D.o \
acousticsVolume2D.o \
acousticsMain2D.o \
meshSetupTriP2D.o \
meshConnectFaceNodesP2D.o \
meshPhysicalNodesTriP2D.o \
meshPlotVTU2DP.o \
meshHaloSetupP.o \
meshSurfaceGeometricFactorsTriP2D.o \
meshLoadReferenceNodesTriP2D.o 

# library objects
LOBJS = \
../../src/meshParallelReaderTri2D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPrint2D.o \
../../src/meshReader2D.o \
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o

acousticsMain2D:$(AOBJS) $(LOBJS)
	$(LD)  $(LDFLAGS)  -o acousticsMain2D $(AOBJS) $(LOBJS) $(paths) $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r $(ABJS) $(LOBJS)


