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
# list of objects to be compiled
AOBJS    = \
insMain2D.o \
insSetup2D.o \
insRun2D.o \
insMakePeriodic2D.o \
insPlotVTU2D.o \
insError2D.o \
insReport2D.o \
insHelmholtzStep2D.o \
insPressureStep2D.o  

LOBJS    = \
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

insMain2D:$(AOBJS) $(LOBJS)
	$(LD)  $(LDFLAGS)  -o insMain2D $(AOBJS) $(LOBJS) $(paths) $(LIBS)
# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS)
 

