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
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -D DHOLMES='"${CURDIR}/../.."' 

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) 

# libraries to be linked in
LIBS	=  $(links)

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
OBJS    = \
boltzmannMainTri3D.o \
boltzmannRunTri3D.o \
boltzmannSetupTri3D.o \
boltzmannPlotVTUTri3DV2.o \
boltzmannOccaSetupTri3D.o \
../../src/meshParallelReaderTri3D.o \
../../src/meshSetupTri3D.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshSurfaceGeometricFactorsTri3D.o \
../../src/meshPhysicalNodesTri3D.o \
../../src/meshParallelConnectNodes.o\
../../src/meshGeometricFactorsTri3D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o \
../../src/hash.o

boltzmannMainTri3D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o boltzmannMainTri3D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)


