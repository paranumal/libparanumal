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
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -D DHOLMES='"${CURDIR}/../.."' -O3 -g

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -O3

# libraries to be linked in
LIBS	=  $(links)

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
OBJS    = \
boltzmannMainQuad3D.o \
boltzmannRunMRSAABQuad3D.o \
boltzmannSetupMRQuad3D.o \
boltzmannPlotVTUQuad3D.o \
boltzmannPlotVTUQuad3DV2.o \
boltzmannOccaSetupQuad3D.o \
boltzmannPlotLevels.o \
../../src/meshParallelReaderQuad3D.o \
../../src/meshSetupQuad3D.o \
../../src/meshLoadReferenceNodesQuad3D.o \
../../src/meshSurfaceGeometricFactorsQuad3D.o \
../../src/meshPhysicalNodesQuad3D.o \
../../src/meshParallelConnectNodes.o\
../../src/meshGeometricFactorsQuad3D.o \
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
../../src/hash.o \
../../src/meshMRABWeightedPartitionQuad3D.o \
../../src/meshMRABSetupQuad3D.o \
../../src/meshClusteredGeometricPartition3D.o \
../../src/meshBuildMRABClusters3D.o

boltzmannMainQuad3D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o boltzmannMainQuad3D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)


