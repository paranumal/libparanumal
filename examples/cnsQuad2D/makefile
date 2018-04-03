ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = -g

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR)  -D DHOLMES='"${CURDIR}/../.."'

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
cnsMainQuad2D.o \
cnsError2D.o \
cnsUpdate2D.o \
cnsRunQuad2D.o \
cnsSetupQuad2D.o \
cnsSurfaceQuad2D.o \
cnsVolumeQuad2D.o \
cnsCavitySolution2D.o \
cnsGaussianPulse2D.o \
cnsPlotVTUQuad2D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshOccaSetup2D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderQuad2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshVTU2D.o \
../../src/mysort.o \
../../src/parallelSort.o

cnsMainQuad2D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o cnsMainQuad2D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS)


