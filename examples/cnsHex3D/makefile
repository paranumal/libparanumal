ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = -O3

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC  = mpic++.openmpi
LD  = mpic++.openmpi

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g  -D DHOLMES='"${CURDIR}/../.."'

# link flags to be used 
LDFLAGS = $(compilerFlags) $(flags) -g

# libraries to be linked in
LIBS  =  $(links)

INCLUDES = cnsHex3D.h

DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
OBJS    = \
cnsEstimateHex3D.o \
cnsStepHex3D.o \
cnsMainHex3D.o \
cnsError3D.o \
cnsRunHex3D.o \
cnsSetupHex3D.o \
cnsCavitySolution3D.o \
cnsGaussianPulse3D.o \
cnsPlotVTUHex3D.o \
cnsReportHex3D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshOccaSetup3D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshSetupHex3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshVTU3D.o \
../../src/mysort.o \
../../src/parallelSort.o \
../../src/matrix.o \
../../src/setupAide.o \
../../src/trace.o

cnsMainHex3D:$(OBJS) 
	$(LD)  $(LDFLAGS)  -o cnsMainHex3D $(OBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(OBJS) cnsMainHex3D


