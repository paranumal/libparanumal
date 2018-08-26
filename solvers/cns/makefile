ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = -O3

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib
OGSDIR  = ../../libs/gatherScatter

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(OGSDIR) -g  -D DHOLMES='"${CURDIR}/../.."' -D DCNS='"${CURDIR}"'

# link flags to be used 
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g

# libraries to be linked in
LIBS	=   -L$(OGSDIR) -logs -L$(GSDIR)/lib  -lgs \
			-L$(OCCA_DIR)/lib $(links)

INCLUDES = cns.h

DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h \
$(OGSDIR)/ogs.hpp 

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
OBJS    = \
./src/cnsEstimate.o \
./src/cnsBodyForce.o \
./src/cnsStep.o \
./src/cnsMain.o \
./src/cnsError.o \
./src/cnsForces.o \
./src/cnsRun.o \
./src/cnsSetup.o \
./src/cnsGaussianPulse.o \
./src/cnsPlotVTU.o \
./src/cnsReport.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshLoadReferenceNodesTet3D.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshOccaSetup2D.o \
../../src/meshOccaSetup3D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderTri2D.o \
../../src/meshParallelReaderQuad2D.o \
../../src/meshParallelReaderTet3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshPlotVTU2D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint2D.o \
../../src/meshPrint3D.o \
../../src/meshSetupTri2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshSetupTet3D.o \
../../src/meshSetupHex3D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshVTU2D.o \
../../src/meshVTU3D.o \
../../src/mysort.o \
../../src/parallelSort.o \
../../src/setupAide.o \
../../src/trace.o \
../../src/readArray.o \
../../src/occaDeviceConfig.o \
../../src/occaHostMallocPinned.o \
../../src/timer.o


cnsMain:$(OBJS) libogs
	$(LD)  $(LDFLAGS)  -o cnsMain $(OBJS) $(paths) $(LIBS) 

libogs:
	cd ../../libs/gatherScatter; make -j lib; cd ../../solvers/cns

# what to do if user types "make clean"
clean :
	rm -r $(OBJS) cnsMain


