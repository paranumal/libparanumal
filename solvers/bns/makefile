ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib.github/src

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -g  -D DHOLMES='"${CURDIR}/../.."' -D DBNS='"${CURDIR}"'

# link flags to be used 
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) -g -L../../3rdParty/gslib.github  -lgs -fopenmp

# libraries to be linked in
LIBS	=   -L$(OCCA_DIR)/lib $(links)

INCLUDES = bns.h 
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(HDRDIR)/ogs_t.h 

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
AOBJS    = \
./src/bnsMain.o \
./src/bnsSetup.o \
./src/bnsRun.o \
./src/bnsBodyForce.o \
./src/bnsPmlSetup.o \
./src/bnsMRABPmlSetup.o \
./src/bnsTimeStepperCoefficients.o \
./src/bnsSAADRKCoefficients.o \
./src/bnsPlotVTU.o \
./src/bnsReport.o \
./src/bnsError.o \
./src/bnsForces.o \
./src/bnsLSERKStep.o \
./src/bnsSARKStep.o \
./src/bnsMRSAABStep.o \
./src/bnsIsoPlotVTU.o \
./src/bnsIsoWeldPlotVTU.o \
./src/bnsRunEmbedded.o \
./src/bnsWeldTriVerts.o \
./src/bnsIsoPlotGmsh.o \
./src/bnsRestart.o    

# library objects
LOBJS = \
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
../../src/meshMRABSetup2D.o \
../../src/meshMRABSetup3D.o \
../../src/meshBuildMRABClusters2D.o \
../../src/meshBuildMRABClusters3D.o \
../../src/meshClusteredGeometricPartition2D.o \
../../src/meshClusteredGeometricPartition3D.o \
../../src/meshMRABWeightedPartition2D.o \
../../src/meshMRABWeightedPartition3D.o \
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
../../src/meshParallelConsecutiveGlobalNumbering.o\
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/occaDeviceConfig.o\
../../src/occaHostMallocPinned.o \
../../src/timer.o

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\


bnsMain:$(AOBJS) $(LOBJS) ./src/bnsMain.o gslibInterface
	$(LD)  $(LDFLAGS)  -o bnsMain $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS) 



gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../solvers/bns
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 


# what to do if user types "make clean"
clean :
	cd ../../3rdParty/gslib.github; make clean; cd ../../solvers/bns
	rm -r $(AOBJS) $(LOBJS) bnsMain

