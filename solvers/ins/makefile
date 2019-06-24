
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

# define LIBP_OPT_FLAGS
include ../../make.inp

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib
OGSDIR  = ../../libs/gatherScatter
ALMONDDIR = ../../libs/parAlmond
ELLIPTICDIR = ../elliptic

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(OGSDIR) -I$(ELLIPTICDIR) -I$(ALMONDDIR) -D DHOLMES='"${CURDIR}/../.."' -D DINS='"${CURDIR}"' $(LIBP_OPT_FLAGS) 

# -I../../3rdParty/parmetis-4.0.3/metis/GKlib -I../../3rdParty/parmetis-4.0.3/include -I../../3rdParty/parmetis-4.0.3/metis/include -I../../3rdParty/parmetis-4.0.3/libparmetis

# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) $(LIBP_OPT_FLAGS)

#PARMETISDIR = ../../3rdParty/parmetis-4.0.3/build/Linux-x86_64/libparmetis/
#METISDIR = ../../3rdParty/parmetis-4.0.3/build/Linux-x86_64/libmetis

# libraries to be linked in
LIBS	=  -L$(ELLIPTICDIR) -lelliptic -L$(ALMONDDIR) -lparAlmond  \
		   -L$(OGSDIR) -logs -L$(GSDIR)/lib  -lgs \
		   -L$(OCCA_DIR)/lib -locca $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran \

#-L$(PARMETISDIR) -lparmetis -L$(METISDIR) -lmetis


INCLUDES = ins.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(OGSDIR)/ogs.hpp \
$(ALMONDDIR)/parAlmond.hpp \
$(ELLIPTICDIR)/elliptic.h \
$(ELLIPTICDIR)/ellipticPrecon.h \
$(ELLIPTICDIR)/ellipticMultiGrid.h

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
./src/insPlotWallsVTUHex3D.o \
./src/insPlotVTUHex3D.o \
./src/insSetup.o \
./src/insPlotVTU.o \
./src/insError.o \
./src/insForces.o \
./src/insComputeDt.o \
./src/insReport.o \
./src/insRunARK.o \
./src/insRunEXTBDF.o \
./src/insAdvection.o \
./src/insDiffusion.o \
./src/insGradient.o \
./src/insDivergence.o \
./src/insSubCycle.o \
./src/insVelocityRhs.o \
./src/insVelocitySolve.o \
./src/insVelocityUpdate.o \
./src/insPressureRhs.o \
./src/insPressureSolve.o \
./src/insPressureUpdate.o \
./src/insRestart.o \
./src/insWeldTriVerts.o \
./src/insIsoPlotVTU.o \
./src/insBrownMinionQuad3D.o\
./src/insFilterSetup.o \
./src/insAddVelocityRhs.o

# library objects
LOBJS = \
../../src/meshConnectPeriodicFaceNodes3D.o \
../../src/meshConnectPeriodicFaceNodes2D.o \
../../src/meshSetupBoxHex3D.o \
../../src/meshSetupBoxQuad2D.o \
../../src/meshRecursiveSpectralBisectionPartition.o \
../../src/meshApplyElementMatrix.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshGeometricFactorsTri3D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshGeometricFactorsQuad3D.o \
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
../../src/meshOccaSetupQuad3D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelConsecutiveGlobalNumbering.o\
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelReaderTri2D.o \
../../src/meshParallelReaderTri3D.o \
../../src/meshParallelReaderQuad2D.o \
../../src/meshParallelReaderQuad3D.o \
../../src/meshParallelReaderTet3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshPhysicalNodesTri3D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshPhysicalNodesQuad3D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshPlotVTU2D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint2D.o \
../../src/meshPrint3D.o \
../../src/meshSetup.o \
../../src/meshSetupTri2D.o \
../../src/meshSetupTri3D.o \
../../src/meshSetupQuad2D.o \
../../src/meshSetupQuad3D.o \
../../src/meshSetupTet3D.o \
../../src/meshSetupHex3D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshSurfaceGeometricFactorsTri3D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshSurfaceGeometricFactorsQuad3D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshVTU2D.o \
../../src/meshVTU3D.o \
../../src/matrixInverse.o \
../../src/matrixConditionNumber.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
../../src/setupAide.o \
../../src/readArray.o\
../../src/occaDeviceConfig.o\
../../src/occaHostMallocPinned.o \
../../src/timer.o


insMain:$(AOBJS) $(LOBJS) ./src/insMain.o libblas libogs libparAlmond libelliptic
	$(LD)  $(LDFLAGS)  -o insMain ./src/insMain.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libins.a $(AOBJS)

libogs:
	cd ../../libs/gatherScatter; make -j lib; cd ../../solvers/ins

libblas:
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../solvers/ins

libparAlmond:
	cd ../../libs/parAlmond; make -j lib; cd ../../solvers/ins

libelliptic:
	cd ../elliptic; make -j lib; cd ../ins

all: lib insMain

# what to do if user types "make clean"
clean:
	cd ../elliptic; make clean; cd ../ins
	cd ../../src; rm *.o; cd ../solvers/ins
	rm ./src/*.o insMain libins.a

realclean:
	cd ../elliptic; make realclean; cd ../ins
	cd ../../src; rm *.o; cd ../solvers/ins
	rm ./src/*.o insMain libins.a
