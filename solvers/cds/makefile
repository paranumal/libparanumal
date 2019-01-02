
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib
OGSDIR  = ../../libs/gatherScatter
ALMONDDIR = ../../libs/parAlmond
ELLIPTICDIR = ../elliptic
INSDIR = ../ins

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(OGSDIR) -I$(ELLIPTICDIR) -I$(INSDIR) -I$(ALMONDDIR) -g  -D DHOLMES='"${CURDIR}/../.."' -D DCDS='"${CURDIR}"'  -D DINS='"${INSDIR}"'

# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g

# libraries to be linked in
LIBS	=  -L$(ELLIPTICDIR) -lelliptic -L$(INS) -lins -L$(ALMONDDIR) -lparAlmond  \
		   -L$(OGSDIR) -logs -L$(GSDIR)/lib  -lgs \
		   -L$(OCCA_DIR)/lib $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran \


INCLUDES = cds.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(OGSDIR)/ogs.hpp \
$(ALMONDDIR)/parAlmond.hpp \
$(ELLIPTICDIR)/elliptic.h \
$(ELLIPTICDIR)/ellipticPrecon.h \
$(ELLIPTICDIR)/ellipticMultiGrid.h \
$(INSDIR)/ins.h

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
./src/cdsSetup.o \
./src/cdsRunEXTBDF.o \
./src/cdsAdvection.o \
./src/cdsHelmholtzRhs.o \
./src/cdsHelmholtzSolve.o \
./src/cdsError.o \
./src/cdsReport.o \
./src/cdsPlotVTU.o \
./src/cdsSubCycle.o \
./src/cdsSolveSetup.o \
./src/cdsSolveStep.o

# library objects
LOBJS = \
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


cdsMain:$(AOBJS) $(LOBJS) ./src/cdsMain.o libblas libogs libparAlmond libelliptic
	$(LD)  $(LDFLAGS)  -o cdsMain ./src/cdsMain.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libcds.a $(AOBJS)

libogs:
	cd ../../libs/gatherScatter; make -j lib; cd ../../solvers/cds

libblas:
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../solvers/cds

libparAlmond:
	cd ../../libs/parAlmond; make -j lib; cd ../../solvers/cds

libelliptic:
	cd ../elliptic; make -j lib; cd ../cds

all: lib cdsMain

# what to do if user types "make clean"
clean:
	cd ../elliptic; make clean; cd ../cds
	cd ../../src; rm *.o; cd ../solvers/cds
	rm ./src/*.o cdsMain libcds.a

realclean:
	cd ../elliptic; make realclean; cd ../ins
	cd ../../src; rm *.o; cd ../solvers/ins
	rm ./src/*.o cdsMain libcds.a
