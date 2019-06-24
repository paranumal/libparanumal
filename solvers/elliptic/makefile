ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ../../make.inp

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR = ../../include
GSDIR  = ../../3rdParty/gslib
OGSDIR  = ../../libs/gatherScatter
ALMONDDIR = ../../libs/parAlmond

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
FC	= mpif77  
#-std=legacy -fdefault-real-8 -fdefault-double-8 
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(OGSDIR) -I$(ALMONDDIR) -D DHOLMES='"${CURDIR}/../.."' -D DELLIPTIC='"${CURDIR}"'  $(LIBP_OPT_FLAGS) -I$(OGSDIR)/include

# -I../../3rdParty/parmetis-4.0.3/metis/GKlib -I../../3rdParty/parmetis-4.0.3/include -I../../3rdParty/parmetis-4.0.3/metis/include -I../../3rdParty/parmetis-4.0.3/libparmetis

# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) $(LIBP_OPT_FLAGS)

#PARMETISDIR = ../../3rdParty/parmetis-4.0.3/build/Linux-x86_64/libparmetis/
#METISDIR = ../../3rdParty/parmetis-4.0.3/build/Linux-x86_64/libmetis

# libraries to be linked in
LIBS	=    -L$(ALMONDDIR) -lparAlmond  -L$(OGSDIR) -logs -L$(GSDIR)/lib -lgs -L$(OCCA_DIR)/lib -locca $(links) -L../../3rdParty/BlasLapack  $(LIBP_EXT_LIBS) -lgfortran

#-L$(PARMETISDIR) -lparmetis -L$(METISDIR) -lmetis

INCLUDES = elliptic.h ellipticPrecon.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(OGSDIR)/ogs.hpp \
$(ALMONDDIR)/parAlmond.hpp \

# types of files we are going to construct rules for
.SUFFIXES: .c .f

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

.f.o: $(DEPS)
	$(FC) $(CFLAGS) -o $*.o -c $*.f $(paths)

# add this object for USE_STEFAN_MXM ./src/ellipticSerialLocalAxHex3D.o\


# list of objects to be compiled
AOBJS    = \
./src/ellipticThinOas.o \
./src/ellipticThinOasSetup.o \
./src/ellipticZeroMean.o \
./src/NBFPCG.o \
./src/ellipticUpdateNBFPCG.o \
./src/NBPCG.o \
./src/ellipticUpdateNBPCG.o \
./src/ellipticKernelInfo.o \
./src/PCG.o \
./src/ellipticSerialAxHex3D.o \
./src/ellipticSerialPartialAxHex3D.o \
./src/ellipticPlotVTUHex3D.o \
./src/ellipticBuildContinuous.o \
./src/ellipticBuildIpdg.o \
./src/ellipticBuildJacobi.o \
./src/ellipticBuildLocalPatches.o \
./src/ellipticBuildMultigridLevel.o \
./src/ellipticHaloExchange.o\
./src/ellipticOperator.o \
./src/ellipticPreconditioner.o\
./src/ellipticPreconditionerSetup.o\
./src/ellipticSetup.o \
./src/ellipticSolve.o\
./src/ellipticSolveSetup.o\
./src/ellipticScaledAdd.o \
./src/ellipticUpdatePCG.o \
./src/ellipticWeightedInnerProduct.o \
./src/ellipticWeightedNorm2.o \
./src/ellipticVectors.o \
./src/ellipticSEMFEMSetup.o\
./src/ellipticMultiGridSetup.o \
./src/ellipticMultiGridLevel.o \
./src/ellipticMultiGridLevelSetup.o \

# library objects
LOBJS = \
../../src/meshConnectPeriodicFaceNodes2D.o \
../../src/meshConnectPeriodicFaceNodes3D.o \
../../src/meshSetupBoxHex3D.o \
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
../../src/meshRecursiveSpectralBisectionPartition.o\
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
../../src/meshOccaSetupTri3D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelReaderTri2D.o \
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
../../src/parallelSort.o \
../../src/setupAide.o \
../../src/readArray.o\
../../src/occaDeviceConfig.o\
../../src/occaHostMallocPinned.o \
../../src/timer.o

ellipticMain:$(AOBJS) $(LOBJS) ./src/ellipticMain.o libblas libogs libparAlmond
	$(LD)  $(LDFLAGS)  -o ellipticMain ./src/ellipticMain.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libelliptic.a $(AOBJS)

libogs:
	cd ../../libs/gatherScatter; make -j lib; cd ../../solvers/elliptic

libblas:
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../solvers/elliptic

libparAlmond:
	cd ../../libs/parAlmond; make -j lib; cd ../../solvers/elliptic

all: lib ellipticMain

# what to do if user types "make clean"
clean:
	cd ../../libs/parAlmond; make clean; cd ../../solvers/elliptic
	cd ../../src; rm *.o; cd ../solvers/elliptic
	cd ../../libs/gatherScatter; make clean; cd ../../solvers/elliptic
	rm src/*.o ellipticMain libelliptic.a

realclean:
	cd ../../3rdParty/BlasLapack; make clean; cd ../../solvers/elliptic
	cd ../../libs/gatherScatter; make realclean; cd ../../solvers/elliptic
	cd ../../libs/parAlmond; make clean; cd ../../solvers/elliptic
	cd ../../src; rm *.o; cd ../solvers/elliptic
	rm src/*.o ellipticMain libelliptic.a

