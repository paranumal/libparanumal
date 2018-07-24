
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib.github/src
ALMONDDIR = ../parALMOND
ELLIPTICDIR = ../elliptic

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -I$(ELLIPTICDIR) -g  -D DHOLMES='"${CURDIR}/../.."' -D DMPPF='"${CURDIR}"'

# link flags to be used 
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ELLIPTICDIR) -lelliptic -L$(ALMONDDIR) -lparALMOND 

# libraries to be linked in
LIBS	=  -L$(OCCA_DIR)/lib $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

INCLUDES = mppf.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh2D.h \
$(HDRDIR)/mesh3D.h \
$(HDRDIR)/ogs_t.h \
$(ALMONDDIR)/parALMOND.h \
$(ELLIPTICDIR)/elliptic.h \
$(ELLIPTICDIR)/ellipticPrecon.h  

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
./src/mppfSetup.o \
#./src/insPlotVTU.o \
#./src/insError.o \
#./src/insForces.o \
#./src/insComputeDt.o \
#./src/insReport.o \
#./src/insRunARK.o \
#./src/insRunEXTBDF.o \
#./src/insAdvection.o \
#./src/insDiffusion.o \
#./src/insGradient.o \
#./src/insDivergence.o \
#./src/insSubCycle.o \
#./src/insVelocityRhs.o \
#./src/insVelocitySolve.o \
#./src/insVelocityUpdate.o \
#./src/insPressureRhs.o \
#./src/insPressureSolve.o \
#./src/insPressureUpdate.o \
#./src/insRestart.o \
#./src/insWeldTriVerts.o \
#./src/insIsoPlotVTU.o    

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
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelConsecutiveGlobalNumbering.o\
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
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

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

mppfMain:$(AOBJS) $(LOBJS) ./src/mppfMain.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../solvers/mppf
	cd ../parALMOND; make -j lib; cd ../mppf
	cd ../elliptic; make -j lib; cd ../mppf
	$(LD)  $(LDFLAGS)  -o mppfMain ./src/mppfMain.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libmppf.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../solvers/mppf
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib mppfMain

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../mppf
	cd ../elliptic; make clean; cd ../mppf
	cd ../../src; rm *.o; cd ../solvers/mppf
	rm ./src/*.o mppfMain libmppf.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../solvers/mppf
	cd ../../3rdParty/BlasLapack; make clean; cd ../../solvers/mppf
	cd ../parALMOND; make clean; cd ../mppf
	cd ../elliptic; make clean; cd ../mppf
	cd ../../src; rm *.o; cd ../solvers/mppf
	rm ./src/*.o mppfMain libmppf.a
