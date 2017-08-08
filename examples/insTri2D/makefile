
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I../../3rdParty/gslib.github/src -D DHOLMES='"${CURDIR}/../.."'  \
		 	-I../../utilities/parALMOND/include -I../ellipticTri2D

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L../../utilities/parALMOND -lparALMOND 

# libraries to be linked in
LIBS	=  $(links) -L../../utilities/BlasLapack -lBlasLapack -lgfortran


# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
insSetup2D.o \
insRun2D.o \
insMakePeriodic2D.o \
insPlotVTU2D.o \
insError2D.o \
insErrorNorms2D.o \
insReport2D.o \
insAdvectionStep2D.o \
insAdvectionSubCycleStep2D.o \
insHelmholtzStep2D.o \
insPoissonStep2D.o \
insUpdateStep2D.o \
insBuildVectorIpdgTri2D.o \
../ellipticTri2D/ellipticSetupTri2D.o \
../ellipticTri2D/ellipticPreconditionerSetupTri2D.o\
../ellipticTri2D/ellipticCoarsePreconditionerSetupTri2D.o\
../ellipticTri2D/ellipticPreconditionerTri2D.o\
../ellipticTri2D/ellipticSolveTri2D.o\
../ellipticTri2D/ellipticSolveSetupTri2D.o\
../ellipticTri2D/ellipticHaloExchange2D.o\
../ellipticTri2D/ellipticParallelGatherScatterTri2D.o \
../ellipticTri2D/ellipticParallelGatherScatterSetup.o \
../ellipticTri2D/ellipticBuildIpdgTri2D.o \
../ellipticTri2D/ellipticBuildContinuousTri2D.o \
../ellipticTri2D/ellipticBuildPatchesIpdgTri2D.o 

# library objects
LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderTri2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/meshSetupTri2D.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshParallelGather.o \
../../src/meshParallelScatter.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherSetup.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup2D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
../../src/timer.o\

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

insMainTri2D:$(AOBJS) $(LOBJS) insMainTri2D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o insMainTri2D insMainTri2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) insMainTri2D.o
