
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
		 	-I../../utilities/parALMOND/include -I../ellipticTet3D

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
insSetup3D.o \
insRun3D.o \
insError3D.o \
insReport3D.o \
insPlotVTU3D.o \
insAdvectionStep3D.o \
insHelmholtzStep3D.o \
insPoissonStep3D.o \
insUpdateStep3D.o \
insErrorNorms3D.o \
../ellipticTet3D/ellipticSetupTet3D.o \
../ellipticTet3D/ellipticPreconditionerSetupTet3D.o\
../ellipticTet3D/ellipticCoarsePreconditionerSetupTet3D.o\
../ellipticTet3D/ellipticSolveTet3D.o\
../ellipticTet3D/ellipticSolveSetupTet3D.o\
../ellipticTet3D/ellipticHaloExchange3D.o\
../ellipticTet3D/ellipticParallelGatherScatterTet3D.o \
../ellipticTet3D/ellipticBuildContinuousTet3D.o \
../ellipticTet3D/ellipticBuildIpdgTet3D.o

LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint3D.o \
../../src/meshParallelReaderTet3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/meshSetupTet3D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshLoadReferenceNodesTet3D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshParallelGather.o \
../../src/meshParallelScatter.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherSetup.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup3D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o

COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o

ins3D:$(AOBJS) $(LOBJS) insMainTet3D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o ins3D insMainTet3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ins3D insMainTet3D.o


#insMakePeriodic2D.o \
#insErrorNorms2D.o \
#insAdvectionSubCycleStep2D.o \

#insUpdateStep2D.o \
#insBuildVectorIpdgTri2D.o \