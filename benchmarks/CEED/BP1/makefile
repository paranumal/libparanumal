ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -O3 -I../../../3rdParty/gslib.github/src -D LIBP_DIR='"${CURDIR}/../../.."'  -I../../../utilities/parALMOND/include

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -O3 -L../../../3rdParty/gslib.github -lgs\
			-L../../../utilities/parALMOND -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../../utilities/BlasLapack -lBlasLapack

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
massSetupHex3D.o \
massSolveSetupHex3D.o \
massHaloExchange3D.o \
massParallelGatherScatter.o \
massParallelGatherScatterSetup.o \
massSolveHex3D.o
#massPipelinedSolveHex3D.o 

# library objects
LOBJS = \
../../../src/meshConnect.o \
../../../src/meshConnectBoundary.o \
../../../src/meshConnectFaceNodes3D.o \
../../../src/meshGeometricPartition3D.o \
../../../src/meshHaloExchange.o \
../../../src/meshHaloSetup.o \
../../../src/meshParallelConnectOpt.o \
../../../src/meshParallelPrint3D.o \
../../../src/meshParallelReaderHex3D.o \
../../../src/meshPartitionStatistics.o \
../../../src/meshParallelConnectNodes.o \
../../../src/meshPlotVTU3D.o \
../../../src/meshPrint3D.o \
../../../src/meshVTU3D.o \
../../../src/meshSetupHex3D.o \
../../../src/meshPhysicalNodesHex3D.o \
../../../src/meshGeometricFactorsHex3D.o \
../../../src/meshLoadReferenceNodesHex3D.o \
../../../src/meshSurfaceGeometricFactorsHex3D.o \
../../../src/meshParallelGather.o \
../../../src/meshParallelScatter.o \
../../../src/meshParallelGatherScatter.o \
../../../src/meshParallelGatherSetup.o \
../../../src/meshParallelGatherScatterSetup.o \
../../../src/meshParallelConsecutiveGlobalNumbering.o \
../../../src/meshOccaSetup3D.o \
../../../src/mysort.o \
../../../src/parallelSort.o\
../../../src/hash.o \
../../../src/timer.o 
#../../../src/mpistubs.o

COBJS = \
../../../src/gsParallelGatherScatter.o\
../../../src/gsParallelGatherScatterSetup.o\
../../../src/xxtCoarseSolve.o 

massMainHex3D:$(AOBJS) $(LOBJS) massMainHex3D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o massMainHex3D massMainHex3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../../src/gsParallelGatherScatter.o ../../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../../src/gsParallelGatherScatterSetup.o ../../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../../src/xxtCoarseSolve.o ../../../src/xxtCoarseSolve.c $(paths)

BP3:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o BP3 $(POBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) massMainHex3D.o


