ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR)  -D DHOLMES='"${CURDIR}/../.."'

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags)

# libraries to be linked in
LIBS	=  $(links)

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
AOBJS    = \
boltzmannMain3D.o \
boltzmannSetup3D.o \
boltzmannRun3D.o \
boltzmannReport3D.o \
boltzmannError3D.o \
boltzmannComputeVorticity3D.o \
boltzmannRampFunction3D.o \
boltzmannLserkStep3D.o \
boltzmannSark3Step3D.o \
boltzmannSaab3Step3D.o \
boltzmannLsimexStep3D.o 

# library objects
LOBJS = \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshLoadReferenceNodesTet3D.o \
../../src/meshParallelReaderTet3D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshSetupTet3D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshOccaSetup3D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/mysort.o \
../../src/parallelSort.o

boltzmannMain3D:$(AOBJS) $(LOBJS)
	$(LD)  $(LDFLAGS)  -o boltzmannMain3D $(AOBJS) $(LOBJS) $(paths) $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS)


