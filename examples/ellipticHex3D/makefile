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
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I../../3rdParty/gslib.github/src

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github -lgs

# libraries to be linked in
LIBS	=  $(links)

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
ellipticSetupHex3D.o \
ellipticPreconditionerSetupHex3D.o

# library objects
LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary3D.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange3D.o \
../../src/meshHaloSetup3D.o \
../../src/meshParallelConnectOpt3D.o \
../../src/meshParallelPrint3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics3D.o \
../../src/meshParallelConnectNodesHex3D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint3D.o \
../../src/meshVTU3D.o \
../../src/meshSetupHex3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshParallelGatherScatter3D.o \
../../src/meshParallelGatherScatterSetup3D.o \
../../src/meshOccaSetup3D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o


COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o

ellipticMainHex3D:$(AOBJS) $(LOBJS) ellipticMainHex3D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o ellipticMainHex3D ellipticMainHex3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 

ellipticMainSplitPml3D:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o ellipticMainSplitPml3D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

test3D: $(LOBJS) $(AOBJS) testMainHex3D.o
	$(LD)  $(LDFLAGS)  -o testMain3D testMainHex3D.o $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ellipticMainHex3D.o


