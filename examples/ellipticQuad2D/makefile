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
ellipticSetupQuad2D.o \
ellipticPreconditionerSetupQuad2D.o

# library objects
LOBJS = \
../../src/meshConnect.o \
../../src/meshConnectBoundary2D.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloSetup.o \
../../src/meshParallelConnectOpt2D.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderQuad2D.o \
../../src/meshPartitionStatistics2D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshParallelGatherScatter2D.o \
../../src/meshParallelGatherScatterSetup2D.o \
../../src/meshOccaSetup2D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
./diagnostic.o


COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o

ellipticMainQuad2D:$(AOBJS) $(LOBJS) ellipticMainQuad2D.o gslibInterface
	$(LD)  $(LDFLAGS)  -o ellipticMainQuad2D ellipticMainQuad2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 

ellipticMainSplitPml2D:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o ellipticMainSplitPml2D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

test2D: $(LOBJS) $(AOBJS) testMainQuad2D.o
	$(LD)  $(LDFLAGS)  -o testMain2D testMainQuad2D.o $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ellipticMainQuad2D.o


