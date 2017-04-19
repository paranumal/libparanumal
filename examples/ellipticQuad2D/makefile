ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include

AMGDIR = ../../../AMG2013

AMGINCL  = -I$(AMGDIR) \
		   -I$(AMGDIR)/utilities \
		   -I$(AMGDIR)/parcsr_ls \
		   -I$(AMGDIR)/parcsr_mv \
		   -I$(AMGDIR)/krylov \
		   -I$(AMGDIR)/seq_mv \
		   -I$(AMGDIR)/sstruct_mv \
		   -I$(AMGDIR)/struct_mv \
		   -I$(AMGDIR)/IJ_mv \

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) $(AMGINCL) -O3 -g -I../../3rdParty/gslib.github/src -D DHOLMES='"${CURDIR}/../.."' -I../../utilities/parALMOND/include

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -O3 -L../../3rdParty/gslib.github -lgs \
			-L../../3rdParty/amg2013 -lamg2013 \
			-L../../utilities/parALMOND -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -llapack -lblas

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
AOBJS    = \
ellipticSetupQuad2D.o \
ellipticPreconditionerSetupQuad2D.o\
ellipticCoarsePreconditionerSetupQuad2D.o\
ellipticCoarsePreconditionerQuad2D.o\
../../src/almondCoarseSolve.o

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
../../src/meshParallelReaderQuad2D.o \
../../src/meshPartitionStatistics.o \
../../src/meshParallelConnectNodes.o \
../../src/meshPlotVTU2D.o \
../../src/meshPrint2D.o \
../../src/meshVTU2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshParallelGatherScatter.o \
../../src/meshParallelGatherScatterSetup.o \
../../src/meshParallelConsecutiveGlobalNumbering.o \
../../src/meshOccaSetup2D.o \
../../src/mysort.o \
../../src/parallelSort.o\
../../src/hash.o\
./diagnostic.o\


COBJS = \
../../src/gsParallelGatherScatter.o\
../../src/gsParallelGatherScatterSetup.o\
../../src/xxtCoarseSolve.o \
../../src/amg2013CoarseSolve.o 

ellipticMainQuad2D:$(AOBJS) $(LOBJS) ellipticMainQuad2D.o gslibInterface amg2013Interface
	$(LD)  $(LDFLAGS)  -o ellipticMainQuad2D ellipticMainQuad2D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

gslibInterface:
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

amg2013Interface:
	$(cc) $(CFLAGS) -c -o ../../src/amg2013CoarseSolve.o ../../src/amg2013CoarseSolve.c $(paths)

ellipticMainSplitPml2D:$(POBJS) $(LOBJS) 
	$(LD)  $(LDFLAGS)  -o ellipticMainSplitPml2D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

test2D: $(LOBJS) $(AOBJS) testMainQuad2D.o
	$(LD)  $(LDFLAGS)  -o testMain2D testMainQuad2D.o $(AOBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS) $(COBJS) ellipticMainQuad2D.o


