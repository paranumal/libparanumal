
ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = 

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR  = ../../include
GSDIR  = ../../3rdParty/gslib.github/src
ALMONDDIR = ../parALMOND
ELLIPTICDIR = ../ellipticTet3D

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -g -I$(GSDIR) -I$(ELLIPTICDIR) -D DHOLMES='"${CURDIR}/../.."'  

# link flags to be used 
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND -L$(ELLIPTICDIR) -lellipticTet3D

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran


INCLUDES = ins3D.h
DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h \
$(HDRDIR)/ogs_t.h \
$(HDRDIR)/hgs_t.h \
$(ALMONDDIR)/parALMOND.h \
$(ELLIPTICDIR)/ellipticTet3D.h \
$(ELLIPTICDIR)/ellipticPreconTet3D.h  

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
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

insMainTet3D:$(AOBJS) $(LOBJS) insMainTet3D.o gslibInterface
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../examples/insTet3D
	cd ../parALMOND; make -j lib; cd ../insTet3D
	cd ../ellipticTet3D; make -j lib; cd ../insTet3D
	$(LD)  $(LDFLAGS)  -o insMainTet3D insMainTet3D.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libinsTet3D.a $(AOBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../examples/insTet3D
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatter.o ../../src/gsParallelGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ../../src/gsParallelGatherScatterSetup.o ../../src/gsParallelGatherScatterSetup.c $(paths) 
	$(cc) $(CFLAGS) -c -o ../../src/xxtCoarseSolve.o ../../src/xxtCoarseSolve.c $(paths)

all: lib insMainTet3D

# what to do if user types "make clean"
clean:
	cd ../parALMOND; make clean; cd ../insTet3D
	cd ../ellipticTet3D; make clean; cd ../insTet3D
	cd ../../src; rm *.o; cd ../examples/insTet3D
	rm *.o insMainTet3D libinsTet3D.a

reallyclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../examples/insTet3D
	cd ../../3rdParty/BlasLapack; make clean; cd ../../examples/insTet3D
	cd ../parALMOND; make clean; cd ../insTet3D
	cd ../ellipticTet3D; make clean; cd ../insTet3D
	cd ../../src; rm *.o; cd ../examples/insTet3D
	rm *.o insMainTet3D libinsTet3D.a
