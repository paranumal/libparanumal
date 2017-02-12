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
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR)

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
acousticsMain3D.o \
acousticsError3D.o \
acousticsRun3D.o \
acousticsSetup3D.o \
acousticsSurface3D.o \
acousticsUpdate3D.o \
acousticsVolume3D.o \
acousticsCavitySolution3D.o


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

acousticsMain3D:$(AOBJS) $(LOBJS)
	$(LD)  $(LDFLAGS)  -o acousticsMain3D $(AOBJS) $(LOBJS) $(paths) $(LIBS)

acousticsMainSplitPml3D:$(POBJS) $(LOBJS)
	$(LD)  $(LDFLAGS)  -o acousticsMainSplitPml3D $(POBJS) $(LOBJS) $(paths) $(LIBS) 

# what to do if user types "make clean"
clean :
	rm -r $(AOBJS) $(LOBJS) $(POBJS)


