ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ../../make.inp

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR = ../../include
GSDIR  = ../../3rdParty/gslib
OGSDIR  = ../../libs/gatherScatter
ALMONDDIR = ../../libs/parAlmond

# Building p4est:
# $ tar xzvf p4est-2.0.tar.gz
# $ cd p4est-2.0
# $ ./configure CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 --enable-mpi --without-blas
# $ make -j install
P4ESTDIR = ../../3rdParty/p4est-2.0
include $(P4ESTDIR)/local/etc/Makefile.p4est.mk

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
FC	= mpif77  
#-std=legacy -fdefault-real-8 -fdefault-double-8 
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I./include -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) $(P4EST_CPPFLAGS) -I../../include -D DHOLMES='"${CURDIR}/../.."' -D DADAPTIVE='"${CURDIR}"' $(LIBP_OPT_FLAGS)

# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) $(LIBP_OPT_FLAGS)

# libraries to be linked in
LIBS	= -Wl,-rpath=$(P4ESTDIR)/local/lib,--enable-new-dtags $(P4EST_LDFLAGS) $(P4EST_LIBS) -Wl,-rpath=$(OCCA_DIR)/lib,--enable-new-dtags -L$(OCCA_DIR)/lib -locca $(links) -L../../3rdParty/BlasLapack  $(LIBP_EXT_LIBS) -lgfortran

INCLUDES = adaptive.h \
           app.h \
           asd.h \
           connectivity.h \
           myoccautil.h \
           level.h \
           mesh.h \
           number_types.h \
           operators.h \
           solver_info.h
DEPS = $(INCLUDES)

# types of files we are going to construct rules for
.SUFFIXES: .c .f

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

.f.o: $(DEPS)
	$(FC) $(CFLAGS) -o $*.o -c $*.f $(paths)

# add this object for USE_STEFAN_MXM ./src/adaptiveSerialLocalAxHex3D.o\


# list of objects to be compiled
AOBJS    = \
./src/adaptiveOutputVTUHex3D.o \
./src/asd.o \
./src/app.o \
./src/connectivity.o \
./src/level.o\
./src/mesh.o\
./src/operators.o\
./src/occaDeviceConfig.o

# library objects
LOBJS = \
../../src/setupAide.o \
../../src/timer.o

adaptiveMain:$(AOBJS) $(LOBJS) ./src/adaptiveMain.o libblas
	$(LD)  $(LDFLAGS)  -o adaptiveMain ./src/adaptiveMain.o $(COBJS) $(AOBJS) $(LOBJS) $(paths) $(LIBS)

lib:$(AOBJS)
	ar -cr libadaptive.a $(AOBJS)

libogs:
	cd ../../libs/gatherScatter; make -j lib; cd ../../solvers/adaptive

libblas:
	cd ../../3rdParty/BlasLapack; make -j lib; cd ../../solvers/adaptive

libparAlmond:
	cd ../../libs/parAlmond; make -j lib; cd ../../solvers/adaptive

all: lib adaptiveMain

# what to do if user types "make clean"
clean:
	cd ../../libs/parAlmond; make clean; cd ../../solvers/adaptive
	cd ../../src; rm *.o; cd ../solvers/adaptive
	cd ../../libs/gatherScatter; make clean; cd ../../solvers/adaptive
	rm src/*.o adaptiveMain libadaptive.a

realclean:
	cd ../../3rdParty/BlasLapack; make clean; cd ../../solvers/adaptive
	cd ../../libs/gatherScatter; make realclean; cd ../../solvers/adaptive
	cd ../../libs/parAlmond; make clean; cd ../../solvers/adaptive
	cd ../../src; rm *.o; cd ../solvers/adaptive
	rm src/*.o adaptiveMain libadaptive.a

