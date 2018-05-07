ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/makefile

# define variables
HDRDIR = ../../include

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc.openmpi
CC	= mpic++.openmpi
LD	= mpic++.openmpi

# compiler flags to be used (set to compile with debugging on)
CFLAGS = $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -g  -D DHOLMES='"${CURDIR}/../.."'


# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs \
			-L$(ALMONDDIR) -lparALMOND

# libraries to be linked in
LIBS	=  $(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

#-llapack -lblas

INCLUDES = autoTester.h
DEPS = $(INCLUDES) 

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# library objects
LOBJS = \
../../src/setupAide.o \
../../src/timer.o\
../../src/trace.o 

autoTester:../../src/matrix.o  ../../src/setupAide.o ../../src/trace.o autoTester.o
	$(LD) $(LDFLAGS) -o autoTester autoTester.o ../../src/matrix.o  ../../src/trace.o ../../src/setupAide.o $(paths) $(LIBS)

all: autoTester

# what to do if user types "make clean"
clean:
	rm *.o autoTester

