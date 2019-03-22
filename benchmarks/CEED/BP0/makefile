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

testInnerProduct: testInnerProduct.o
	$(LD)  $(LDFLAGS)  -o testInnerProduct testInnerProduct.o $(paths) $(LIBS)

testScaledAdd: testScaledAdd.o
	$(LD)  $(LDFLAGS)  -o testScaledAdd testScaledAdd.o $(paths) $(LIBS)

# what to do if user types "make clean"
clean :
	rm -r testInnerProduct testScaledAdd testInnerProduct.o testScaledAdd.o


