ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR = ../../include
OGSDIR  = ./include
GSDIR  = ../../3rdParty/gslib.github/src

# set options for this machine
# specify which compilers to use for c, fortran and linking
cc	= mpicc
CC	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR) -g  -D DHOLMES='"${CURDIR}/../.."' -D DOGS='"${CURDIR}"'


# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g -L../../3rdParty/gslib.github  -lgs 

# libraries to be linked in
LIBS	=   -L$(OCCA_DIR)/lib  $(links) 

#-llapack -lblas

INCLUDES = ogs.h ogsKernels.h ogsInterface.h
DEPS = $(INCLUDES) \
$(HDRDIR)/types.h 

# types of files we are going to construct rules for
.SUFFIXES: .c

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

# list of objects to be compiled
OBJS    = \
./src/ogsGatherScatter.o \
./src/ogsGather.o \
./src/ogsScatter.o \
./src/ogsSetup.o \
./src/ogsKernels.o 

COBJS = \
./src/ogsHostGatherScatter.o \
./src/ogsHostGather.o \
./src/ogsHostScatter.o \
./src/ogsHostSetup.o


lib:$(OBJS) gslibInterface
	ar -cr libogs.a $(OBJS) $(COBJS)

gslibInterface:
	cd ../../3rdParty/gslib.github; make -j libgs.a; cd ../../libs/gatherScatter
	$(cc) $(CFLAGS) -c -o ./src/ogsHostGatherScatter.o ./src/ogsHostGatherScatter.c $(paths)
	$(cc) $(CFLAGS) -c -o ./src/ogsHostGather.o        ./src/ogsHostGather.c        $(paths)
	$(cc) $(CFLAGS) -c -o ./src/ogsHostScatter.o       ./src/ogsHostScatter.c       $(paths)
	$(cc) $(CFLAGS) -c -o ./src/ogsHostSetup.o         ./src/ogsHostSetup.c         $(paths)

all: lib

# what to do if user types "make clean"
clean:
	rm src/*.o libogs.a

realclean:
	cd ../../3rdParty/gslib.github; make clean; cd ../../libs/gatherScatter
	rm src/*.o libogs.a

