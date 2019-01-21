ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS = -ftree-loop-distribution -floop-nest-optimize -floop-block -floop-interchange -fstrict-aliasing -ftree-vectorize -funroll-loops -O3 

include ${OCCA_DIR}/scripts/Makefile
MEMCOPY_ON = 0

# define variables
HDRDIR = ../../include
OGSDIR  = ./include
GSDIR  = ../../3rdParty/gslib

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpicc
CXX	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(GSDIR)/src -g  -D DHOLMES='"${CURDIR}/../.."' -D DOGS='"${CURDIR}"'

CFLAGS+= -DMEMCOPY=$(MEMCOPY_ON) 
# link flags to be used
LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g 

# libraries to be linked in
LIBS	=   -L$(OCCA_DIR)/lib  $(links) -L$(GSDIR)/lib  -lgs 

#-llapack -lblas

INCLUDES = ogs.hpp ogsKernels.hpp ogsInterface.h
DEPS = $(INCLUDES) \
$(HDRDIR)/types.h 

# types of files we are going to construct rules for
.SUFFIXES: .c .cpp

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths)

.cpp.o: $(DEPS)
	$(CXX) $(CFLAGS) -o $*.o -c $*.cpp $(paths)

# list of objects to be compiled
OBJS    = \
./src/ogsGatherScatter.o \
./src/ogsGatherScatterVec.o \
./src/ogsGatherScatterMany.o \
./src/ogsGather.o \
./src/ogsGatherVec.o \
./src/ogsGatherMany.o \
./src/ogsScatter.o \
./src/ogsScatterVec.o \
./src/ogsScatterMany.o \
./src/ogsSetup.o \
./src/ogsKernels.o 

COBJS = \
./src/ogsHostGatherScatter.o \
./src/ogsHostGatherScatterVec.o \
./src/ogsHostGatherScatterMany.o \
./src/ogsHostGather.o \
./src/ogsHostGatherVec.o \
./src/ogsHostGatherMany.o \
./src/ogsHostScatter.o \
./src/ogsHostScatterVec.o \
./src/ogsHostScatterMany.o \
./src/ogsHostSetup.o


lib:$(OBJS) gslibInterface
	ar -cr libogs.a $(OBJS) $(COBJS)

gslibInterface:
	cd ../../3rdParty/gslib; CC=$(CC) make -j install; cd ../../libs/gatherScatter
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGatherScatter.o     ./src/ogsHostGatherScatter.c     $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGatherScatterVec.o  ./src/ogsHostGatherScatterVec.c  $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGatherScatterMany.o ./src/ogsHostGatherScatterMany.c $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGather.o        ./src/ogsHostGather.c        $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGatherVec.o     ./src/ogsHostGatherVec.c     $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostGatherMany.o    ./src/ogsHostGatherMany.c    $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostScatter.o       ./src/ogsHostScatter.c       $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostScatterVec.o    ./src/ogsHostScatterVec.c    $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostScatterMany.o   ./src/ogsHostScatterMany.c   $(paths)
	$(CC) $(CFLAGS) -c -o ./src/ogsHostSetup.o         ./src/ogsHostSetup.c         $(paths)

all: lib

# what to do if user types "make clean"
clean:
	rm src/*.o libogs.a

realclean:
	cd ../../3rdParty/gslib; make clean; cd ../../libs/gatherScatter
	rm src/*.o libogs.a

