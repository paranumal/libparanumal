ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

CXXFLAGS =

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR = ../../
OGSDIR  = ../gatherScatter
GSDIR  = ../../3rdParty/gslib

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpicc
CXX	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -I./include/ $(compilerFlags) $(flags) -I$(HDRDIR)/include/ -I$(OGSDIR) -g -D DPARALMOND='"${CURDIR}"'

# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) -g

# libraries to be linked in
LIBS	=   -L$(OCCA_DIR)/lib  $(links) -L$(OGSDIR) -logs -L$(GSDIR)/lib -lgs \
						$(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

INCLUDES = $(includes) \
include/agmg.hpp \
include/coarse.hpp \
include/defines.hpp \
include/exact.hpp \
include/kernels.hpp \
include/level.hpp \
include/matrix.hpp \
include/oas.hpp \
include/solver.hpp \
include/utils.hpp \
include/vector.hpp

DEPS = $(INCLUDES) \
$(OGSDIR)/ogs.hpp \
$(HDRDIR)/include/types.h

# types of files we are going to construct rules for
.SUFFIXES: .c .cpp

# rule for .c files
.cpp.o: $(DEPS)
	$(CXX) $(CFLAGS) -o $*.o -c $*.cpp $(paths)

# list of objects to be compiled
OBJS = \
./src/agmgLevel.o \
./src/agmgSmoother.o \
./src/exactSolver.o \
./src/kernels.o \
./src/level.o \
./src/matrix.o \
./src/multigrid.o \
./src/oasSolver.o \
./src/parAlmond.o \
./src/pcg.o \
./src/pgmres.o \
./src/solver.o \
./src/SpMV.o \
./src/utils.o \
./src/vector.o \
./src/agmgSetup/agmgSetup.o \
./src/agmgSetup/constructProlongation.o \
./src/agmgSetup/formAggregates.o \
./src/agmgSetup/galerkinProd.o \
./src/agmgSetup/strongGraph.o \
./src/agmgSetup/transpose.o \


all: lib

lib: $(OBJS) $(DEPS)
	ar -cr libparAlmond.a $(OBJS)

clean:
	rm -f libparAlmond.a
	rm -f ./src/*.o
	rm -f ./src/agmgSetup/*.o
