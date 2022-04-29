#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

define LIB_HELP_MSG

Libraries makefile targets:

	 make all (default)
	 make lib{}
	 make clean
	 make clean-kernels
	 make realclean
	 make info
	 make help

Usage:

make all
	 Build core, ogs, parAdogs, mesh, linAlg, timeStepper, linearSovler, and parAlmond libraries.
make {lib}
	 Build only {lib} library.
make clean
	 Clean all libraries, and object files.
make clean-kernels
	 In addition to "make clean", also cleans the cached OCCA kernels.
make realclean
	 In addition to "make clean-kernels", also clean 3rd party libraries.
make info
	 List directories and compiler flags in use.
make help
	 Display this help message.

Can use "make verbose=true" for verbose output.

endef

ifeq (,$(filter core linAlg linearSolver parAdogs mesh ogs parAlmond timeStepper \
								clean clean-libs clean-kernels realclean info help, $(MAKECMDGOALS)))
ifneq (,$(MAKECMDGOALS))
$(error ${LIB_HELP_MSG})
endif
endif

ifndef LIBP_MAKETOP_LOADED
ifeq (,$(wildcard ../make.top))
$(error cannot locate ${PWD}/../make.top)
else
include ../make.top
endif
endif

#libraries
LIBCORE_DIR=${LIBP_LIBS_DIR}/core
LIBLINALG_DIR=${LIBP_LIBS_DIR}/linAlg
LIBLINEARSOLVER_DIR=${LIBP_LIBS_DIR}/linearSolver
LIBPARADOGS_DIR=${LIBP_LIBS_DIR}/parAdogs
LIBMESH_DIR=${LIBP_LIBS_DIR}/mesh
LIBOGS_DIR=${LIBP_LIBS_DIR}/ogs
LIBPARALMOND_DIR=${LIBP_LIBS_DIR}/parAlmond
LIBTIMESTEPPER_DIR=${LIBP_LIBS_DIR}/timeStepper

#includes
INCLUDES=${LIBP_INCLUDES}

#defines
LIBCORE_DEFINES=-DLIBP_DIR='"${LIBP_DIR}"'
LIBLINALG_DEFINES=-DLINALG_DIR='"${LIBLINALG_DIR}"'
LIBLINEARSOLVER_DEFINES=-DLINEARSOLVER_DIR='"${LIBLINEARSOLVER_DIR}"'
LIBPARADOGS_DEFINES=-DPARADOGS_DIR='"${PARADOGS_DIR}"'
LIBMESH_DEFINES=-DMESH_DIR='"${LIBMESH_DIR}"'
LIBOGS_DEFINES=-DLIBP_DIR='"${LIBP_DIR}"'  -DOGS_DIR='"${LIBOGS_DIR}"'
LIBPARALMOND_DEFINES=-DPARALMOND_DIR='"${LIBPARALMOND_DIR}"'
LIBTIMESTEPPER_DEFINES=-DTIMESTEPPER_DIR='"${LIBTIMESTEPPER_DIR}"'

ifeq (true,${gpu-aware-mpi})
  LIBOGS_DEFINES+= -DGPU_AWARE_MPI
endif

#.cpp compilation flags
LIB_CXXFLAGS=${LIBP_DEFINES} ${LIBP_CXXFLAGS} ${INCLUDES}

#object dependancies
LIB_DEPS=$(wildcard $(LIBP_INCLUDE_DIR)/*.h)  \
				 $(wildcard ${LIBP_INCLUDE_DIR}/*.hpp)

LIBMESH_DEPS=$(wildcard ${LIBP_INCLUDE_DIR}/mesh/*.h)   \
						 $(wildcard ${LIBP_INCLUDE_DIR}/mesh/*.hpp) \
						 ${LIB_DEPS}

LIBOGS_DEPS=$(wildcard ${LIBP_INCLUDE_DIR}/ogs/*.h)   \
						$(wildcard ${LIBP_INCLUDE_DIR}/ogs/*.hpp) \
						${LIB_DEPS}

LIBPARADOGS_DEPS=$(wildcard ${LIBP_INCLUDE_DIR}/parAdogs/*.hpp) \
								 ${LIB_DEPS}

LIBPARALMOND_DEPS=$(wildcard ${LIBP_INCLUDE_DIR}/parAlmond/*.h)   \
				  				$(wildcard ${LIBP_INCLUDE_DIR}/parAlmond/*.hpp) \
								  ${LIB_DEPS}

LIBCORE_SRC =$(wildcard core/*.cpp)
LIBLINALG_SRC =$(wildcard linAlg/*.cpp)
LIBLINEARSOLVER_SRC =$(wildcard linearSolver/*.cpp)
LIBPARADOGS_SRC =$(wildcard parAdogs/*.cpp)
LIBMESH_SRC =$(wildcard mesh/*.cpp)
LIBOGS_SRC =$(wildcard ogs/*.cpp)
LIBPARALMOND_SRC =$(wildcard parAlmond/*.cpp)
LIBTIMESTEPPER_SRC =$(wildcard timeStepper/*.cpp)

LIBCORE_OBJS=$(LIBCORE_SRC:.cpp=.o)
LIBLINALG_OBJS=$(LIBLINALG_SRC:.cpp=.o)
LIBLINEARSOLVER_OBJS=$(LIBLINEARSOLVER_SRC:.cpp=.o)
LIBPARADOGS_OBJS=$(LIBPARADOGS_SRC:.cpp=.o)
LIBMESH_OBJS=$(LIBMESH_SRC:.cpp=.o)
LIBOGS_OBJS=$(LIBOGS_SRC:.cpp=.o)
LIBPARALMOND_OBJS=$(LIBPARALMOND_SRC:.cpp=.o)
LIBTIMESTEPPER_OBJS=$(LIBTIMESTEPPER_SRC:.cpp=.o)

.PHONY: all core linAlg linearSolver mesh ogs parAdogs parAlmond timeStepper \
				clean realclean silentUpdate-core             \
				silentUpdate-linAlg silentUpdate-linearSolver \
				silentUpdate-ogs silentUpdate-mesh            \
				silentUpdate-parAdogs                         \
				silentUpdate-parAlmond silentUpdate-timeStepper

all: core linAlg linearSolver parAdogs mesh ogs parAlmond timeStepper

core: libcore.a silentUpdate-core
linAlg: liblinAlg.a silentUpdate-linAlg
linearSolver: liblinearSolver.a silentUpdate-linearSolver
mesh: libmesh.a silentUpdate-mesh
ogs: libogs.a silentUpdate-ogs
parAdogs: libparAdogs.a silentUpdate-parAdogs
parAlmond: libparAlmond.a silentUpdate-parAlmond
timeStepper: libtimeStepper.a silentUpdate-timeStepper

libcore.a: $(LIBCORE_OBJS)
ifneq (,${verbose})
	ar -cr libcore.a $(LIBCORE_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libcore.a $(LIBCORE_OBJS)
endif

libogs.a: $(LIBOGS_OBJS) | libcore.a
ifneq (,${verbose})
	ar -cr libogs.a $(LIBOGS_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libogs.a $(LIBOGS_OBJS)
endif

liblinAlg.a: $(LIBLINALG_OBJS) | libcore.a
ifneq (,${verbose})
	ar -cr liblinAlg.a $(LIBLINALG_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr liblinAlg.a $(LIBLINALG_OBJS)
endif

libparAdogs.a: $(LIBPARADOGS_OBJS) | liblinAlg.a libogs.a
ifneq (,${verbose})
	ar -cr libparAdogs.a $(LIBPARADOGS_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libparAdogs.a $(LIBPARADOGS_OBJS)
endif

libmesh.a: $(LIBMESH_OBJS) | libparAdogs.a
ifneq (,${verbose})
	ar -cr libmesh.a $(LIBMESH_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libmesh.a $(LIBMESH_OBJS)
endif

libtimeStepper.a: $(LIBTIMESTEPPER_OBJS) | libmesh.a
ifneq (,${verbose})
	ar -cr libtimeStepper.a $(LIBTIMESTEPPER_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libtimeStepper.a $(LIBTIMESTEPPER_OBJS)
endif

liblinearSolver.a: $(LIBLINEARSOLVER_OBJS) | libcore.a
ifneq (,${verbose})
	ar -cr liblinearSolver.a $(LIBLINEARSOLVER_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr liblinearSolver.a $(LIBLINEARSOLVER_OBJS)
endif

libparAlmond.a: $(LIBPARALMOND_OBJS) | liblinearSolver.a libogs.a
ifneq (,${verbose})
	ar -cr libparAlmond.a $(LIBPARALMOND_OBJS)
else
	@printf "%b" "$(LIB_COLOR)Building library $(@F)$(NO_COLOR)\n";
	@ar -cr libparAlmond.a $(LIBPARALMOND_OBJS)
endif

silentUpdate-core:
	@true

silentUpdate-linAlg:
	@true

silentUpdate-linearSolver:
	@true

silentUpdate-ogs:
	@true

silentUpdate-mesh:
	@true

silentUpdate-parAdogs:
	@true

silentUpdate-parAlmond:
	@true

silentUpdate-timeStepper:
	@true

${OCCA_DIR}/lib/libocca.so:
	${MAKE} -C ${OCCA_DIR}

# rule for .cpp files
core/%.o: core/%.cpp $(LIB_DEPS) ${OCCA_DIR}/lib/libocca.so
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBCORE_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBCORE_DEFINES} $(LIB_CXXFLAGS)
endif

linAlg/%.o: linAlg/%.cpp $(LIB_DEPS) | libcore.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBLINALG_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBLINALG_DEFINES} $(LIB_CXXFLAGS)
endif

parAdogs/%.o: parAdogs/%.cpp $(LIB_DEPS) | libogs.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBPARADOGS_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBPARADOGS_DEFINES} $(LIB_CXXFLAGS)
endif

linearSolver/%.o: linearSolver/%.cpp $(LIB_DEPS) | libcore.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBLINEARSOLVER_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBLINEARSOLVER_DEFINES} $(LIB_CXXFLAGS)
endif

mesh/%.o: mesh/%.cpp $(LIBMESH_DEPS) | libparAdogs.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBMESH_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBMESH_DEFINES} $(LIB_CXXFLAGS)
endif

ogs/%.o: ogs/%.cpp $(LIBOGS_DEPS) | liblinAlg.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBOGS_DEFINES} -I${GS_DIR}/src $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBOGS_DEFINES} -I${GS_DIR}/src $(LIB_CXXFLAGS)
endif

parAlmond/%.o: parAlmond/%.cpp $(LIBPARALMOND_DEPS) | liblinearSolver.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBPARALMOND_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBPARALMOND_DEFINES} $(LIB_CXXFLAGS)
endif

timeStepper/%.o: timeStepper/%.cpp $(LIB_DEPS) | libmesh.a
ifneq (,${verbose})
	$(LIBP_CXX) -o $@ -c $< ${LIBTIMESTEPPER_DEFINES} $(LIB_CXXFLAGS)
else
	@printf "%b" "$(OBJ_COLOR)Compiling $(@F)$(NO_COLOR)\n";
	@$(LIBP_CXX) -o $@ -c $< ${LIBTIMESTEPPER_DEFINES} $(LIB_CXXFLAGS)
endif

#cleanup
clean:
	rm -f core/*.o linAlg/*.o linearSolver/*.o mesh/*.o ogs/*.o parAdogs/*.o parAlmond/*.o timeStepper/*.o *.a

clean-kernels: clean
	rm -rf ${LIBP_DIR}/.occa/

realclean: clean
	${MAKE} -C ${OCCA_DIR} clean

help:
	$(info $(value LIB_HELP_MSG))
	@true

info:
	$(info OCCA_DIR  = $(OCCA_DIR))
	$(info LIBP_DIR  = $(LIBP_DIR))
	$(info LIBP_ARCH = $(LIBP_ARCH))
	$(info CXXFLAGS  = $(LIB_CXXFLAGS))
	$(info LIBS      = $(LIBS))
	@true
