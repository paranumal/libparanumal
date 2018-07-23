sDir = ./src
iDir = ./include
objDir = ./

include ${OCCA_DIR}/scripts/Makefile

sources = $(wildcard $(sDir)/*c)
includes = $(wildcard $(iDir)/*h)
objects = $(subst $(sDir)/,$(objDir)/,$(sources:.c=.o))
deps = $(includes) \
../../include/mesh.h \
../../include/ogs_t.h \
../../include/hgs_t.h \
../../include/parAlmond.h

flags = -DOCCA_VERSION_1_0 -I${OCCA_DIR}/include -I$(iDir) -I../../include
libs   =  -L${OCCA_DIR}/lib -locca -llapack -lblas

flags  += -D DPWD='"${CURDIR}"'
CC =  mpic++
#flags += -fopenmp

# Debug Option
ifeq ($(DEBUG), 1)
flags += -g
else
# <> For debugging purposes
flags += -O3 -DNDEBUG  -fopenmp
endif

#flags += -DINS_MPI=$(INS_MPI) -DINS_RENDER=$(INS_RENDER) -DINS_CLUSTER=$(INS_CLUSTER)

all: lib

lib: $(objects) $(deps)
	ar -cr libparALMOND.a $(objects)

$(objDir)/%.o:$(sDir)/%.c $(deps) #$(wildcard $(subst $(sDir),$(iDir),$(<:.cpp=.hpp)))
	$(CC) $(flags) -o $@ $(libs) -c $< $(paths)

clean:
	rm -f libparALMOND.a
	rm -f $(objDir)/*.o
