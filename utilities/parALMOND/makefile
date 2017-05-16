sDir = ./src
iDir = ./include
objDir = ./obj


sources = $(wildcard $(sDir)/*c)
includes = $(wildcard $(iDir)/*h)
objects = $(subst $(sDir)/,$(objDir)/,$(sources:.c=.o))

flags = -fopenmp -I${OCCA_DIR}/include -I$(iDir) -I../../include
libs   =  -L${OCCA_DIR}/lib -locca -llapack -lblas 

flags  += -D DPWD='"${CURDIR}"'
CC =  mpic++ #${CXX}#mpicxx

# Debug Option
ifeq ($(DEBUG), 1)
flags += -g 
else
# <> For debugging purposes
flags += -O3 -DNDEBUG -g
endif

#flags += -DINS_MPI=$(INS_MPI) -DINS_RENDER=$(INS_RENDER) -DINS_CLUSTER=$(INS_CLUSTER)

all: $(objects) $(includes)
	ar -cr libparALMOND.a $(objects) 


$(objDir)/%.o:$(sDir)/%.c #$(wildcard $(subst $(sDir),$(iDir),$(<:.cpp=.hpp)))
	$(CC) $(flags) -o $@ $(libs) -c $<

clean:
	rm -f libparALMOND.a
	rm -f $(objDir)/*.o
