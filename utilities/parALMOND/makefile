sDir = ./src
iDir = ./include
objDir = ./obj


sources = $(wildcard $(sDir)/*c)
includes = $(wildcard $(iDir)/*h)
objects = $(subst $(sDir)/,$(objDir)/,$(sources:.c=.o))

flags = -I${OCCA_DIR}/include -I$(iDir) 
libs   =  -L${OCCA_DIR}/lib -locca -llapack -lblas 

ifeq ($(shell uname -s),Darwin)

##libs   += -framework Accelerate -framework OpenCL -framework OpenGL -framework GLUT
flags += -D 'OS_OSX'
CC = ${CXX}#mpicxx
# LINUX
else
libs   += -lcuda -lOpenCL
# -lglut -lGL -lGLU -lGLEW
flags  += -D 'OS_LINUX' -D DPWD='"${CURDIR}"'
CC =  mpic++ #${CXX}#mpicxx
endif


# Debug Option
ifeq ($(DEBUG), 1)
flags += -g 
else
# <> For debugging purposes
flags += -O3 -DNDEBUG -g
endif

# Float/Double
ifeq ($(DOUBLE), 1)
flags += -D 'datafloat=double'
else
flags += -D 'datafloat=float'
endif

#flags += -DINS_MPI=$(INS_MPI) -DINS_RENDER=$(INS_RENDER) -DINS_CLUSTER=$(INS_CLUSTER)

all: $(objects) $(includes)
	ar -cr libparALMOND.a $(objects) 


$(objDir)/%.o:$(sDir)/%.c #$(wildcard $(subst $(sDir),$(iDir),$(<:.cpp=.hpp)))
	$(CC) $(flags) -o $@ $(libs) -c $<

clean:
	rm -f libparALMOND.a
	rm -f $(objDir)/*.o
