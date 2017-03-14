CC=mpicc
AR=ar

CPP=-DUSE_NAIVE_BLAS=1 -DNO_NEK_EXITT=1 -DMPI=1 -DUNDERSCORE=1 -DGLOBAL_LONG_LONG  -DPREFIX=jl_

SOURCES = gs.c sort.c sarray_transfer.c sarray_sort.c \
gs_local.c crystal.c comm.c tensor.c fail.c fcrystal.c \
xxt.c sparse_cholesky.c

OBJECTS = $(SOURCES:%.c=obj/%.o)

libgs.a: $(OBJECTS)
	$(AR) -rs $@ $(OBJECTS)

obj/%.o: src/%.c 
	mkdir -p obj
	$(CC) $(CPP) -g -c $< -o $@

clean:
	rm -rf obj
	rm libgs.a
