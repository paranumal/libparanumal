CC=mpicc #-std=c99 --pedantic
CFLAGS+=-DMPI 
#CFLAGS+= -acc -ta=nvidia -Minfo=accel
CFLAGS+=-DPREFIX=jl_
CFLAGS+=-DNO_NEK_EXITT
CFLAGS+=-DGLOBAL_LONG
LDFLAGS+=-lm #-acc
SRCLOC=../../src
#CFLAGS+=-DPRINT_MALLOCS=1

CFLAGS+=-DUSE_NAIVE_BLAS
#CFLAGS+=-DUSE_CBLAS
#LDFLAGS+=-lcblas

#CFLAGS+=-DAMG_DUMP
CFLAGS+=-DGS_TIMING -DGS_BARRIER

CFLAGS+=-O0 -g
#CFLAGS+=-O3 #-march=native

#CFLAGS+=-W #-Wall -Wno-unused-function -Wno-unused-parameter
#CFLAGS+=-Minform=warn

CCCMD=$(CC) $(G) $(CFLAGS)
LINKCMD=$(CC) $(G) $(LDFLAGS)
#RLINKCMD = $(LD) -r
all: gslib_test gs_alltoall gs_allreduce
clean: ; @$(RM) $(TESTS) *.o *.s

cmds: ; @echo CC = $(CCCMD); echo LINK = $(LINKCMD);

%.o: %.c ; @echo CC $<; $(CCCMD) -c $< -o $@
%.s: %.c ; @echo CC -S $<; $(CCCMD) -S $< -o $@
objects: $(OBJECTS) ;

GS_OBJECTS=$(SRCLOC)/gs.o $(SRCLOC)/sort.o $(SRCLOC)/sarray_transfer.o $(SRCLOC)/sarray_sort.o \
           $(SRCLOC)/gs_local.o $(SRCLOC)/fail.o $(SRCLOC)/crystal.o $(SRCLOC)/comm.o $(SRCLOC)/tensor.o

gslib_test: gslib_test.o $(GS_OBJECTS);		@echo LINK $@; $(LINKCMD) $^ -lm -o $@
gs_alltoall: gs_alltoall.o $(GS_OBJECTS);		@echo LINK $@; $(LINKCMD) $^ -o $@
gs_allreduce: gs_allreduce.o $(GS_OBJECTS);		@echo LINK $@; $(LINKCMD) $^ -o $@

