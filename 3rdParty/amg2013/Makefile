#BHEADER**********************************************************************
# Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of HYPRE.  See file COPYRIGHT for details.
#
# HYPRE is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
# $Revision: 2.4 $
#EHEADER**********************************************************************

.SUFFIXES:
.SUFFIXES: .c .o
                                                                                               
srcdir = .

CC = mpicc

AMGDIR= ../../../AMG2013
AR=ar
                                                                                               
CFLAGS = \
 -I$(AMGDIR) \
 -I$(AMGDIR)/utilities\
 -I$(AMGDIR)/struct_mv\
 -I$(AMGDIR)/sstruct_mv\
 -I$(AMGDIR)/IJ_mv\
 -I$(AMGDIR)/seq_mv\
 -I$(AMGDIR)/parcsr_mv\
 -I$(AMGDIR)/parcsr_ls\
 -I$(AMGDIR)/krylov\
 ${INCLUDE_CFLAGS}\
 -DHYPRE_TIMING
                                                                                               
LFLAGS =\
 -L.\
 -L$(AMGDIR)/parcsr_ls\
 -L$(AMGDIR)/parcsr_mv\
 -L$(AMGDIR)/IJ_mv\
 -L$(AMGDIR)/seq_mv\
 -L$(AMGDIR)/sstruct_mv\
 -L$(AMGDIR)/struct_mv\
 -L$(AMGDIR)/krylov\
 -L$(AMGDIR)/utilities\
 -lparcsr_ls\
 -lparcsr_mv\
 -lseq_mv\
 -lsstruct_mv\
 -lIJ_mv\
 -lHYPRE_struct_mv\
 -lkrylov\
 -lHYPRE_utilities\
 ${INCLUDE_LFLAGS}

libs = $(AMGDIR)/parcsr_ls/libparcsr_ls.a\
 $(AMGDIR)/parcsr_mv/libparcsr_mv.a\
 $(AMGDIR)/seq_mv/libseq_mv.a\
 $(AMGDIR)/sstruct_mv/libsstruct_mv.a\
 $(AMGDIR)/IJ_mv/libIJ_mv.a\
 $(AMGDIR)/struct_mv/libHYPRE_struct_mv.a\
 $(AMGDIR)/krylov/libkrylov.a\
 $(AMGDIR)/utilities/libHYPRE_utilities.a

##################################################################
# Targets
##################################################################
                                                                                                                                         
all: libamg2013.a
                                                                                                                                         
clean:
	@rm -f *.o *.a
                                                                                                                                         
veryclean: clean
	@rm -f amg2013
                                                                                                                                         
##################################################################
# Rules
##################################################################
                                                                                                                                         
libamg2013.a: $(AMGDIR) amg2013CoarseSolve.o $(libs)
	@echo  "Linking... "
	${AR} -rcT libamg2013.a amg2013CoarseSolve.o $(libs)                                                                                                                       
                 
$(AMGDIR):
	cd ../../../ && wget https://codesign.llnl.gov/amg2013/amg2013.tgz && \
	tar xvzf amg2013.tgz && rm amg2013.tgz && \
	cd AMG2013  && $(MAKE)

##################################################################
# Generic rules
##################################################################
                                                                                                                                         
.c.o:
	${CC} -o $@ -c ${CFLAGS} $<
                                                               