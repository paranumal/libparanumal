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
                                                                                               
srcdir = $(shell pwd)

CC = mpicc

AMGDIR= ../../../AMG2013
AR=ar
             


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
                                                                                                                                         
libamg2013.a: $(AMGDIR)
	@echo  "Linking... " 
	${AR} -x $(AMGDIR)/parcsr_ls/libparcsr_ls.a; \
	${AR} -x $(AMGDIR)/parcsr_mv/libparcsr_mv.a; \
	${AR} -x $(AMGDIR)/seq_mv/libseq_mv.a;       \
	${AR} -x $(AMGDIR)/sstruct_mv/libsstruct_mv.a; \
	${AR} -x $(AMGDIR)/IJ_mv/libIJ_mv.a; \
	${AR} -x $(AMGDIR)/struct_mv/libHYPRE_struct_mv.a; \
	${AR} -x $(AMGDIR)/krylov/libkrylov.a; \
	${AR} -x $(AMGDIR)/utilities/libHYPRE_utilities.a; \
	${AR} rc libamg2013.a *.o; \
	rm *.o

$(AMGDIR):
	cd ../../../ && wget https://codesign.llnl.gov/amg2013/amg2013.tgz && \
	tar xvzf amg2013.tgz && rm amg2013.tgz; \
	cp $(srcdir)/Makefile.include AMG2013/; \
	cd AMG2013  && $(MAKE)

##################################################################
# Generic rules
##################################################################
                                                                                                                                         
.c.o:
	${CC} -o $@ -c ${CFLAGS} $<
                                                               