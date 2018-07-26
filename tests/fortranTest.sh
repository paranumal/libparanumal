#!/bin/bash

: ${HOLMES_DIR:=`cd ..; pwd`}
ELLIPTIC_DIR=${HOLMES_DIR}/solvers/elliptic
SRC_DIR=${HOLMES_DIR}/src
: ${OMP_FLAG:=-fopenmp}

cd ${ELLIPTIC_DIR}

# Just to build setupAide.o and fortranInterface.o
echo "Building Holmes ..."
make -j

cd -

echo "Compiling fortranInterfaceTest ..."
mpif77 -static-libgfortran -c fortranInterfaceTest.f

mpicxx ${OMP_FLAG} -o fortranInterfaceTest fortranInterfaceTest.o \
	${SRC_DIR}/setupAide.o ${SRC_DIR}/fortranInterface.o  -lgfortran

echo "Running the interfaceTest ..."
./fortranInterfaceTest
