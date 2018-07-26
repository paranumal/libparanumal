#!/bin/bash

: ${HOLMES_DIR:=`cd ..; pwd`}
ELLIPTIC_DIR=${HOLMES_DIR}/solvers/elliptic
SRC_DIR=${HOLMES_DIR}/src

cd ${ELLIPTIC_DIR}

# Just to build setupAide.o and fortranInterface.o
echo "Building Holmes ..."
make -j

cd -

echo "Compiling fortranInterfaceTest ..."
mpif77 -static-libgfortran -c fortranInterfaceTest.f

mpicxx -fopenmp -o fortranInterfaceTest fortranInterfaceTest.o \
	${HOLMES_DIR}/src/setupAide.o ${HOLMES_DIR}/src/fortranInterface.o  -lgfortran

echo "Running the interfaceTest ..."
./fortranInterfaceTest
