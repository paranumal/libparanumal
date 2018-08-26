#!/bin/bash

: ${HOLMES_DIR:=`cd ..; pwd`}
: ${OCCA_DIR:="${HOLMES_DIR}/../occa"}
: ${OMP_FLAG:=-fopenmp}

ELLIPTIC_DIR=${HOLMES_DIR}/solvers/elliptic
SRC_DIR=${HOLMES_DIR}/src

cd ${OCCA_DIR}; make -j; cd -

# Just to build setupAide.o and fortranInterface.o
echo "Building Holmes ..."
cd ${ELLIPTIC_DIR}; make -j; cd -

echo "Compiling fortranInterfaceTest ..."
mpif77 -static-libgfortran -c fortranInterfaceTest.f

mpicxx ${OMP_FLAG} -o fortranInterfaceTest fortranInterfaceTest.o \
       ${SRC_DIR}/setupAide.o ${SRC_DIR}/fortranInterface.o ${SRC_DIR}/fortranMeshSetup.o \
       ${SRC_DIR}/meshPartitionStatistics.o ${SRC_DIR}/meshConnectBoundary.o \
       ${SRC_DIR}/meshLoadReferenceNodesHex3D.o ${SRC_DIR}/meshHaloSetup.o \
       ${SRC_DIR}/meshConnectFaceNodes3D.o ${SRC_DIR}/meshParallelConnectNodes.o \
       ${SRC_DIR}/meshParallelConnectOpt.o ${SRC_DIR}/meshConnect.o \
       ${SRC_DIR}/meshSurfaceGeometricFactorsHex3D.o ${SRC_DIR}/readArray.o \
       ${SRC_DIR}/meshHaloExchange.o ${SRC_DIR}/mysort.o -lgfortran

echo "Running the interfaceTest ..."
./fortranInterfaceTest
