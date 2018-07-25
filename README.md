# holmes
It's elementary.

## Clone: Holmes
git clone https://github.com/tcew/holmes

## OCCA dependency (currently OCCA 1.0 forked by Noel Chalmers)
git clone https://github.com/noelchalmers/occa

## Build OCCA
cd occa

export OCCA_DIR=`pwd`

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib

make -j

cd ../

## build holmes elliptic example
cd holmes

cd solvers/elliptic

make -j

## run elliptic example: 2 MPI processes with settings in supplied setup file
mpiexec -n 2 ./ellipticMain setups/setupQuad2D.rc
