# holmes
It's elementary.

# installation
git clone https://github.com/tcew/holmes

# OCCA dependency
git clone https://github.com/libocca/occa

# build OCCA
cd occa
export OCCA_DIR=`pwd`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OCCA_DIR/lib
make -j
cd ../

# build holmes SEM elliptic example
cd holmes
cd examples/ellipticQuad2D
make -j
cd ../../

# run SEM elliptic example: 2 MPI processes, degree 3, on sample mesh
mpiexec -n 2 ./examples/ellipticQuad2D/ellipticMainQuad2D meshes/cavityQuadH02.msh 3
