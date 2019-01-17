To run in Serial mode with 12 MPI processes on a dual hexa core CPU system:

1. make sure setups/setupHex3D.rc requests Serial THREAD MODEL

2.

OCCA_CXXFLAGS='-fstrict-aliasing -funroll-loops -ftree-vectorize -mavx2 -O3'   mpiexec.mpich -n 12  -bind-to core -map-by core  ./ellipticMain setups/setupHex3D.rc