/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"


void elliptic_t::SetupNewCoefficient(elliptic_t &elliptic){

  mesh_t &meshC = elliptic.mesh; 

  const dlong Nall   = meshC.Np *(meshC.Nelements+meshC.totalHaloPairs);     
  // currently scalar coefficients are supported
  elliptic.coeff     = (dfloat *) calloc(2*Nall, sizeof(dfloat)); 
  elliptic.o_coeff   = meshC.device.malloc(2*Nall*sizeof(dfloat), elliptic.coeff);

  dfloat *cR; occa::memory o_cR; 
 
  // printf("Interpolating Coefficient from  %d  to %d\n", mesh.N, meshC.N); 
  meshC.BuildInterpolation(&cR, o_cR, mesh.N, meshC.N);

  //build kernels
    occa::properties kernelInfo = elliptic.props;

    // set kernel name suffix
    char *suffix;
    if(meshC.elementType==TRIANGLES)
      suffix = strdup("Tri2D");
    else if(meshC.elementType==QUADRILATERALS)
      suffix = strdup("Quad2D");
    else if(meshC.elementType==TETRAHEDRA)
      suffix = strdup("Tet3D");
    else if(meshC.elementType==HEXAHEDRA)
      suffix = strdup("Hex3D");

    char fileName[BUFSIZ], kernelName[BUFSIZ];

    kernelInfo["defines/" "p_NqFine"]    = mesh.N+1;
    kernelInfo["defines/" "p_NqCoarse"]  = meshC.N+1;

    kernelInfo["defines/" "p_NpFine"]    = mesh.Np;
    kernelInfo["defines/" "p_NpCoarse"]  = meshC.Np;

    int NblockVFine   = 512/mesh.Np;
    int NblockVCoarse = 512/meshC.Np;
    kernelInfo["defines/" "p_NblockVFine"]= NblockVFine;
    kernelInfo["defines/" "p_NblockVCoarse"]= NblockVCoarse;

    sprintf(fileName, DELLIPTIC "/okl/ellipticCoarsen%s.okl", suffix);
    sprintf(kernelName, "ellipticCoarsenCoefficient%s", suffix);
    
    elliptic.coarsenCoefficientKernel = buildKernel(meshC.device, fileName, kernelName, kernelInfo, meshC.comm);

    elliptic.coarsenCoefficientKernel(meshC.Nelements, o_cR, o_coeff, elliptic.o_coeff); 

    elliptic.o_coeff.copyTo(elliptic.coeff);

}
