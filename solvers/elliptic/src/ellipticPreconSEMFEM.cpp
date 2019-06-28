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

#include "ellipticPrecon.hpp"

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
SEMFEMPrecon::SEMFEMPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings) {

}

void SEMFEMPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  if (mesh.elementType==TRIANGLES||mesh.elementType==TETRAHEDRA) {

    dlong Ntotal = mesh.Np*mesh.Nelements;

    // Mr = invDegree.*r
    elliptic.linAlg.amxpy(Ntotal, 1.0, elliptic.o_weight, o_r, 0.0, o_Mr);

    SEMFEMInterpKernel(mesh.Nelements,mesh.o_SEMFEMAnterp,o_Mr,o_rFEM);
    ogsGather(o_GrFEM, o_rFEM, ogsDfloat, ogsAdd, FEMogs);

    parAlmond::Precon(parAlmondHandle, o_GzFEM, o_GrFEM);

    ogsScatter(o_zFEM, o_GzFEM, ogsDfloat, ogsAdd, FEMogs);
    SEMFEMAnterpKernel(mesh.Nelements,mesh.o_SEMFEMAnterp,o_zFEM,o_Mr);

    // Mr = invDegree.*Mr
    elliptic.linAlg.amx(Ntotal, 1.0, elliptic.o_weight, o_Mr);

    ogsGatherScatter(o_Mr, ogsDfloat, ogsAdd, elliptic.ogsMasked);

  } else {

    ogsGather(o_rhsG, o_r, ogsDfloat, ogsAdd, FEMogs);

    dlong N = FEMogs->Ngather;
    elliptic.linAlg.amx(N, 1.0, FEMogs->o_gatherInvDegree, o_rhsG);

    parAlmond::Precon(parAlmondHandle, o_xG, o_rhsG);

    ogsScatter(o_Mr, o_xG, ogsDfloat, ogsAdd, FEMogs);
  }

  if (elliptic.Nmasked)
      elliptic.maskKernel(elliptic.Nmasked, elliptic.o_maskIds, o_Mr);

#if USE_NULL_PROJECTION==1
  if(elliptic.allNeumann) // zero mean of RHS
    elliptic.ZeroMean(o_Mr);
#endif
}