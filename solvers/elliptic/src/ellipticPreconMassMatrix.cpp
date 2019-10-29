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

// Inverse Mass Matrix preconditioner
MassMatrixPrecon::MassMatrixPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings) {

  //sanity checking
  if (mesh.elementType!=TRIANGLES && mesh.elementType!=TETRAHEDRA )
    LIBP_ABORT(string("MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead."));

  if (elliptic.lambda==0)
    LIBP_ABORT(string("MASSMATRIX preconditioner is unavailble when lambda=0."));

  dlong Ntotal = mesh.Np*mesh.Nelements;
  o_rtmp = mesh.device.malloc(Ntotal*sizeof(dfloat));
  o_invMM = mesh.device.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.invMM);

  // OCCA build stuff
  occa::properties kernelInfo = elliptic.props; //copy base occa properties

  int NblockV = mymax(1,512/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    blockJacobiKernel = buildKernel(mesh.device, DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "blockJacobi", kernelInfo, mesh.comm);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    partialBlockJacobiKernel = buildKernel(mesh.device, DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "partialBlockJacobi", kernelInfo, mesh.comm);
  }
}

void MassMatrixPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {
  dfloat invLambda = 1./elliptic.lambda;

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {

    blockJacobiKernel(mesh.Nelements, invLambda, mesh.o_vgeo, o_invMM, o_r, o_Mr);

  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {

    dlong Ntotal = mesh.Np*mesh.Nelements;

    // rtmp = invDegree.*r
    elliptic.linAlg.amxpy(Ntotal, 1.0, elliptic.o_weight, o_r, 0.0, o_rtmp);

    if(mesh.NglobalGatherElements)
      partialBlockJacobiKernel(mesh.NglobalGatherElements,
                               mesh.o_globalGatherElementList,
                               invLambda, mesh.o_vgeo, o_invMM, o_rtmp, o_Mr);

    elliptic.ogsMasked->GatherScatterStart(o_Mr, ogs_dfloat, ogs_add, ogs_sym);

    if(mesh.NlocalGatherElements)
      partialBlockJacobiKernel(mesh.NlocalGatherElements,
                               mesh.o_localGatherElementList,
                               invLambda, mesh.o_vgeo, o_invMM, o_rtmp, o_Mr);

    elliptic.ogsMasked->GatherScatterFinish(o_Mr, ogs_dfloat, ogs_add, ogs_sym);

    // Mr = invDegree.*Mr
    elliptic.linAlg.amx(Ntotal, 1.0, elliptic.o_weight, o_Mr);

    //post-mask
    if (elliptic.Nmasked)
      elliptic.maskKernel(elliptic.Nmasked, elliptic.o_maskIds, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MassMatrixPrecon::~MassMatrixPrecon(){
  blockJacobiKernel.free();
  partialBlockJacobiKernel.free();
}