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

  o_invMM = elliptic.platform.malloc(mesh.Np*mesh.Np*sizeof(dfloat), mesh.invMM);

  // OCCA build stuff
  occa::properties kernelInfo = elliptic.mesh.props; //copy base occa properties

  int blockMax = 256;
  if (elliptic.platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = mymax(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    blockJacobiKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "blockJacobi", kernelInfo);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    dlong Ntotal = elliptic.ogsMasked->Ngather + elliptic.ogsMasked->NgatherHalo;
    o_rtmp = elliptic.platform.malloc(Ntotal*sizeof(dfloat));
    o_MrL  = elliptic.platform.malloc(mesh.Np*mesh.Nelements*sizeof(dfloat));

    partialBlockJacobiKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "partialBlockJacobi", kernelInfo);
  }
}

void MassMatrixPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {
  dfloat invLambda = 1./elliptic.lambda;

  if (elliptic.disc_c0) {//C0

    // rtmp = invDegree.*r
    elliptic.linAlg.amxpy(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_r, 0.0, o_rtmp);

    elliptic.ogsMasked->GatheredHaloExchangeStart(o_rtmp, 1, ogs_dfloat);

    if(mesh.NlocalGatherElements)
      partialBlockJacobiKernel(mesh.NlocalGatherElements,
                               mesh.o_localGatherElementList,
                               elliptic.ogsMasked->o_GlobalToLocal,
                               invLambda, mesh.o_vgeo, o_invMM,
                               o_rtmp, o_MrL);

    elliptic.ogsMasked->GatheredHaloExchangeFinish(o_rtmp, 1, ogs_dfloat);

    if(mesh.NglobalGatherElements)
      partialBlockJacobiKernel(mesh.NglobalGatherElements,
                               mesh.o_globalGatherElementList,
                               elliptic.ogsMasked->o_GlobalToLocal,
                               invLambda, mesh.o_vgeo, o_invMM,
                               o_rtmp, o_MrL);

    //gather result to Aq
    elliptic.ogsMasked->Gather(o_Mr, o_MrL, ogs_dfloat, ogs_add, ogs_trans);

    // Mr = invDegree.*Mr
    elliptic.linAlg.amx(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_Mr);

  } else {
    //IPDG
    blockJacobiKernel(mesh.Nelements, invLambda, mesh.o_vgeo, o_invMM, o_r, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MassMatrixPrecon::~MassMatrixPrecon(){
  blockJacobiKernel.free();
  partialBlockJacobiKernel.free();
}