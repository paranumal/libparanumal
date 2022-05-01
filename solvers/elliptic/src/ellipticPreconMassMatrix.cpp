/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
  LIBP_ABORT("MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead.",
             mesh.elementType!=Mesh::TRIANGLES && mesh.elementType!=Mesh::TETRAHEDRA);

  LIBP_ABORT("MASSMATRIX preconditioner is unavailble when lambda=0.",
             elliptic.lambda==0);

  o_invMM = elliptic.platform.malloc<dfloat>(mesh.invMM);

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  int blockMax = 256;
  if (elliptic.platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    blockJacobiKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "blockJacobi", kernelInfo);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    dlong Ntotal = elliptic.ogsMasked.Ngather + elliptic.gHalo.Nhalo;
    o_rtmp = elliptic.platform.malloc<dfloat>(Ntotal);
    o_MrL  = elliptic.platform.malloc<dfloat>(mesh.Np*mesh.Nelements);

    partialBlockJacobiKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticPreconBlockJacobi.okl",
                                     "partialBlockJacobi", kernelInfo);
  }
}

void MassMatrixPrecon::Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr) {
  dfloat invLambda = 1./elliptic.lambda;

  linAlg_t& linAlg = elliptic.platform.linAlg();

  if (elliptic.disc_c0) {//C0

    // rtmp = invDegree.*r
    linAlg.amxpy(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_r, 0.0, o_rtmp);

    elliptic.gHalo.ExchangeStart(o_rtmp, 1);

    if(mesh.NlocalGatherElements/2)
      partialBlockJacobiKernel(mesh.NlocalGatherElements/2,
                               mesh.o_localGatherElementList,
                               elliptic.o_GlobalToLocal,
                               invLambda, mesh.o_vgeo, o_invMM,
                               o_rtmp, o_MrL);

    // finalize halo exchange
    elliptic.gHalo.ExchangeFinish(o_rtmp, 1);

    if(mesh.NglobalGatherElements)
      partialBlockJacobiKernel(mesh.NglobalGatherElements,
                               mesh.o_globalGatherElementList,
                               elliptic.o_GlobalToLocal,
                               invLambda, mesh.o_vgeo, o_invMM,
                               o_rtmp, o_MrL);

    //gather result to Aq
    elliptic.ogsMasked.GatherStart(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);

    if((mesh.NlocalGatherElements+1)/2){
      partialBlockJacobiKernel((mesh.NlocalGatherElements+1)/2,
                               mesh.o_localGatherElementList+mesh.NlocalGatherElements/2,
                               elliptic.o_GlobalToLocal,
                               invLambda, mesh.o_vgeo, o_invMM,
                               o_rtmp, o_MrL);
    }

    elliptic.ogsMasked.GatherFinish(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);

    // Mr = invDegree.*Mr
    linAlg.amx(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_Mr);

  } else {
    //IPDG
    blockJacobiKernel(mesh.Nelements, invLambda, mesh.o_vgeo, o_invMM, o_r, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}
