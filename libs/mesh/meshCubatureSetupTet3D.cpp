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

#include "mesh.hpp"

namespace libp {

void mesh_t::CubatureSetupTet3D(){

  /* Cubature data */
  cubN = 2*N; //cubature order
  CubatureNodesTet3D(cubN, cubNp, cubr, cubs, cubt, cubw);

  InterpolationMatrixTet3D(N, r, s, t, cubr, cubs, cubt, cubInterp);

  //cubature project cubProject = M^{-1} * cubInterp^T
  // Defined such that cubProject * cubW * cubInterp = Identity
  CubaturePmatrixTet3D(N, r, s, t, cubr, cubs, cubt, cubProject);

  //cubature derivates matrices, cubD: differentiate on cubature nodes
  // we dont use cubD on Tris/Tets  so skip computing

  // Instead, it's cheaper to:
  // make weak cubature derivatives cubPDT = cubProject * cubD^T
  CubatureWeakDmatricesTet3D(N, r, s, t,
                             cubr, cubs, cubt,
                             cubPDT);
  cubPDrT = cubPDT + 0*cubNp*Np;
  cubPDsT = cubPDT + 1*cubNp*Np;
  cubPDtT = cubPDT + 2*cubNp*Np;

  // Surface cubature nodes
  CubatureNodesTri2D(cubN, intNfp, intr, ints, intw);
  cubNfp = intNfp;

  CubatureSurfaceMatricesTet3D(N, r, s, t, faceNodes,
                               intr, ints, intw,
                               intInterp, intLIFT);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

  // build transposes (we hold matrices as column major on device)
  memory<dfloat> cubProjectT(cubNp*Np);
  memory<dfloat> cubInterpT(cubNp*Np);
  linAlg_t::matrixTranspose(cubNp, Np, cubInterp, Np, cubInterpT, cubNp);
  linAlg_t::matrixTranspose(Np, cubNp, cubProject, cubNp, cubProjectT, Np);

  //pre-multiply cubProject by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubProjectT[m+n*Np] *= cubw[n];
    }
  }

  memory<dfloat> cubPDTT(3*cubNp*Np);
  memory<dfloat> cubPDrTT = cubPDTT + 0*cubNp*Np;
  memory<dfloat> cubPDsTT = cubPDTT + 1*cubNp*Np;
  memory<dfloat> cubPDtTT = cubPDTT + 2*cubNp*Np;
  linAlg_t::matrixTranspose(Np, cubNp, cubPDrT, cubNp, cubPDrTT, Np);
  linAlg_t::matrixTranspose(Np, cubNp, cubPDsT, cubNp, cubPDsTT, Np);
  linAlg_t::matrixTranspose(Np, cubNp, cubPDtT, cubNp, cubPDtTT, Np);

  //pre-multiply cubPDT by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubPDrTT[m+n*Np] *= cubw[n];
      cubPDsTT[m+n*Np] *= cubw[n];
      cubPDtTT[m+n*Np] *= cubw[n];
    }
  }

  // build surface integration matrix transposes
  memory<dfloat> intLIFTT(Np*Nfaces*intNfp);
  memory<dfloat> intInterpT(Nfp*Nfaces*intNfp);
  linAlg_t::matrixTranspose(Np, Nfaces*intNfp, intLIFT, Nfaces*intNfp, intLIFTT, Np);
  linAlg_t::matrixTranspose(Nfaces*intNfp, Nfp, intInterp, Nfp, intInterpT, Nfaces*intNfp);

  o_cubInterp  = platform.malloc<dfloat>(Np*cubNp, cubInterpT);
  o_cubProject = platform.malloc<dfloat>(Np*cubNp, cubProjectT);

  o_cubPDT = platform.malloc<dfloat>(3*Np*cubNp, cubPDTT);

  o_intInterp = platform.malloc<dfloat>(Nfp*Nfaces*intNfp, intInterpT);
  o_intLIFT   = platform.malloc<dfloat>(Np*Nfaces*intNfp, intLIFTT);

  o_cubwJ = o_wJ;
  o_cubvgeo = o_vgeo;
  o_cubggeo = o_ggeo;
  o_cubsgeo = o_sgeo;
}

} //namespace libp
