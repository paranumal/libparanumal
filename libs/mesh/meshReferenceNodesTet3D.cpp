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

void mesh_t::ReferenceNodesTet3D(){

  Nfp = ((N+1)*(N+2))/2;
  Np = ((N+1)*(N+2)*(N+3))/6;

  /* Nodal Data */
  NodesTet3D(N, r, s, t);
  FaceNodesTet3D(N, r, s, t, faceNodes);
  VertexNodesTet3D(N, r, s, t, vertexNodes);

  memory<dfloat> V;
  VandermondeTet3D(N, r, s, t, V);

  //Mass matrix
  MassMatrixTet3D(Np, V, MM);
  invMassMatrixTet3D(Np, V, invMM);
  o_MM = platform.malloc<dfloat>(MM); //MM is symmetric

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_MM = o_MM;
    pfloat_invMM = invMM;
  } else {
    memory<pfloat> pfloat_MM(Np*Np);
    pfloat_invMM.malloc(Np*Np);
    for(int n=0;n<Np*Np;++n){
      pfloat_MM[n] = MM[n];
      pfloat_invMM[n] = invMM[n];
    }
    o_pfloat_MM = platform.malloc<pfloat>(pfloat_MM);
  }


  //packed D matrices
  DmatrixTet3D(N, r, s, t, D);
  Dr = D + 0*Np*Np;
  Ds = D + 1*Np*Np;
  Dt = D + 2*Np*Np;

  memory<dfloat> DT(3*Np*Np);
  memory<dfloat> DrT = DT + 0*Np*Np;
  memory<dfloat> DsT = DT + 1*Np*Np;
  memory<dfloat> DtT = DT + 2*Np*Np;
  linAlg_t::matrixTranspose(Np, Np, Dr, Np, DrT, Np);
  linAlg_t::matrixTranspose(Np, Np, Ds, Np, DsT, Np);
  linAlg_t::matrixTranspose(Np, Np, Dt, Np, DtT, Np);
  o_D = platform.malloc<dfloat>(DT);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_D = o_D;
  } else {
    memory<pfloat> pfloat_DT(Np*Np*dim);
    for(int n=0;n<Np*Np*dim;++n) {
      pfloat_DT[n] = DT[n];
    }
    o_pfloat_D = platform.malloc<pfloat>(pfloat_DT);
  }


  LIFTmatrixTet3D(N, faceNodes, r, s, t, LIFT);
  SurfaceMassMatrixTet3D(N, MM, LIFT, sM);

  memory<dfloat> LIFTT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, LIFT, Nfp*Nfaces, LIFTT, Np);

  memory<dfloat> sMT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, sM, Nfp*Nfaces, sMT, Np);

  o_sM = platform.malloc<dfloat>(sMT);
  o_LIFT = platform.malloc<dfloat>(LIFTT);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_LIFT = o_LIFT;
  } else {
    memory<pfloat> pfloat_LIFTT(Np*Np*dim);
    for(int n=0;n<Np*Nfaces*Nfp;++n) {
      pfloat_LIFTT[n] = LIFTT[n];
    }
    o_pfloat_LIFT = platform.malloc<pfloat>(pfloat_LIFTT);
  }

  //packed stiffness matrices
  SmatrixTet3D(N, Dr, Ds, Dt, MM, S);
  Srr = S + 0*Np*Np;
  Srs = S + 1*Np*Np;
  Srt = S + 2*Np*Np;
  Sss = S + 3*Np*Np;
  Sst = S + 4*Np*Np;
  Stt = S + 5*Np*Np;

  memory<dfloat> ST(6*Np*Np);
  memory<dfloat> SrrT = ST + 0*Np*Np;
  memory<dfloat> SrsT = ST + 1*Np*Np;
  memory<dfloat> SrtT = ST + 2*Np*Np;
  memory<dfloat> SssT = ST + 3*Np*Np;
  memory<dfloat> SstT = ST + 4*Np*Np;
  memory<dfloat> SttT = ST + 5*Np*Np;
  linAlg_t::matrixTranspose(Np, Np, Srr, Np, SrrT, Np);
  linAlg_t::matrixTranspose(Np, Np, Srs, Np, SrsT, Np);
  linAlg_t::matrixTranspose(Np, Np, Srt, Np, SrtT, Np);
  linAlg_t::matrixTranspose(Np, Np, Sss, Np, SssT, Np);
  linAlg_t::matrixTranspose(Np, Np, Sst, Np, SstT, Np);
  linAlg_t::matrixTranspose(Np, Np, Stt, Np, SttT, Np);

  o_S = platform.malloc<dfloat>(ST);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_S = o_S;
  } else {
    memory<pfloat> pfloat_ST(Np*Np*( (dim)*(dim+1)/2 ));
    for(int n=0;n<(Np*Np*dim*(dim+1))/2;++n) {
      pfloat_ST[n] = ST[n];
    }
    o_pfloat_S = platform.malloc<pfloat>(pfloat_ST);
  }

  // packed strong form stiffness matrices
  //packed stiffness matrices
  StrongSmatrixTet3D(N, Dr, Ds, Dt, MM, strongS);
  strongSrr = strongS + 0*Np*Np;
  strongSrs = strongS + 1*Np*Np;
  strongSrt = strongS + 2*Np*Np;
  strongSss = strongS + 3*Np*Np;
  strongSst = strongS + 4*Np*Np;
  strongStt = strongS + 5*Np*Np;

  memory<dfloat> strongST(6*Np*Np);
  memory<dfloat> strongSrrT = strongST + 0*Np*Np;
  memory<dfloat> strongSrsT = strongST + 1*Np*Np;
  memory<dfloat> strongSrtT = strongST + 2*Np*Np;
  memory<dfloat> strongSssT = strongST + 3*Np*Np;
  memory<dfloat> strongSstT = strongST + 4*Np*Np;
  memory<dfloat> strongSttT = strongST + 5*Np*Np;
  linAlg_t::matrixTranspose(Np, Np, strongSrr, Np, strongSrrT, Np);
  linAlg_t::matrixTranspose(Np, Np, strongSrs, Np, strongSrsT, Np);
  linAlg_t::matrixTranspose(Np, Np, strongSrt, Np, strongSrtT, Np);
  linAlg_t::matrixTranspose(Np, Np, strongSss, Np, strongSssT, Np);
  linAlg_t::matrixTranspose(Np, Np, strongSst, Np, strongSstT, Np);
  linAlg_t::matrixTranspose(Np, Np, strongStt, Np, strongSttT, Np);

  o_strongS = platform.malloc<dfloat>(strongST);



  
  /* Plotting data */
  plotN = N + 1; //enriched interpolation space for plotting
  plotNp = (plotN+1)*(plotN+2)*(plotN+3)/6;

  /* Plotting nodes */
  EquispacedNodesTet3D(plotN, plotR, plotS, plotT);

  plotNelements = plotN*plotN*plotN;
  plotNverts = 4;
  EquispacedEToVTet3D(plotN, plotEToV);
  InterpolationMatrixTet3D(N, r, s, t, plotR, plotS, plotT, plotInterp);

  props["defines/" "p_N"]= N;
  props["defines/" "p_Np"]= Np;
  props["defines/" "p_Nfp"]= Nfp;
  props["defines/" "p_NfacesNfp"]= Nfp*Nfaces;
}

} //namespace libp
