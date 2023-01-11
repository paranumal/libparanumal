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

  LIFTmatrixTet3D(N, faceNodes, r, s, t, LIFT);
  SurfaceMassMatrixTet3D(N, MM, LIFT, sM);

  memory<dfloat> LIFTT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, LIFT, Nfp*Nfaces, LIFTT, Np);

  memory<dfloat> sMT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, sM, Nfp*Nfaces, sMT, Np);

  o_sM = platform.malloc<dfloat>(sMT);
  o_LIFT = platform.malloc<dfloat>(LIFTT);

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

  /* Plotting data */
  plotN = N + 3; //enriched interpolation space for plotting
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
