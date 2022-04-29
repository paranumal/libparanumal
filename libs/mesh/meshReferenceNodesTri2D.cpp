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

void mesh_t::ReferenceNodesTri2D(){

  Nfp = N+1;
  Np = (N+1)*(N+2)/2;

  /* Nodal Data */
  NodesTri2D(N, r, s);
  FaceNodesTri2D(N, r, s, faceNodes);
  VertexNodesTri2D(N, r, s, vertexNodes);

  memory<dfloat> V;
  VandermondeTri2D(N, r, s, V);

  //Mass matrix
  MassMatrixTri2D(Np, V, MM);
  invMassMatrixTri2D(Np, V, invMM);
  o_MM = platform.malloc<dfloat>(MM); //MM is symmetric

  //packed D matrices
  DmatrixTri2D(N, r, s, D);
  Dr = D + 0*Np*Np;
  Ds = D + 1*Np*Np;

  memory<dfloat> DT(2*Np*Np);
  memory<dfloat> DrT = DT + 0*Np*Np;
  memory<dfloat> DsT = DT + 1*Np*Np;
  linAlg_t::matrixTranspose(Np, Np, Dr, Np, DrT, Np);
  linAlg_t::matrixTranspose(Np, Np, Ds, Np, DsT, Np);
  o_D = platform.malloc<dfloat>(DT);

  LIFTmatrixTri2D(N, faceNodes, r, s, LIFT);
  SurfaceMassMatrixTri2D(N, MM, LIFT, sM);

  memory<dfloat> LIFTT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, LIFT, Nfp*Nfaces, LIFTT, Np);

  memory<dfloat> sMT(Np*Nfaces*Nfp);
  linAlg_t::matrixTranspose(Np, Nfp*Nfaces, sM, Nfp*Nfaces, sMT, Np);

  o_sM = platform.malloc<dfloat>(sMT);
  o_LIFT = platform.malloc<dfloat>(LIFTT);

  //packed stiffness matrices
  SmatrixTri2D(N, Dr, Ds, MM, S);
  Srr = S + 0*Np*Np;
  Srs = S + 1*Np*Np;
  Sss = S + 2*Np*Np;

  memory<dfloat> ST(3*Np*Np);
  memory<dfloat> SrrT = ST + 0*Np*Np;
  memory<dfloat> SrsT = ST + 1*Np*Np;
  memory<dfloat> SssT = ST + 2*Np*Np;
  linAlg_t::matrixTranspose(Np, Np, Srr, Np, SrrT, Np);
  linAlg_t::matrixTranspose(Np, Np, Srs, Np, SrsT, Np);
  linAlg_t::matrixTranspose(Np, Np, Sss, Np, SssT, Np);

  o_S = platform.malloc<dfloat>(ST);

  /* Plotting data */
  plotN = N + 3; //enriched interpolation space for plotting
  plotNp = (plotN+1)*(plotN+2)/2;

  /* Plotting nodes */
  EquispacedNodesTri2D(plotN, plotR, plotS);

  plotNelements = plotN*plotN;
  plotNverts = 3;
  EquispacedEToVTri2D(plotN, plotEToV);
  InterpolationMatrixTri2D(N, r, s, plotR, plotS, plotInterp);

  props["defines/" "p_N"]= N;
  props["defines/" "p_Np"]= Np;
  props["defines/" "p_Nfp"]= Nfp;
  props["defines/" "p_NfacesNfp"]= Nfp*Nfaces;
}

} //namespace libp
