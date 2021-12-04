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

#include "mesh.hpp"
#include "mesh/mesh3D.hpp"

void meshTet3D::ReferenceNodes(int N_){

  N = N_;
  Nfp = ((N+1)*(N+2))/2;
  Np = ((N+1)*(N+2)*(N+3))/6;

  /* Nodal Data */
  r = (dfloat *) malloc(Np*sizeof(dfloat));
  s = (dfloat *) malloc(Np*sizeof(dfloat));
  t = (dfloat *) malloc(Np*sizeof(dfloat));
  NodesTet3D(N, r, s, t);

  faceNodes = (int *) malloc(Nfaces*Nfp*sizeof(int));
  FaceNodesTet3D(N, r, s, t, faceNodes);

  vertexNodes = (int*) calloc(Nverts, sizeof(int));
  VertexNodesTet3D(N, r, s, t, vertexNodes);

  dfloat *V = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  VandermondeTet3D(N, Np, r, s, t, V);

  //Mass matrix
  MM = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  invMM = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  MassMatrixTet3D(Np, V, MM);
  invMassMatrixTet3D(Np, V, invMM);
  free(V);

  //packed D matrices
  D  = (dfloat *) malloc(3*Np*Np*sizeof(dfloat));
  Dr = D + 0*Np*Np;
  Ds = D + 1*Np*Np;
  Dt = D + 2*Np*Np;
  DmatrixTet3D(N, Np, r, s, t, Dr, Ds, Dt);

  //packed D matrices
  DW  = (dfloat *) malloc(3*Np*Np*sizeof(dfloat));
  DWr = DW + 0*Np*Np;
  DWs = DW + 1*Np*Np;
  DWt = DW + 2*Np*Np;
  DWmatrixTet3D(N, Np, r, s, t, MM, DWr, DWs, DWt);

  LIFT = (dfloat *) malloc(Np*Nfaces*Nfp*sizeof(dfloat));
  LIFTmatrixTet3D(N, faceNodes, r, s, t, LIFT);

  sM = (dfloat *) calloc(Np*Nfaces*Nfp,sizeof(dfloat));
  SurfaceMassMatrixTet3D(N, MM, LIFT, sM);

  //packed stiffness matrices
  S = (dfloat*) calloc(6*Np*Np, sizeof(dfloat));
  Srr = S + 0*Np*Np;
  Srs = S + 1*Np*Np;
  Srt = S + 2*Np*Np;
  Sss = S + 3*Np*Np;
  Sst = S + 4*Np*Np;
  Stt = S + 5*Np*Np;
  SmatrixTet3D(N, Dr, Ds, Dt, MM, Srr, Srs, Srt, Sss, Sst, Stt);

  /* Plotting data */
  plotN = N + 3; //enriched interpolation space for plotting
  plotNp = (plotN+1)*(plotN+2)*(plotN+3)/6;

  /* Plotting nodes */
  plotR = (dfloat *) malloc(plotNp*sizeof(dfloat));
  plotS = (dfloat *) malloc(plotNp*sizeof(dfloat));
  plotT = (dfloat *) malloc(plotNp*sizeof(dfloat));
  EquispacedNodesTet3D(plotN, plotR, plotS, plotT);

  plotNelements = plotN*plotN*plotN;
  plotNverts = 4;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));
  EquispacedEToVTet3D(plotN, plotEToV);

  plotInterp = (dfloat *) malloc(Np*plotNp*sizeof(dfloat));
  InterpolationMatrixTet3D(N, Np, r, s, t, plotNp, plotR, plotS, plotT, plotInterp);
}


