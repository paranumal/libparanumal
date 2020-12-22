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

void meshHex3D::ReferenceNodes(int N_){

  N = N_;
  Nq = N+1;
  Nfp = Nq*Nq;
  Np = Nq*Nq*Nq;

  /* Nodal Data */
  r = (dfloat *) malloc(Np*sizeof(dfloat));
  s = (dfloat *) malloc(Np*sizeof(dfloat));
  t = (dfloat *) malloc(Np*sizeof(dfloat));
  NodesHex3D(N, r, s, t);

  faceNodes = (int *) malloc(Nfaces*Nfp*sizeof(int));
  FaceNodesHex3D(N, r, s, t, faceNodes);

  vertexNodes = (int*) calloc(Nverts, sizeof(int));
  VertexNodesHex3D(N, r, s, t, vertexNodes);

  //GLL quadrature
  dfloat *gllz = (dfloat *) malloc((N+1)*sizeof(dfloat));
  w = (dfloat *) malloc((N+1)*sizeof(dfloat));
  JacobiGLL(N, gllz, w);

  //Lumped Mass matrix
  MM    = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  invMM = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  LumpedMassMatrixHex3D(N, w, MM);
  invLumpedMassMatrixHex3D(N, w, invMM);

  // D matrix
  D = (dfloat *) malloc(Nq*Nq*sizeof(dfloat));
  Dmatrix1D(N, Nq, gllz, Nq, gllz, D);

  /* Plotting data */
  plotN = N_ + 3; //enriched interpolation space for plotting
  plotNq = plotN + 1;
  plotNp = plotNq*plotNq*plotNq;

  /* Plotting nodes */
  plotR = (dfloat *) malloc(plotNp*sizeof(dfloat));
  plotS = (dfloat *) malloc(plotNp*sizeof(dfloat));
  plotT = (dfloat *) malloc(plotNp*sizeof(dfloat));
  EquispacedNodesHex3D(plotN, plotR, plotS, plotT);

  plotNelements = 6*plotN*plotN*plotN;
  plotNverts = 4;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));
  EquispacedEToVHex3D(plotN, plotEToV);

  dfloat *plot1D = (dfloat *) malloc(plotNq*sizeof(dfloat));
  EquispacedNodes1D(plotN, plot1D);

  plotInterp = (dfloat *) malloc(Nq*plotNq*sizeof(dfloat));
  InterpolationMatrix1D(N, Nq, gllz, plotNq, plot1D, plotInterp);

  free(gllz);
  free(plot1D);
}

