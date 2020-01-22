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
#include "mesh2D.hpp"
#include "mesh3D.hpp"

void meshTri3D::ReferenceNodes(int N_){
  mesh_t *mesh_p = (mesh_t*) this;
  meshTri2D* trimesh = (meshTri2D*) mesh_p;
  trimesh->meshTri2D::ReferenceNodes(N);
}

void meshTri2D::ReferenceNodes(int N_){

  N = N_;
  Nfp = N+1;
  Np = (N+1)*(N+2)/2;

  /* Nodal Data */
  r = (dfloat *) malloc(Np*sizeof(dfloat));
  s = (dfloat *) malloc(Np*sizeof(dfloat));
  NodesTri2D(N, r, s);

  faceNodes = (int *) malloc(Nfaces*Nfp*sizeof(int));
  FaceNodesTri2D(N, r, s, faceNodes);

  vertexNodes = (int*) calloc(Nverts, sizeof(int));
  VertexNodesTri2D(N, r, s, vertexNodes);

  dfloat *V = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  VandermondeTri2D(N, Np, r, s, V);

  //Mass matrix
  MM = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  MassMatrixTri2D(Np, V, MM);
  free(V);

  //D matrices
  Dr = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  Ds = (dfloat *) malloc(Np*Np*sizeof(dfloat));
  DmatrixTri2D(N, Np, r, s, Dr, Ds);

  LIFT = (dfloat *) malloc(Np*Nfaces*Nfp*sizeof(dfloat));
  LIFTmatrixTri2D(N, faceNodes, r, s, LIFT);

  /* Plotting data */
  int plotN = N + 3; //enriched interpolation space for plotting
  plotNp = (plotN+1)*(plotN+2)/2;

  /* Plotting nodes */
  plotR = (dfloat *) malloc(plotNp*sizeof(dfloat));
  plotS = (dfloat *) malloc(plotNp*sizeof(dfloat));
  EquispacedNodesTri2D(plotN, plotR, plotS);

  plotNelements = plotN*plotN;
  plotNverts = 3;
  plotEToV = (int*) malloc(plotNelements*plotNverts*sizeof(int));
  EquispacedEToVTri2D(plotN, plotEToV);

  plotInterp = (dfloat *) malloc(Np*plotNp*sizeof(dfloat));
  InterpolationMatrixTri2D(N, Np, r, s, plotNp, plotR, plotS, plotInterp);
}
