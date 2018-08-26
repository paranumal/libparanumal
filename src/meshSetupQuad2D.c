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

#include "mesh2D.h"

mesh2D *meshSetupQuad2D(char *filename, int N){

  // read chunk of elements
  mesh2D *mesh = meshParallelReaderQuad2D(filename);
  
  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition2D(mesh);

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);
  
  // load reference (r,s) element nodes
  meshLoadReferenceNodesQuad2D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesQuad2D(mesh);

  // compute geometric factors
  meshGeometricFactorsQuad2D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);
  
  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsQuad2D(mesh);
  
  // global nodes
  meshParallelConnectNodes(mesh);
  
  // initialize LSERK4 time stepping coefficients
  int Nrk = 5;

  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  dfloat rkc[6] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0,
		   1.}; 

  mesh->Nrk = Nrk;
  memcpy(mesh->rka, rka, Nrk*sizeof(dfloat));
  memcpy(mesh->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(mesh->rkc, rkc, (Nrk+1)*sizeof(dfloat));
    
  return mesh;
}
