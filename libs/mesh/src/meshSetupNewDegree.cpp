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

#include "core.hpp"
#include "mesh.hpp"
#include "mesh2D.hpp"
#include "mesh3D.hpp"

//build a new mesh object from another with a different degree.
mesh_t& mesh_t::SetupNewDegree(int Nf){

  //just reuse the current mesh if the degree isnt changing.
  if (Nf==N) return *this;

  //make a new occa properties list
  occa::properties *newProps = new occa::properties();
  occaDeviceProperties(device, *newProps);

  mesh_t *mesh=NULL;
  switch(elementType){
  case TRIANGLES:
    if(dim==2)
      mesh = new meshTri2D(device, comm, settings, *newProps);
    else
      mesh = new meshTri3D(device, comm, settings, *newProps);
    break;
  case QUADRILATERALS:
    if(dim==2)
      mesh = new meshQuad2D(device, comm, settings, *newProps);
    else
      mesh = new meshQuad3D(device, comm, settings, *newProps);
    break;
  case TETRAHEDRA:
    mesh = new meshTet3D(device, comm, settings, *newProps);
    break;
  case HEXAHEDRA:
    mesh = new meshHex3D(device, comm, settings, *newProps);
    break;
  }

  //shallow copy of base mesh geometry
  mesh->dim           = dim;
  mesh->Nverts        = Nverts;
  mesh->Nfaces        = Nfaces;
  mesh->NfaceVertices = NfaceVertices;
  mesh->faceVertices  = faceVertices;

  mesh->elementType = elementType;

  mesh->Nnodes = Nnodes;
  mesh->EX = EX; // coordinates of vertices for each element
  mesh->EY = EY;
  mesh->EZ = EZ;

  mesh->Nelements = Nelements;
  mesh->NelementsGlobal = NelementsGlobal;
  mesh->EToV = EToV; // element-to-vertex connectivity
  mesh->EToE = EToE; // element-to-element connectivity
  mesh->EToF = EToF; // element-to-(local)face connectivity
  mesh->EToP = EToP; // element-to-partition/process connectivity
  mesh->EToB = EToB; // element-to-boundary condition type

  mesh->elementInfo = elementInfo;

  mesh->NboundaryFaces = NboundaryFaces;
  mesh->boundaryInfo = boundaryInfo;

  mesh->halo = halo;
  mesh->NinternalElements = NinternalElements;
  mesh->NhaloElements = NhaloElements;
  mesh->totalHaloPairs = totalHaloPairs;
  mesh->internalElementIds = internalElementIds;
  mesh->haloElementIds = haloElementIds;
  mesh->o_internalElementIds = o_internalElementIds;
  mesh->o_haloElementIds     = o_haloElementIds;

  mesh->ogs = ogs;
  mesh->globalIds = globalIds;

  mesh->NglobalGatherElements = NglobalGatherElements;
  mesh->globalGatherElementList = globalGatherElementList;
  mesh->o_globalGatherElementList = o_globalGatherElementList;

  mesh->NlocalGatherElements = NlocalGatherElements;
  mesh->localGatherElementList = localGatherElementList;
  mesh->o_localGatherElementList = o_localGatherElementList;

  mesh->defaultStream = defaultStream;

  // load reference (r,s) element nodes
  mesh->ReferenceNodes(Nf);

  // compute physical (x,y) locations of the element nodes
  mesh->PhysicalNodes();

  // compute geometric factors
  mesh->GeometricFactors();

  // connect face nodes (find trace indices)
  mesh->ConnectFaceNodes();

  // compute surface geofacs
  mesh->SurfaceGeometricFactors();

  // make a global indexing
  mesh->ParallelConnectNodes();

  // make an ogs operator and label local/global gather elements
  mesh->ParallelGatherScatterSetup();

  mesh->OccaSetup();

  return *mesh;
}
