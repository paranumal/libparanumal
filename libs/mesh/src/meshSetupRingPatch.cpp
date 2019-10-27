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

//build a new mesh object consisting of the orignal mesh with an
// 1-element overlap with neighboring meshes
mesh_t& mesh_t::SetupRingPatch(){

  //setup the 1-ring halo exchange
  HaloRingSetup();

  //just reuse the current mesh if there are no neighbors
  if (size==1) return *this;

  // single process communicator for new mesh
  MPI_Comm* splitComm = new MPI_Comm;
  MPI_Comm_split(comm, rank, rank, splitComm);

  //make a new occa properties list
  occa::properties *newProps = new occa::properties();
  occaDeviceProperties(device, *newProps);

  mesh_t *mesh=NULL;
  switch(elementType){
  case TRIANGLES:
    if(dim==2)
      mesh = new meshTri2D(device, *splitComm, settings, *newProps);
    else
      mesh = new meshTri3D(device, *splitComm, settings, *newProps);
    break;
  case QUADRILATERALS:
    if(dim==2)
      mesh = new meshQuad2D(device, *splitComm, settings, *newProps);
    else
      mesh = new meshQuad3D(device, *splitComm, settings, *newProps);
    break;
  case TETRAHEDRA:
    mesh = new meshTet3D(device, *splitComm, settings, *newProps);
    break;
  case HEXAHEDRA:
    mesh = new meshHex3D(device, *splitComm, settings, *newProps);
    break;
  }

  //shallow copy of base mesh geometry
  mesh->dim           = dim;
  mesh->Nverts        = Nverts;
  mesh->Nfaces        = Nfaces;
  mesh->NfaceVertices = NfaceVertices;
  mesh->faceVertices  = faceVertices;

  mesh->elementType = elementType;

  mesh->Nnodes = Nnodes; //not really correct, but unused
  mesh->Nelements = Nelements+totalRingElements;
  mesh->NelementsGlobal = Nelements+totalRingElements;

  //populate mesh vertices
  mesh->EX = (dfloat*) calloc(mesh->Nelements*Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nelements*Nverts, sizeof(dfloat));
  if(dim==3)
    mesh->EZ = (dfloat*) calloc(mesh->Nelements*Nverts, sizeof(dfloat));

  memcpy(mesh->EX, EX, Nelements*Nverts*sizeof(dfloat));
  memcpy(mesh->EY, EY, Nelements*Nverts*sizeof(dfloat));
  if(dim==3)
    memcpy(mesh->EZ, EZ, Nelements*Nverts*sizeof(dfloat));

  ringHalo->Exchange(mesh->EX, Nverts, ogs_dfloat);
  ringHalo->Exchange(mesh->EY, Nverts, ogs_dfloat);
  if(dim==3)
    ringHalo->Exchange(mesh->EZ, Nverts, ogs_dfloat);

  mesh->EToV = (hlong*) calloc(mesh->Nelements*Nverts, sizeof(hlong));
  memcpy(mesh->EToV, EToV, Nelements*Nverts*sizeof(hlong));
  ringHalo->Exchange(mesh->EToV, Nverts, ogs_hlong);

  mesh->elementInfo = (hlong*) calloc(mesh->Nelements, sizeof(hlong));
  memcpy(mesh->elementInfo, elementInfo, Nelements*sizeof(hlong));
  ringHalo->Exchange(mesh->elementInfo, 1, ogs_hlong);

  // connect elements using parallel sort
  mesh->ParallelConnect();

  mesh->NboundaryFaces = NboundaryFaces;
  mesh->boundaryInfo = boundaryInfo;

  // element-to-boundary condition type
  mesh->EToB = (int*) calloc(mesh->Nelements*Nfaces, sizeof(int));
  memcpy(mesh->EToB, EToB, Nelements*Nfaces*sizeof(int));
  ringHalo->Exchange(mesh->EToB, Nfaces, ogs_int);

  // correct bcs (replaces unconnected faces with Dirichlet)
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong id = e*Nfaces+f;
      if(mesh->EToE[id]==-1 && mesh->EToB[id]==-1){
        mesh->EToB[id] = 1; // hack to 1 assume Dirichlet
        mesh->EToE[id] = e; // hack to 1 assume Dirichlet
      }
    }
  }

  //Reference Nodes
  mesh->N = N;
  mesh->Np = Np;
  mesh->Nq = Nq;
  mesh->Nfp = Nfp;

  mesh->r = r;
  mesh->s = s;
  mesh->t = t;

  mesh->Dr = Dr;
  mesh->Ds = Ds;
  mesh->Dt = Dt;
  mesh->MM = MM;
  mesh->faceNodes = faceNodes;
  mesh->LIFT = LIFT;
  mesh->interpRaise = interpRaise;
  mesh->interpLower = interpLower;

  mesh->gllz = gllz;
  mesh->gllw = gllw;
  mesh->D = D;
  mesh->gjr = gjr;
  mesh->gjw = gjw;
  mesh->gjI = gjI;
  mesh->gjD = gjD;
  mesh->gjD2 = gjD2;

  mesh->plotNp = plotNp;
  mesh->plotNelements = plotNelements;
  mesh->plotNverts = plotNverts;
  mesh->plotR = plotR;
  mesh->plotS = plotS;
  mesh->plotT = plotT;
  mesh->plotInterp = plotInterp;
  mesh->plotEToV = plotEToV;

  mesh->cubNp = cubNp;
  mesh->cubNq = cubNq;
  mesh->cubNfp = cubNfp;
  mesh->cubr = cubr;
  mesh->cubs = cubs;
  mesh->cubt = cubt;
  mesh->cubw = cubw;
  mesh->cubInterp = cubInterp;
  mesh->cubDrW = cubDrW;
  mesh->cubDsW = cubDsW;
  mesh->cubDtW = cubDtW;
  mesh->cubD = cubD;
  mesh->cubDW = cubDW;
  mesh->cubProject = cubProject;
  mesh->intNfp = intNfp;
  mesh->intInterp = intInterp;
  mesh->intLIFT = intLIFT;

  mesh->NpFEM = NpFEM;
  mesh->NelFEM = NelFEM;
  mesh->rFEM = rFEM;
  mesh->sFEM = sFEM;
  mesh->tFEM = tFEM;
  mesh->SEMFEMInterp = SEMFEMInterp;
  mesh->FEMEToV = FEMEToV;

  mesh->vertexNodes = vertexNodes;

  //Halo
  mesh->HaloSetup();

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
