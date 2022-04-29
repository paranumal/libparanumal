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

//build a new mesh object consisting of the orignal mesh with an
// 1-element overlap with neighboring meshes
mesh_t mesh_t::SetupRingPatch(){

  //setup the 1-ring halo exchange
  HaloRingSetup();

  /*Copy underlying mesh object*/
  mesh_t mesh = *this;

  //just reuse the current mesh if there are no neighbors
  if (size==1) return mesh;

  // single process communicator for new mesh
  mesh.comm = comm.Split(rank, rank);
  mesh.rank = mesh.comm.rank();
  mesh.size = mesh.comm.size();

  mesh.Nelements = Nelements+totalRingElements;
  mesh.NelementsGlobal = Nelements+totalRingElements;

  //populate mesh vertices
  mesh.EX.malloc(mesh.Nelements*Nverts);
  mesh.EY.malloc(mesh.Nelements*Nverts);
  if(dim==3)
    mesh.EZ.malloc(mesh.Nelements*Nverts);

  mesh.EX.copyFrom(EX, Nelements*Nverts);
  mesh.EY.copyFrom(EY, Nelements*Nverts);
  if(dim==3)
    mesh.EZ.copyFrom(EZ, Nelements*Nverts);

  ringHalo.Exchange(mesh.EX, Nverts);
  ringHalo.Exchange(mesh.EY, Nverts);
  if(dim==3)
    ringHalo.Exchange(mesh.EZ, Nverts);

  mesh.EToV.malloc(mesh.Nelements*Nverts);
  mesh.EToV.copyFrom(EToV, Nelements*Nverts);
  ringHalo.Exchange(mesh.EToV, Nverts);

  mesh.elementInfo.malloc(mesh.Nelements);
  mesh.elementInfo.copyFrom(elementInfo, Nelements);
  ringHalo.Exchange(mesh.elementInfo, 1);

  // connect elements using parallel sort
  mesh.Connect();

  mesh.NboundaryFaces = NboundaryFaces;
  mesh.boundaryInfo = boundaryInfo;

  // element-to-boundary condition type
  mesh.EToB.malloc(mesh.Nelements*Nfaces);
  mesh.EToB.copyFrom(EToB, Nelements*Nfaces);
  ringHalo.Exchange(mesh.EToB, Nfaces);

  // correct bcs (replaces unconnected faces with Dirichlet)
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong id = e*Nfaces+f;
      if(mesh.EToE[id]==-1 && mesh.EToB[id]==-1){
        mesh.EToB[id] = 1; // hack to 1 assume Dirichlet
        mesh.EToE[id] = e; // hack to 1 assume Dirichlet
      }
    }
  }

  //Halo
  mesh.HaloSetup();

  // connect face vertices
  mesh.ConnectFaceVertices();

  // connect face nodes (find trace indices)
  mesh.ConnectFaceNodes();

  // make global indexing
  mesh.ConnectNodes();

  // compute physical (x,y) locations of the element nodes
  mesh.PhysicalNodes();

  // compute geometric factors
  mesh.GeometricFactors();

  // compute surface geofacs
  mesh.SurfaceGeometricFactors();

  // label local/global gather elements
  mesh.GatherScatterSetup();

  return mesh;
}

} //namespace libp
