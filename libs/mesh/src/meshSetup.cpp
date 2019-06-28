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

mesh_t& mesh_t::Setup(occa::device& device, MPI_Comm& comm,
                     settings_t& settings, occa::properties& props){

  string fileName;
  int N, dim, elementType;

  settings.getSetting("MESH FILE", fileName);
  settings.getSetting("POLYNOMIAL DEGREE", N);
  settings.getSetting("ELEMENT TYPE", elementType);
  settings.getSetting("MESH DIMENSION", dim);

  mesh_t *mesh=NULL;
  switch(elementType){
  case TRIANGLES:
    if(dim==2)
      mesh = new meshTri2D(device, comm, settings, props);
    else
      mesh = new meshTri3D(device, comm, settings, props);
    break;
  case QUADRILATERALS:
    if(dim==2)
      mesh = new meshQuad2D(device, comm, settings, props);
    else
      mesh = new meshQuad3D(device, comm, settings, props);
    break;
  case TETRAHEDRA:
    mesh = new meshTet3D(device, comm, settings, props);
    break;
  case HEXAHEDRA:
    mesh = new meshHex3D(device, comm, settings, props);
    break;
  }

  mesh->elementType = elementType;

  if (settings.compareSetting("MESH FILE","BOX")) {
    //build a box mesh
    mesh->SetupBox();
  } else {
    // read chunk of elements from file
    mesh->ParallelReader(fileName.c_str());

    // partition elements using Morton ordering & parallel sort
    mesh->GeometricPartition();
  }

  // connect elements using parallel sort
  mesh->ParallelConnect();

  // print out connectivity statistics
  mesh->PrintPartitionStatistics();

  // connect elements to boundary faces
  mesh->ConnectBoundary();

  // load reference (r,s) element nodes
  mesh->LoadReferenceNodes(N);

  // set up halo exchange info for MPI (do before connect face nodes)
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
