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

#include "subcell.hpp"
#include "subcell2D.hpp"
// #include "subcell3D.hpp"

// subcell_t:: subcell_t(mesh_t &_mesh, settings_t& _settings):
//     mesh(_mesh),
//     comm(_mesh.comm),
//     device(_mesh.device),
//     settings(_settings),
//     props(_mesh.props){};

subcell_t:: subcell_t(solver_t &_solver):
    mesh(_solver.mesh),
    comm(_solver.comm),
    device(_solver.device),
    settings(_solver.settings),
    props(_solver.props),
    solver(_solver){};


  // subcell_t& subcell_t::Setup(mesh_t& _mesh, settings_t& _settings){
  subcell_t& subcell_t::Setup(solver_t& _solver){

  subcell_t *subcell=NULL; 

  switch(_solver.mesh.elementType){
  case TRIANGLES:
    if(_solver.mesh.dim==2)
      subcell = new subcellTri2D(_solver);
      // subcell = new subcellTri2D(_mesh, _settings);
    else
      // subcell = new subcellTri3D(device, comm, settings, props);
    break;
  // case QUADRILATERALS:
  //   if(_mesh.dim==2)
  //     subcell = new subcellQuad2D(_mesh, _settings); 
  //   else
  //     // mesh = new meshQuad3D(device, comm, settings, props);
  //   break;
  // case TETRAHEDRA:
    // mesh = new meshTet3D(device, comm, settings, props);
    // break;
  // case HEXAHEDRA:
    // mesh = new meshHex3D(device, comm, settings, props);
    // break;
  }

  // subcell->settings.getSetting("SUBCELL NUMBER", subcell->N);
  // printf("%d\n", subcell->N);
  
  // Create Minor
  subcell->CreateMinorGrid(); 

  // Create local connection in subcell
  subcell->LocalConnect(); 

  // Global connection of FV elements
  subcell->GlobalConnect(); 

  // Compute Geometric factors for FV dicretization
  subcell->GeometricFactors(); 

  // Setup projection and recontruction operators
  subcell->SetupOperators(); 

  // Setup detector i.e. Skyline for now
  subcell->SetupDetector(); 

  // Create device info
  subcell->OccaSetup();
  
  #if 0
 
  // connect elements using parallel sort
  mesh->ParallelConnect();

  // print out connectivity statistics
  mesh->PrintPartitionStatistics();

  // connect elements to boundary faces
  mesh->ConnectBoundary();

  // load reference (r,s) element nodes
  mesh->ReferenceNodes(N);

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
#endif

  return *subcell;
}

//
void subcell_t::LeastSquaresFit(int _N, dfloat *_LSF){

dfloat *temp = (dfloat *)malloc(2*_N*sizeof(dfloat)); 
for(int n=0; n<_N; n++){
  const dfloat logmode = log10(n+1); 
  temp[2*n + 0] = logmode; 
  temp[2*n + 1] = 1.0; 
}

matrixPInverse(_N, 2, temp);

for(int n=0; n<_N; n++){
  _LSF[n] = temp[n];
}
free(temp);

}

// 
void subcell_t::BaseLineDecay(int _N, dfloat *_BLD){

dfloat bsum = 0.0; 
for(int j=1; j<_N+1; j++)
    bsum +=1.0/pow(j, 2*_N); 

bsum = 1.0/sqrt(bsum); 
  
BLD[0] = 0.0; 
// baseline decay: AK: squared !!!!
for(int n=1; n<_N+1; n++){
    const dfloat bdecay = bsum*1.0/(pow(n,_N));
    _BLD[n] = bdecay*bdecay;
 }
}


// void subcell_t::InvVandermonde1D(int _N, dfloat *invV1D){
  
//   int _Nq = _N+1;

//   dfloat *r1D = (dfloat*) malloc(_Nq*sizeof(dfloat));
  
//   mesh.JacobiGLL(_N, r1D); //Gauss-Legendre-Lobatto nodes

//   mesh.Vandermonde1D(_N, _Nq, r1D, invV1D);

//   matrixInverse(_Nq, invV1D); 

// }


subcell_t::~subcell_t() {
  
}


