#include "mesh3D.h"

mesh3D *meshSetupHex3D(char *filename, int N){

  // read chunk of elements
  mesh3D *mesh = meshParallelReaderHex3D(filename);

  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition3D(mesh); 
  
  // print out connectivity statistics
  meshPartitionStatistics3D(mesh);

  // connect elements using parallel sort
  meshParallelConnect3D(mesh);

  // connect elements to boundary faces
  meshConnectBoundary3D(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesHex3D(mesh);

  // compute geometric factors
  meshGeometricFactorsHex3D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup3D(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);
  
  // compute surface geofacs (including halo)
  meshSurfaceGeometricFactorsHex3D(mesh);
  
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
