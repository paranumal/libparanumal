#include "mesh3D.h"

mesh_t *meshSetupQuad3D(char *filename, int N, dfloat sphereRadius,char *mode){

  // read chunk of elements
  mesh_t *mesh = meshParallelReaderQuad3D(filename);

  // set sphere radius (will be used later in building physical nodes)
  mesh->sphereRadius = sphereRadius;

  // partition elements using Morton ordering & parallel sort
  //meshGeometricPartition3D(mesh); // need to double check this

  // connect elements using parallel sort
  meshParallelConnect(mesh);


  
  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);

#if 0
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      printf("%d ", mesh->EToB[e*mesh->Nfaces+f]);
    }
    printf("\n");
  }
#endif
  
  // load reference (r,s) element nodes
  void meshLoadReferenceNodesQuad3D(mesh_t *mesh, int N);
  meshLoadReferenceNodesQuad3D(mesh, N);

  // compute physical (x,y,z) locations of the element nodes
  //meshPhysicalNodesQuad3D(mesh);
  //meshSphericalNodesQuad3D(mesh);
  meshEquiSphericalNodesQuad3D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  // compute geometric factors
  meshGeometricFactorsQuad3D(mesh);
  
  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);

  for(int n=0;n<mesh->Nfp*mesh->Nelements*mesh->Nfaces;++n){
    if(mesh->vmapM[n]==mesh->vmapP[n]){
      printf("node %d matches self \n");
    }
  }
      
  
  // compute surface geofacs
  meshSurfaceGeometricFactorsQuad3D(mesh);
  
  // global nodes
  meshParallelConnectNodes(mesh);
     
  return mesh;
}
