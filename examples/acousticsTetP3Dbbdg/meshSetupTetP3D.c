#include "mesh3D.h"

mesh3D *meshSetupTetP3D(char *filename, int N){

  // read chunk of elements
  mesh3D *mesh = meshParallelReaderTet3D(filename);

  // partition elements using Morton ordering & parallel sort
  meshGeometricPartition3D(mesh);
  
  // connect elements using parallel sort
  meshParallelConnect(mesh);
  
  // print out connectivity statistics
  meshPartitionStatistics(mesh);

  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesTetP3D(mesh, N);

  // compute physical (x,y,z) locations of the element nodes
  meshPhysicalNodesTetP3D(mesh);

  // compute geometric factors
  meshGeometricFactorsTet3D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetupP(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodesP3D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTetP3D(mesh);


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
  dfloat rkc[5] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0}; 

  mesh->Nrk = Nrk;
  memcpy(mesh->rka, rka, Nrk*sizeof(dfloat));
  memcpy(mesh->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(mesh->rkc, rkc, Nrk*sizeof(dfloat));

#if 0
  for(iint eM=0;eM<mesh->Nelements;++eM){
    for(iint fM=0;fM<mesh->Nfaces;++fM){

      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      iint fP = mesh->EToF[eM*mesh->Nfaces+fM];

      if(eP<mesh->Nelements){
	if(eP<0){ eP = eM; fP = fM; }
      
	dfloat sc = 1;
	if(eM==eP) sc = -1;
      
	iint baseM = (eM*mesh->Nfaces + fM)*mesh->Nsgeo;
	iint baseP = (eP*mesh->Nfaces + fP)*mesh->Nsgeo;

	dfloat avenx = mesh->sgeo[baseM+NXID] + sc*mesh->sgeo[baseP+NXID];
	dfloat aveny = mesh->sgeo[baseM+NYID] + sc*mesh->sgeo[baseP+NYID];
	dfloat avenz = mesh->sgeo[baseM+NZID] + sc*mesh->sgeo[baseP+NZID];
	dfloat dsJ   = mesh->sgeo[baseM+SJID] - mesh->sgeo[baseP+SJID];
	dfloat aven  = norm(avenx,aveny,avenz);

	if(aven>1e-4)
	  printf("(%d,%d=>%d,%d) aven = %g,%g,%g", eM,fM, eP,fP, avenx,aveny,avenz);

	if(fabs(dsJ)>1e-4)
	  printf("sJ mismatch %g\n", dsJ);
	
	for(iint n=0;n<mesh->Nfp;++n){
	  iint id = n + fM*mesh->Nfp + eM*mesh->Nfp*mesh->Nfaces;
	  iint idM = mesh->vmapM[id];
	  iint idP = mesh->vmapP[id];

	  dfloat d = norm(mesh->x[idM]-mesh->x[idP],mesh->y[idM]-mesh->y[idP],mesh->z[idM]-mesh->z[idP]);

	  if(d>1e-4)
	    printf(", %g\n", d);
	}
      }
    }
  }
#endif  
  return mesh;
}
