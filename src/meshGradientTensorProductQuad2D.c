#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

// baseline tensor product mesh Gradient for quadrilateral elements
void meshGradientTensorProductQuad2D(mesh2D *mesh,
				     dfloat * q, 
				     dfloat * dqdx, 
				     dfloat * dqdy){
  
  // loop over elements
  for(int e=0;e<mesh->Nelements;++e){
    
    // compute gradient at each node
    for(int j=0;j<mesh->N+1;++j){
      for(int i=0;i<mesh->N+1;++i){
	
	// local node index 
	int n = i + (mesh->N+1)*j;
	
	// load geometric factors
	int gid = mesh->Np*mesh->Nvgeo*e + n;
	float drdx = vgeo[gid + mesh->Np*RXID];
	float drdy = vgeo[gid + mesh->Np*RYID];
	float dsdx = vgeo[gid + mesh->Np*SXID];
	float dsdy = vgeo[gid + mesh->Np*SYID];
	
	// matrix-vector multiplies
	dfloat dqdr = 0, dqds = 0;
	for(int m=0;m<mesh->N+1;++m){
	  dqdr += mesh->D[i*(mesh->N+1) + m]*q[m + j*(mesh->N+1) + e*mesh->Np];
	  dqds += mesh->D[j*(mesh->N+1) + m]*q[i + m*(mesh->N+1) + e*mesh->Np];
	}
	
	// chain rule
	dqdx[n+e*mesh->Np] = drdx*dqdr + dsdx*dqds;
	dqdy[n+e*mesh->Np] = drdy*dqdr + dsdy*dqds;
      }
    }
  }
}

