#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshGradientQuad2D(mesh2D *mesh,
			dfloat *q,
			dfloat *dqdx,
			dfloat *dqdy
			){
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    
    for(iint n=0;n<mesh->Np;++n){
      iint gid = e*mesh->Nvgeo*mesh->Np + n;
      dfloat drdx = mesh->vgeo[gid + mesh->Np*RXID];
      dfloat drdy = mesh->vgeo[gid + mesh->Np*RYID];
      dfloat dsdx = mesh->vgeo[gid + mesh->Np*SXID];
      dfloat dsdy = mesh->vgeo[gid + mesh->Np*SYID];

      dfloat dqdr = 0, dqds = 0;
      
      for(iint m=0;m<mesh->Np;++m){
	
	dqdr += mesh->Dr[n*mesh->Np + m]*q[m + e*mesh->Np];
	dqds += mesh->Ds[n*mesh->Np + m]*q[m + e*mesh->Np];
      }
      
      dqdx[n+e*mesh->Np] = drdx*dqdr + dsdx*dqds;
      dqdy[n+e*mesh->Np] = drdy*dqdr + dsdy*dqds;
    }
  }
}
