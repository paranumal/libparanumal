#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshGradient2D(mesh2D *mesh,
		    dfloat *q,
		    dfloat *dqdx,
		    dfloat *dqdy
		    ){
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
    
    for(iint n=0;n<mesh->Np;++n){
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
