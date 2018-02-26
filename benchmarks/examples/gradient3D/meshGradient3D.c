#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshGradient3D(mesh3D *mesh,
		    dfloat *q,
		    dfloat *dqdx,
		    dfloat *dqdy,
		    dfloat *dqdz
		    ){
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    
    dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
    dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
    dfloat drdz = mesh->vgeo[e*mesh->Nvgeo + RZID];
    dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
    dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
    dfloat dsdz = mesh->vgeo[e*mesh->Nvgeo + SZID];
    dfloat dtdx = mesh->vgeo[e*mesh->Nvgeo + TXID];
    dfloat dtdy = mesh->vgeo[e*mesh->Nvgeo + TYID];
    dfloat dtdz = mesh->vgeo[e*mesh->Nvgeo + TZID];
    
    for(iint n=0;n<mesh->Np;++n){
      dfloat dqdr = 0, dqds = 0, dqdt = 0;
      
      for(iint m=0;m<mesh->Np;++m){
	
	dqdr += mesh->Dr[n*mesh->Np + m]*q[m + e*mesh->Np];
	dqds += mesh->Ds[n*mesh->Np + m]*q[m + e*mesh->Np];
	dqdt += mesh->Dt[n*mesh->Np + m]*q[m + e*mesh->Np];
      }
      
      dqdx[n+e*mesh->Np] = drdx*dqdr + dsdx*dqds + dtdx*dqdt;
      dqdy[n+e*mesh->Np] = drdy*dqdr + dsdy*dqds + dtdy*dqdt;
      dqdz[n+e*mesh->Np] = drdz*dqdr + dsdz*dqds + dtdz*dqdt;
    }
  }
}
