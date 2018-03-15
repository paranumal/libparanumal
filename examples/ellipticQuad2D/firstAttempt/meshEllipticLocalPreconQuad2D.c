#include <stdlib.h>
#include "mesh2D.h"

void meshEllipticLocalPreconQuad2D(mesh2D *mesh, dfloat *qL, dfloat lambda, dfloat *PqL){
  // fetch element plus some halo stuff
  dfloat *qLe  = (dfloat*) calloc((mesh->Nq+2)*(mesh->Nq+2), sizeof(dfloat));
  dfloat *D    = (dfloat*) calloc(mesh->Nq*mesh->Nq, sizeof(dfloat));
  dfloat *tmpr = (dfloat*) calloc(mesh->Nq*mesh->Nq, sizeof(dfloat));
  dfloat *tmps = (dfloat*) calloc(mesh->Nq*mesh->Nq, sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){

    // prefetch qLe
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){
	qLe[i+j*mesh->Nq] = qL[i+j*mesh->Nq+e*mesh->Np];
	 D[i+j*mesh->Nq] = mesh->D[i+j*mesh->Nq];
      }
    }

    // local 'r' and 's' derivatives
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){
	
	dfloat dqdrL = 0, dqdsL = 0;

	for(iint n=0;n<mesh->Nq;++n){
	  dqdrL += D[i*mesh->Nq+n]*qLe[n+j*mesh->Nq];
	  dqdsL += D[j*mesh->Nq+n]*qLe[i+n*mesh->Nq];
	}

	// load geometric factors
	iint gid = mesh->Np*mesh->Nvgeo*e + (i+mesh->Nq*j);

	dfloat drdx = mesh->vgeo[gid + mesh->Np*RXID];
	dfloat drdy = mesh->vgeo[gid + mesh->Np*RYID];
	dfloat dsdx = mesh->vgeo[gid + mesh->Np*SXID];
	dfloat dsdy = mesh->vgeo[gid + mesh->Np*SYID];
	dfloat Jw   = mesh->vgeo[gid + mesh->Np*JWID];
	
	dfloat dqdxL = drdx*dqdrL+dsdx*dqdsL;
	dfloat dqdyL = drdy*dqdrL+dsdy*dqdsL;
	
	tmpr[j*mesh->Nq+i] = Jw*(drdx*dqdxL + drdy*dqdyL);
	tmps[j*mesh->Nq+i] = Jw*(dsdx*dqdxL + dsdy*dqdyL);
      }
    }

    // local op
    for(iint j=0;j<mesh->Nq;++j){
      for(iint i=0;i<mesh->Nq;++i){
	iint lid = i + mesh->Nq*j;
	iint gid = mesh->Np*mesh->Nvgeo*e + lid + mesh->Np*JWID;
	dfloat Jw = mesh->vgeo[gid];
	dfloat Aq = lambda*Jw*qLe[lid]; // lambda*M*q
	
	for(iint n=0;n<mesh->Nq;++n){
	  Aq += D[n*mesh->Nq+i]*tmpr[j*mesh->Nq+n];
	  Aq += D[n*mesh->Nq+j]*tmps[n*mesh->Nq+i];
	}
	
	AqL[lid+mesh->Np*e] = Aq;
      }
    }
  }
  
  free(qLe); free(D); free(tmpr); free(tmps);
}
