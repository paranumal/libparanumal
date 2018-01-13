#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsTri3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 7;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				mesh->Nsgeo*mesh->Nfp*mesh->Nfaces, 
				sizeof(dfloat));

  dfloat *xr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *yr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zr = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *xs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *ys = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat *zs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat *J  = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

    for(int n=0;n<mesh->Np;++n){
      
      dfloat x = mesh->x[n+e*mesh->Np];
      dfloat y = mesh->y[n+e*mesh->Np];
      dfloat z = mesh->z[n+e*mesh->Np];
      
      dfloat xrn = 0, yrn = 0, zrn = 0;
      dfloat xsn = 0, ysn = 0, zsn = 0;
      
      for(int m=0;m<mesh->Nq;++m){
	
	dfloat Dr = mesh->Dr[n*mesh->Np+m];
	dfloat Ds = mesh->Ds[n*mesh->Np+m];
	
	xrn += Dr*mesh->x[m+e*mesh->Np];
	yrn += Dr*mesh->y[m+e*mesh->Np];
	zrn += Dr*mesh->z[m+e*mesh->Np];
	
	xsn += Ds*mesh->x[m+e*mesh->Np];
	ysn += Ds*mesh->y[m+e*mesh->Np];
	zsn += Ds*mesh->z[m+e*mesh->Np];

      }
      
      dfloat txn = yrn*zsn - zrn*ysn;
      dfloat tyn = zrn*xsn - xrn*zsn;
      dfloat tzn = xrn*ysn - yrn*xsn;
      
      dfloat Gx = txn, Gy = tyn, Gz = tzn;
      
      dfloat Jn = x*txn + y*tyn + z*tzn;
      
      xr[n] = xrn;
      yr[n] = yrn;
      zr[n] = zrn;
      
      xs[n] = xsn;
      ys[n] = ysn;
      zs[n] = zsn;
      
      J[n] = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
    }
    
    
    // face 0
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int id = mesh->faceNodes[n+f*mesh->Nq];
	
	dfloat xid = mesh->x[id+e*mesh->Np];
	dfloat yid = mesh->y[id+e*mesh->Np];
	dfloat zid = mesh->z[id+e*mesh->Np];
	dfloat Jid = J[id];
	
	dfloat nx, ny, nz;
	
	if(f==0){
	  nx = yr[id]*zid - zr[id]*yid;
	  ny = zr[id]*xid - xr[id]*zid;
	  nz = xr[id]*yid - yr[id]*xid;
	}
	
	if(f==1){
	  nx = (ys[id]-yr[id])*zid - (zs[id]-zr[id])*yid;
	  ny = (zs[id]-zr[id])*xid - (xs[id]-xr[id])*zid;
	  nz = (xs[id]-xr[id])*yid - (ys[id]-yr[id])*xid;
	}
	
	if(f==2){
	  nx = -ys[id]*zid + zs[id]*yid;
	  ny = -zs[id]*xid + xs[id]*zid;
	  nz = -xs[id]*yid + ys[id]*xid;
	}
	
	dfloat R = sqrt(xid*xid+yid*yid+zid*zid);
	
	nx /= R;
	ny /= R;
	nz /= R;
	
	dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
	
	nx /= sJ;
	ny /= sJ;
	nz /= sJ;
	
	if(sJ<1e-8) { printf("Negative or small surface Jacobian: %g\n", sJ); exit(-1);}
	
	int base = e*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + n + f*mesh->Nq;
	
	mesh->sgeo[base+mesh->Nq*mesh->Nfaces*NXID] = nx;
	mesh->sgeo[base+mesh->Nq*mesh->Nfaces*NYID] = ny;
	mesh->sgeo[base+mesh->Nq*mesh->Nfaces*NZID] = nz;
	mesh->sgeo[base+mesh->Nq*mesh->Nfaces*SJID] = sJ;
	
	mesh->sgeo[base+mesh->Nq*mesh->Nfaces*IJID] = 1./Jid;
      }
    }
  }


#if 0
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int idM = n+f*mesh->Nq+e*mesh->Nfaces*mesh->Nq;
	int idP = mesh->mapP[idM];
	int eP = idP/(mesh->Nq*mesh->Nfaces);
	int fP = (idP%(mesh->Nq*mesh->Nfaces))/mesh->Nq;
	int nP = (idP%mesh->Nq);
	int baseM = e*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + f*mesh->Nq*mesh->Nsgeo + n;
	int baseP = eP*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + fP*mesh->Nq*mesh->Nsgeo + nP;
	printf("e,f,n=(%d,%d,%d)-(%d,%d,%d): xP-xM=(%g,%g,%g) : norP+norM=%g,%g,%g\n",
	       e,f,n,eP,fP,nP,
	       mesh->x[mesh->vmapP[idM]]-mesh->x[mesh->vmapM[idM]],
	       mesh->y[mesh->vmapP[idM]]-mesh->y[mesh->vmapM[idM]],
	       mesh->z[mesh->vmapP[idM]]-mesh->z[mesh->vmapM[idM]],
	       mesh->sgeo[baseM+NXID*mesh->Nq]+mesh->sgeo[baseP+NXID*mesh->Nq],
	       mesh->sgeo[baseM+NYID*mesh->Nq]+mesh->sgeo[baseP+NYID*mesh->Nq],
	       mesh->sgeo[baseM+NZID*mesh->Nq]+mesh->sgeo[baseP+NZID*mesh->Nq]);

      }
    }
  }
#endif  
  // TW: omit 1/min(h) calculation
}
