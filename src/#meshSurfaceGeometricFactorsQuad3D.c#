#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsQuad3D(mesh_t *mesh){

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

    for(int j=0;j<mesh->Nq;++j){
      for(int i=0;i<mesh->Nq;++i){

	dfloat x = mesh->x[i+j*mesh->Nq+e*mesh->Np];
	dfloat y = mesh->y[i+j*mesh->Nq+e*mesh->Np];
	dfloat z = mesh->z[i+j*mesh->Nq+e*mesh->Np];
	
	dfloat xrij = 0, yrij = 0, zrij = 0;
	dfloat xsij = 0, ysij = 0, zsij = 0;

	for(int n=0;n<mesh->Nq;++n){

	  dfloat Din = mesh->D[i*mesh->Nq+n];
	  dfloat Djn = mesh->D[j*mesh->Nq+n];

	  xrij += Din*mesh->x[n+j*mesh->Nq+e*mesh->Np];
	  yrij += Din*mesh->y[n+j*mesh->Nq+e*mesh->Np];
	  zrij += Din*mesh->z[n+j*mesh->Nq+e*mesh->Np];

	  xsij += Djn*mesh->x[i+n*mesh->Nq+e*mesh->Np];
	  ysij += Djn*mesh->y[i+n*mesh->Nq+e*mesh->Np];
	  zsij += Djn*mesh->z[i+n*mesh->Nq+e*mesh->Np];

	}

	dfloat txij = yrij*zsij - zrij*ysij;
	dfloat tyij = zrij*xsij - xrij*zsij;
	dfloat tzij = xrij*ysij - yrij*xsij;

	dfloat Gx = txij, Gy = tyij, Gz = tzij;
	
	dfloat Jij = x*txij + y*tyij + z*tzij;
	
	xr[i+j*mesh->Nq] = xrij;
	yr[i+j*mesh->Nq] = yrij;
	zr[i+j*mesh->Nq] = zrij;
			       
	xs[i+j*mesh->Nq] = xsij;
	ys[i+j*mesh->Nq] = ysij;
	zs[i+j*mesh->Nq] = zsij;

	J[i+j*mesh->Nq] = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
      }
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
	  nx = ys[id]*zid - zs[id]*yid;
	  ny = zs[id]*xid - xs[id]*zid;
	  nz = xs[id]*yid - ys[id]*xid;
	}

	if(f==2){
	  nx = -yr[id]*zid + zr[id]*yid;
	  ny = -zr[id]*xid + xr[id]*zid;
	  nz = -xr[id]*yid + yr[id]*xid;
	}

	if(f==3){
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
	
	int base = e*mesh->Nq*mesh->Nfaces*mesh->Nsgeo + n + f*mesh->Nq*mesh->Nsgeo;
	
	mesh->sgeo[base+mesh->Nq*NXID] = nx;
	mesh->sgeo[base+mesh->Nq*NYID] = ny;
	mesh->sgeo[base+mesh->Nq*NZID] = nz;
	mesh->sgeo[base+mesh->Nq*SJID] = sJ;

	mesh->sgeo[base+mesh->Nq*IJID] = 1./Jid;
	
	mesh->sgeo[base+mesh->Nq*WSJID] = sJ*mesh->gllw[n];
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
