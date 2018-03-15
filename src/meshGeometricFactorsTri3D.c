#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

// custom geometric factors specialized for 3D tri on sphere

void meshGeometricFactorsTri3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 11; // 
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    for(int n=0;n<mesh->Np;++n){
      
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      dfloat zn = mesh->z[n+e*mesh->Np];
      
      dfloat xr = 0, yr = 0, zr = 0;
      dfloat xs = 0, ys = 0, zs = 0;
      
      for(int m=0;m<mesh->Np;++m){
	
	dfloat Drnm = mesh->Dr[n*mesh->Np+m];
	dfloat Dsnm = mesh->Ds[n*mesh->Np+m];
	
	xr += Drnm*mesh->x[m+e*mesh->Np];
	yr += Drnm*mesh->y[m+e*mesh->Np];
	zr += Drnm*mesh->z[m+e*mesh->Np];
	
	xs += Dsnm*mesh->x[m+e*mesh->Np];
	ys += Dsnm*mesh->y[m+e*mesh->Np];
	zs += Dsnm*mesh->z[m+e*mesh->Np];
	
      }
      
      dfloat rx = ys*zn - zs*yn; // dXds x X
      dfloat ry = zs*xn - xs*zn;
      dfloat rz = xs*yn - ys*xn;
      
      dfloat sx = zr*yn - yr*zn; // -dXdr x X
      dfloat sy = xr*zn - zr*xn;
      dfloat sz = yr*xn - xr*yn;
      
      dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
      dfloat ty = zr*xs - xr*zs;
      dfloat tz = xr*ys - yr*xs;
      
      dfloat Gx = tx, Gy = ty, Gz = tz;
      
      dfloat J = xn*tx + yn*ty + zn*tz;
      
      if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
      
      rx /= J;
      ry /= J;
      rz /= J;
      
      sx /= J;
      sy /= J;
      sz /= J;
      
      tx /= J;
      ty /= J;
      tz /= J;
      
      // use this for "volume" Jacobian
      J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
      
      if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}
      
      /* store geometric factors */
      int base = mesh->Nvgeo*mesh->Np*e + n;

      mesh->vgeo[base + mesh->Np*RXID] = rx;
      mesh->vgeo[base + mesh->Np*RYID] = ry;
      mesh->vgeo[base + mesh->Np*RZID] = rz;
      mesh->vgeo[base + mesh->Np*SXID] = sx;
      mesh->vgeo[base + mesh->Np*SYID] = sy;
      mesh->vgeo[base + mesh->Np*SZID] = sz;
      mesh->vgeo[base + mesh->Np*TXID] = tx;
      mesh->vgeo[base + mesh->Np*TYID] = ty;
      mesh->vgeo[base + mesh->Np*TZID] = tz;
      mesh->vgeo[base + mesh->Np*JID]  = J;
      
    }
  }
}
