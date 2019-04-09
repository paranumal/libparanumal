#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

// custom geometric factors specialized for 3D quad on sphere

void meshGeometricFactorsQuad3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 20;
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    for(int j=0;j<mesh->Nq;++j){
      for(int i=0;i<mesh->Nq;++i){

	dfloat xij = mesh->x[i+j*mesh->Nq+e*mesh->Np];
	dfloat yij = mesh->y[i+j*mesh->Nq+e*mesh->Np];
	dfloat zij = mesh->z[i+j*mesh->Nq+e*mesh->Np];

	dfloat xr = 0, yr = 0, zr = 0;
	dfloat xs = 0, ys = 0, zs = 0;
	dfloat xrm = 0, yrm = 0, zrm = 0;
	dfloat xsm = 0, ysm = 0, zsm = 0;
	
	for(int n=0;n<mesh->Nq;++n){

	  dfloat Din = mesh->D[i*mesh->Nq+n];
	  dfloat Djn = mesh->D[j*mesh->Nq+n];

	  xr += Din*mesh->x[n+j*mesh->Nq+e*mesh->Np];
	  yr += Din*mesh->y[n+j*mesh->Nq+e*mesh->Np];
	  zr += Din*mesh->z[n+j*mesh->Nq+e*mesh->Np];

	  xs += Djn*mesh->x[i+n*mesh->Nq+e*mesh->Np];
	  ys += Djn*mesh->y[i+n*mesh->Nq+e*mesh->Np];
	  zs += Djn*mesh->z[i+n*mesh->Nq+e*mesh->Np];

	}

	for(int n=0;n<mesh->Nq;++n){

	  dfloat Din = mesh->weakD[i*mesh->Nq+n];
	  dfloat Djn = mesh->weakD[j*mesh->Nq+n];

	  xrm += Din*mesh->x[n+j*mesh->Nq+e*mesh->Np];
	  yrm += Din*mesh->y[n+j*mesh->Nq+e*mesh->Np];
	  zrm += Din*mesh->z[n+j*mesh->Nq+e*mesh->Np];

	  xsm += Djn*mesh->x[i+n*mesh->Nq+e*mesh->Np];
	  ysm += Djn*mesh->y[i+n*mesh->Nq+e*mesh->Np];
	  zsm += Djn*mesh->z[i+n*mesh->Nq+e*mesh->Np];

	}

	dfloat rx = ys*zij - zs*yij; // dXds x X
	dfloat ry = zs*xij - xs*zij;
	dfloat rz = xs*yij - ys*xij;

	dfloat sx = zr*yij - yr*zij; // -dXdr x X
	dfloat sy = xr*zij - zr*xij;
	dfloat sz = yr*xij - xr*yij;

	dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
	dfloat ty = zr*xs - xr*zs;
	dfloat tz = xr*ys - yr*xs;

	dfloat rxm = ysm*zij - zsm*yij; // dXds x X
	dfloat rym = zsm*xij - xsm*zij;
	dfloat rzm = xsm*yij - ysm*xij;

	dfloat sxm = zrm*yij - yrm*zij; // -dXdr x X
	dfloat sym = xrm*zij - zrm*xij;
	dfloat szm = yrm*xij - xrm*yij;

	dfloat txm = yrm*zsm - zrm*ysm; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
	dfloat tym = zrm*xsm - xrm*zsm;
	dfloat tzm = xrm*ysm - yrm*xsm;

	dfloat Gx = tx, Gy = ty, Gz = tz;

#if 0
	dfloat foo = xij*tx+yij*ty+zij*tz;
	printf("foo = %g [%g,%g,%g,%g,%g,%g] %d\n", foo,
	       yr,mesh->x[3+j*mesh->Nq+e*mesh->Np],mesh->D[i*mesh->Nq+3],zr,mesh->D[j*mesh->Nq + 3],mesh->x[i+3*mesh->Nq+e*mesh->Np],mesh->cubeFaceNumber[e]);
#endif
	
	dfloat J = xij*tx + yij*ty + zij*tz;

	if(J<1e-8) { printf("Negative or small Jacobian: %g  %d %g %g %g %g %g %g\n", J,mesh->cubeFaceNumber[e],xr,yr,zr,xs,ys,zs); /*exit(-1);*/}
	
	rx /= J;
	ry /= J;
	rz /= J;

	sx /= J;
	sy /= J;
	sz /= J;

	tx /= J;
	ty /= J;
	tz /= J;

	rxm /= J;
	rym /= J;
	rzm /= J;

	sxm /= J;
	sym /= J;
	szm /= J;

	txm /= J;
	tym /= J;
	tzm /= J;
	
	// use this for "volume" Jacobian
	J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
	if(J<1e-8) { printf("Negative or small Jacobian: %g\n", J); exit(-1);}

	dfloat JW = mesh->gllw[i]*mesh->gllw[j]*J;
	
	/* store geometric factors */
	int base = mesh->Nvgeo*mesh->Np*e + j*mesh->Nq + i;

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
	mesh->vgeo[base + mesh->Np*JWID] = JW;
	mesh->vgeo[base + mesh->Np*RXMID] = rxm;
	mesh->vgeo[base + mesh->Np*RYMID] = rym;
	mesh->vgeo[base + mesh->Np*RZMID] = rzm;
	mesh->vgeo[base + mesh->Np*SXMID] = sxm;
	mesh->vgeo[base + mesh->Np*SYMID] = sym;
	mesh->vgeo[base + mesh->Np*SZMID] = szm;
	mesh->vgeo[base + mesh->Np*TXMID] = txm;
	mesh->vgeo[base + mesh->Np*TYMID] = tym;
	mesh->vgeo[base + mesh->Np*TZMID] = tzm;
      }
    }
  }
}
