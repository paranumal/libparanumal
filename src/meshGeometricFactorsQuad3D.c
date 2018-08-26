/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

// custom geometric factors specialized for 3D quad on sphere

void meshGeometricFactorsQuad3D(mesh_t *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 11; // 
  
  /* note that we have volume geometric factors for each node */
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo*mesh->Np, sizeof(dfloat));

  for(int e=0;e<mesh->Nelements;++e){ /* for each element */

    for(int j=0;j<mesh->Nq;++j){
      for(int i=0;i<mesh->Nq;++i){

	dfloat xij = mesh->x[i+j*mesh->Nq+e*mesh->Np];
	dfloat yij = mesh->y[i+j*mesh->Nq+e*mesh->Np];
	dfloat zij = mesh->z[i+j*mesh->Nq+e*mesh->Np];

	dfloat xr = 0, yr = 0, zr = 0;
	dfloat xs = 0, ys = 0, zs = 0;
	
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

	dfloat rx = ys*zij - zs*yij; // dXds x X
	dfloat ry = zs*xij - xs*zij;
	dfloat rz = xs*yij - ys*xij;

	dfloat sx = zr*yij - yr*zij; // -dXdr x X
	dfloat sy = xr*zij - zr*xij;
	dfloat sz = yr*xij - xr*yij;

	dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
	dfloat ty = zr*xs - xr*zs;
	dfloat tz = xr*ys - yr*xs;

	dfloat Gx = tx, Gy = ty, Gz = tz;

#if 0
	dfloat foo = xij*tx+yij*ty+zij*tz;
	printf("foo = %g [%g,%g,%g] [%g,%g,%g]\n", foo,
	       xij,yij,zij,tx-foo*xij,ty-foo*yij,tz-foo*zij);
#endif
	
	dfloat J = xij*tx + yij*ty + zij*tz;

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

	dfloat JW = J*mesh->gllw[i]*mesh->gllw[j];
	
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
      
      }
    }
  }
}
