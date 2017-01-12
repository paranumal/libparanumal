#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshPhysicalNodesHex3D(mesh3D *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];
    dfloat xe5 = mesh->EX[id+4]; 
    dfloat xe6 = mesh->EX[id+5];
    dfloat xe7 = mesh->EX[id+6];
    dfloat xe8 = mesh->EX[id+7];
    
    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    dfloat ye5 = mesh->EY[id+4]; 
    dfloat ye6 = mesh->EY[id+5];
    dfloat ye7 = mesh->EY[id+6];
    dfloat ye8 = mesh->EY[id+7];

    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];
    dfloat ze5 = mesh->EZ[id+4]; 
    dfloat ze6 = mesh->EZ[id+5];
    dfloat ze7 = mesh->EZ[id+6];
    dfloat ze8 = mesh->EZ[id+7];

    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s,t) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];
      dfloat tn = mesh->t[n];

      /* physical coordinate of interpolation node */
      mesh->x[cnt] = 
	+0.125*(1-rn)*(1-sn)*(1-tn)*xe1
	+0.125*(1+rn)*(1-sn)*(1-tn)*xe2
	+0.125*(1+rn)*(1+sn)*(1-tn)*xe3
	+0.125*(1-rn)*(1+sn)*(1-tn)*xe4
	+0.125*(1-rn)*(1-sn)*(1+tn)*xe5
	+0.125*(1+rn)*(1-sn)*(1+tn)*xe6
	+0.125*(1+rn)*(1+sn)*(1+tn)*xe7
	+0.125*(1-rn)*(1+sn)*(1+tn)*xe8;

      mesh->y[cnt] = 
	+0.125*(1-rn)*(1-sn)*(1-tn)*ye1
	+0.125*(1+rn)*(1-sn)*(1-tn)*ye2
	+0.125*(1+rn)*(1+sn)*(1-tn)*ye3
	+0.125*(1-rn)*(1+sn)*(1-tn)*ye4
	+0.125*(1-rn)*(1-sn)*(1+tn)*ye5
	+0.125*(1+rn)*(1-sn)*(1+tn)*ye6
	+0.125*(1+rn)*(1+sn)*(1+tn)*ye7
	+0.125*(1-rn)*(1+sn)*(1+tn)*ye8;

      mesh->z[cnt] = 
	+0.125*(1-rn)*(1-sn)*(1-tn)*ze1
	+0.125*(1+rn)*(1-sn)*(1-tn)*ze2
	+0.125*(1+rn)*(1+sn)*(1-tn)*ze3
	+0.125*(1-rn)*(1+sn)*(1-tn)*ze4
	+0.125*(1-rn)*(1-sn)*(1+tn)*ze5
	+0.125*(1+rn)*(1-sn)*(1+tn)*ze6
	+0.125*(1+rn)*(1+sn)*(1+tn)*ze7
	+0.125*(1-rn)*(1+sn)*(1+tn)*ze8;

      ++cnt;
    }
  }
}
