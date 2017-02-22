#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void meshPhysicalNodesTetP3D(mesh3D *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */
    iint N = mesh->N[e];

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    
    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];
    
    for(iint n=0;n<mesh->Np[N];++n){ /* for each node */
      
      /* (r,s,t) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[N][n]; 
      dfloat sn = mesh->s[N][n];
      dfloat tn = mesh->t[N][n];

      /* physical coordinate of interpolation node */
      mesh->x[e*mesh->NpMax+n] = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
      mesh->y[e*mesh->NpMax+n] = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
      mesh->z[e*mesh->NpMax+n] = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
    }
  }
}
