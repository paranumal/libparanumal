#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshPhysicalNodesTriP2D(mesh2D *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */
    iint N = mesh->N[e];

    iint id = e*mesh->Nverts+0;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    
    for(iint n=0;n<mesh->Np[N];++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[N][n]; 
      dfloat sn = mesh->s[N][n];

      /* physical coordinate of interpolation node */
      mesh->x[e*mesh->NpMax+n] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      mesh->y[e*mesh->NpMax+n] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
    }
  }
}
