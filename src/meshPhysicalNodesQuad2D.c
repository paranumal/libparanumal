#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshPhysicalNodesQuad2D(mesh2D *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    
    for(iint n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      mesh->x[cnt] = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      mesh->y[cnt] = 
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      ++cnt;
    }
  }
}
