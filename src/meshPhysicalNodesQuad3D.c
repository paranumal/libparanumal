#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void meshPhysicalNodesQuad3D(mesh_t *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  mesh->z = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){ /* for each element */

    int id = e*mesh->Nverts;

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

    
    for(int n=0;n<mesh->Np;++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n]; 
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      dfloat xlin = 
	+0.25*(1-rn)*(1-sn)*xe1
	+0.25*(1+rn)*(1-sn)*xe2
	+0.25*(1+rn)*(1+sn)*xe3
	+0.25*(1-rn)*(1+sn)*xe4;

      dfloat ylin =
	+0.25*(1-rn)*(1-sn)*ye1
	+0.25*(1+rn)*(1-sn)*ye2
	+0.25*(1+rn)*(1+sn)*ye3
	+0.25*(1-rn)*(1+sn)*ye4;

      dfloat zlin =
	+0.25*(1-rn)*(1-sn)*ze1
	+0.25*(1+rn)*(1-sn)*ze2
	+0.25*(1+rn)*(1+sn)*ze3
	+0.25*(1-rn)*(1+sn)*ze4;

      //      printf("xlin=%g, ylin=%g, zlin=%g\n", xlin, ylin, zlin);
      
      // project to sphere
      dfloat rlin = sqrt(xlin*xlin+ylin*ylin+zlin*zlin);
      mesh->x[cnt] = xlin/rlin; 
      mesh->y[cnt] = ylin/rlin; 
      mesh->z[cnt] = zlin/rlin; 
	
      ++cnt;
    }
  }
}
