
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshGeometricFactors2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nvgeo = 5;
  mesh->vgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nvgeo, 
				sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    int id = e*mesh->Nverts+0;

    dfloat xe1 = mesh->EX[id+0];
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0];
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));

    if(J<0) printf("bugger: got negative geofac\n");
    dfloat rx =  (0.5/J)*(ye3-ye1);
    dfloat ry = -(0.5/J)*(xe3-xe1);
    dfloat sx = -(0.5/J)*(ye2-ye1);
    dfloat sy =  (0.5/J)*(xe2-xe1);
    
    /* store geometric factors */
    mesh->vgeo[mesh->Nvgeo*e + RXID] = rx;
    mesh->vgeo[mesh->Nvgeo*e + RYID] = ry;
    mesh->vgeo[mesh->Nvgeo*e + SXID] = sx;
    mesh->vgeo[mesh->Nvgeo*e + SYID] = sy;
    mesh->vgeo[mesh->Nvgeo*e +  JID] = J;
  }
}
