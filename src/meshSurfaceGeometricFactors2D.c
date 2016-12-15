
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshSurfaceGeometricFactors2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 4;
  mesh->sgeo = (dfloat*) calloc(mesh->Nelements*mesh->Nsgeo*mesh->Nfaces, 
				sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    int id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0];
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0];
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));
    if(J<0) printf("bugger: got negative geofac\n");
    
    /* face 1 */
    int base = mesh->Nsgeo*mesh->Nfaces*e;
    dfloat nx1 = ye2-ye1;
    dfloat ny1 = -(xe2-xe1);
    dfloat  d1 = norm(nx1,ny1);

    mesh->sgeo[base+NXID] = nx1/d1;
    mesh->sgeo[base+NYID] = ny1/d1;
    mesh->sgeo[base+SJID] = d1/2.;
    mesh->sgeo[base+IJID] = 1./J;

    /* face 2 */
    base += mesh->Nsgeo;
    dfloat nx2 = ye3-ye2;
    dfloat ny2 = -(xe3-xe2);
    dfloat  d2 = norm(nx2,ny2);

    mesh->sgeo[base+NXID] = nx2/d2;
    mesh->sgeo[base+NYID] = ny2/d2;
    mesh->sgeo[base+SJID] = d1/2.;
    mesh->sgeo[base+IJID] = 1./J;

    /* face 3 */
    base += mesh->Nsgeo;
    dfloat nx3 = ye1-ye3;
    dfloat ny3 = -(xe1-xe3);
    dfloat  d3 = norm(nx3,ny3);

    mesh->sgeo[base+NXID] = nx3/d3;
    mesh->sgeo[base+NYID] = ny3/d3;
    mesh->sgeo[base+SJID] = d3/2.;
    mesh->sgeo[base+IJID] = 1./J;
  }
}
