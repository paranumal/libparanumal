
#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshSurfaceGeometricFactorsTri2D(mesh2D *mesh){

  /* unified storage array for geometric factors */
  mesh->Nsgeo = 6;
  mesh->sgeo = (dfloat*) calloc((mesh->Nelements+mesh->totalHaloPairs)*
				mesh->Nsgeo*mesh->Nfaces, 
				sizeof(dfloat));
  
  for(iint e=0;e<mesh->Nelements+mesh->totalHaloPairs;++e){ /* for each element */

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
    mesh->sgeo[base+SJID] = d2/2.; // TW fixed bug d1=>d2
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

  
  for(iint e=0;e<mesh->Nelements;++e){ /* for each non-halo element */
    for(iint f=0;f<mesh->Nfaces;++f){
      iint baseM = e*mesh->Nfaces + f;

      // awkward: (need to find eP,fP relative to bulk+halo)
      iint idP = mesh->vmapP[e*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+0];
      iint eP = (idP>=0) ? (idP/mesh->Np):e;

      iint fP = mesh->EToF[baseM];
      fP = (fP==-1)?f:fP;

      iint baseP = eP*mesh->Nfaces + fP;
      
      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)  A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
      dfloat hinvM = 0.5*mesh->sgeo[baseM*mesh->Nsgeo + SJID]*mesh->sgeo[baseM*mesh->Nsgeo + IJID];
      dfloat hinvP = 0.5*mesh->sgeo[baseP*mesh->Nsgeo + SJID]*mesh->sgeo[baseP*mesh->Nsgeo + IJID];
      
      mesh->sgeo[baseM*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
      mesh->sgeo[baseP*mesh->Nsgeo+IHID] = mymax(hinvM,hinvP);
#if 0
      printf("e=%d f=%d (eP=%d,fP=%d) nx=%5.4f, ny=%5.4f, sJ=%5.4f, invJ=%5.4f, hinv=%f\n"
	     ,e,f,eP,fP,
	     mesh->sgeo[baseM*mesh->Nsgeo+NXID],
	     mesh->sgeo[baseM*mesh->Nsgeo+NYID],
	     mesh->sgeo[baseM*mesh->Nsgeo+SJID],
	     mesh->sgeo[baseM*mesh->Nsgeo+IJID],
	     mesh->sgeo[baseM*mesh->Nsgeo+IHID]);
#endif
    }
  }

}
